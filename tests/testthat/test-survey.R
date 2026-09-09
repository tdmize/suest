survey_fixture <- function(fpc = FALSE) {
  skip_if_not_installed("survey")
  set.seed(8041)
  d <- data.frame(id = 1:120, strata = rep(1:4, each = 30),
    psu = rep(1:20, each = 6), unit = rep(1:6, 20))
  d$x <- rnorm(120) + rep(rnorm(20), each = 6)
  d$z <- rnorm(120)
  d$w <- 1 + (d$id %% 9) / 4
  d$y <- 2 + .7 * d$x - .4 * d$z + rep(rnorm(20), each = 6) + rnorm(120)
  d$y2 <- 1 - .2 * d$x + .8 * d$z + .6 * d$y + rnorm(120)
  d$pop <- rep(c(10, 15, 20, 25), each = 30)
  design <- if (fpc) survey::svydesign(~psu, strata = ~strata, weights = ~w, fpc = ~pop, data = d) else
    survey::svydesign(~psu, strata = ~strata, weights = ~w, data = d)
  list(data = d, design = design)
}

# Independent reference: explicit PSU totals and stratum centering, with no
# survey variance call. Uses native retained influences, unlike the adapter.
survey_manual_joint <- function(models, design) {
  influence <- lapply(models, function(m) {
    z <- matrix(0, nrow(design$variables), length(coef(m)))
    rows <- match(m$survey.design$variables$id, design$variables$id)
    z[rows, ] <- attr(m, "influence")
    z
  })
  z <- do.call(cbind, influence)
  V <- matrix(0, ncol(z), ncol(z))
  for (h in unique(design$variables$strata)) {
    keep <- design$variables$strata == h
    totals <- rowsum(z[keep, , drop = FALSE], design$variables$psu[keep])
    G <- nrow(totals)
    if (G == 1) next
    totals <- sweep(totals, 2, colMeans(totals))
    correction <- if (is.null(design$fpc$popsize)) 1 else 1 - G / design$fpc$popsize[which(keep)[1], 1]
    V <- V + crossprod(totals) * G / (G - 1) * correction
  }
  V
}

test_that("survey covariance reproduces native blocks and independent PSU algebra", {
  withr::local_options(survey.lonely.psu = "fail", survey.adjust.domain.lonely = FALSE)
  for (fpc in c(FALSE, TRUE)) {
    d <- survey_fixture(fpc)$design
    a <- survey::svyglm(y ~ x, d, influence = TRUE)
    b <- survey::svyglm(y2 ~ x + z, d, influence = TRUE, rescale = FALSE)
    s <- suest(a, b, survey_design = d, observation_id = "id")
    expect_equal(unname(vcov(s)), survey_manual_joint(list(a, b), d), tolerance = 1e-10)
    expect_equal(unname(vcov(s)[s$index[[1]], s$index[[1]]]), unname(vcov(a)), tolerance = 1e-10)
    expect_equal(unname(vcov(s)[s$index[[2]], s$index[[2]]]), unname(vcov(b)), tolerance = 1e-10)
    expect_equal(unname(coef(s)), c(unname(coef(a)), unname(coef(b))))
    expect_equal(s$survey$design_df, 16L)
    expect_false(any(grepl("lnvar", names(coef(s)))))
  }
})

test_that("survey overlap, domains, omitted PSUs, and three models align", {
  withr::local_options(survey.lonely.psu = "fail", survey.adjust.domain.lonely = FALSE)
  d <- survey_fixture(TRUE)$design
  subsets <- list(list(quote(unit <= 4), quote(unit >= 3)),
    list(quote(unit <= 3), quote(unit > 3)),
    list(quote(psu %% 5 == 1), quote(psu %% 5 == 2)),
    list(quote(strata <= 2), quote(strata > 2)))
  for (sub in subsets) {
    a <- survey::svyglm(y ~ x, subset(d, eval(sub[[1]])), influence = TRUE)
    b <- survey::svyglm(y2 ~ x + z, subset(d, eval(sub[[2]])), influence = TRUE)
    s <- suest(a, b, survey_design = d, observation_id = "id")
    expect_equal(unname(vcov(s)), survey_manual_joint(list(a, b), d), tolerance = 1e-10)
    expect_equal(s$nobs_overlap, length(intersect(a$survey.design$variables$id, b$survey.design$variables$id)))
    cross <- vcov(s)[s$index[[1]], s$index[[2]]]
    if (identical(sub, subsets[[4]])) expect_equal(unname(cross), matrix(0, 2, 3)) else
      expect_gt(max(abs(cross)), 1e-5)
  }
  a <- survey::svyglm(y ~ x, d, influence = TRUE)
  s <- suest(a, a, a, survey_design = d, observation_id = "id")
  expect_equal(unname(vcov(s)), survey_manual_joint(list(a, a, a), d), tolerance = 1e-10)
  expect_equal(unname(vcov(s)[s$index[[1]], s$index[[2]]]), unname(vcov(a)))
})

test_that("missingness and row reordering preserve survey alignment", {
  withr::local_options(survey.lonely.psu = "fail", survey.adjust.domain.lonely = FALSE)
  d <- survey_fixture()$design
  d$variables$y[d$variables$id %% 7 == 0] <- NA
  a <- survey::svyglm(y ~ x, d, na.action = na.exclude)
  b <- survey::svyglm(y2 ~ x + z, d)
  s <- suest(a, b, survey_design = d, observation_id = "id")
  expect_equal(unname(vcov(s)[s$index[[1]], s$index[[1]]]), unname(vcov(a)), tolerance = 1e-10)
  revd <- d[rev(seq_len(nrow(d$variables))), ]
  revb <- survey::svyglm(y2 ~ x + z, revd)
  r <- suest(a, revb, survey_design = revd, observation_id = "id")
  expect_equal(unname(vcov(s)), unname(vcov(r)), tolerance = 1e-10)
  expect_equal(nrow(suest_newdata(s)), sum(s$nobs_models))
})

test_that("survey predictions, coefficient perturbation, and effects use joint V", {
  withr::local_options(survey.lonely.psu = "fail", survey.adjust.domain.lonely = FALSE)
  d <- survey_fixture()$design
  a <- survey::svyglm(y ~ x + offset(z), d)
  b <- survey::svyglm(y2 ~ x + z, d)
  s <- suest(a, b, survey_design = d, observation_id = "id")
  nd <- d$variables[1:5, ]
  pred <- marginaleffects::predictions(s, newdata = nd)
  # Native predict.svyglm(newdata=) omits formula offsets in survey 4.5.
  # Check the intended conditional mean directly, including the known offset.
  expected <- c(coef(a)[1] + coef(a)[2] * nd$x + nd$z,
    as.numeric(predict(b, nd, se.fit = FALSE)))
  expect_equal(pred$estimate, expected, tolerance = 1e-8)
  perturbed <- coef(s)
  perturbed[1] <- perturbed[1] + 1
  s2 <- set_coef.suest_model(s, perturbed)
  expect_equal(get_predict.suest_model(s2, nd)$estimate - get_predict.suest_model(s, nd)$estimate, c(rep(1, 5), rep(0, 5)))
  slopes <- marginaleffects::avg_slopes(s, variables = "x", newdata = suest_newdata(s), wts = ".suest_weight")
  ix <- c(s$index[[1]][2], s$index[[2]][2])
  expect_equal(slopes$estimate, unname(coef(s)[ix]), tolerance = 1e-7)
  expect_equal(slopes$std.error, unname(sqrt(diag(vcov(s)))[ix]), tolerance = 2e-5)
  contrast <- marginaleffects::hypotheses(slopes, hypothesis = matrix(c(1, -1), ncol = 1))
  expect_equal(as.numeric(contrast$std.error), sqrt(vcov(s)[ix[1], ix[1]] + vcov(s)[ix[2], ix[2]] - 2 * vcov(s)[ix[1], ix[2]]), tolerance = 2e-5)
  expect_equal(suest_newdata(s)$.suest_weight, rep(d$variables$w, 2))
})

test_that("survey restrictions reject ambiguous or unsupported designs", {
  withr::local_options(survey.lonely.psu = "fail", survey.adjust.domain.lonely = FALSE)
  d <- survey_fixture()$design
  a <- survey::svyglm(y ~ x, d)
  run <- function(design = d, ...) suest(a, a, survey_design = design, observation_id = "id", ...)
  expect_error(suest(a, a), "full common survey design")
  expect_error(suest(a, a, survey_design = d), "observation_id")
  expect_error(run(cluster = "psu"), "without cluster")
  expect_error(run(weight_type = "pweight"), "without cluster")
  expect_error(suest(a, lm(y ~ x, d$variables), survey_design = d, observation_id = "id"), "only survey")
  binary <- survey::svyglm(I(y > 2) ~ x, d, family = quasibinomial())
  expect_error(suest(a, binary, survey_design = d, observation_id = "id"), "does not combine Gaussian and binary")
  wrong <- d
  wrong$prob[1] <- wrong$prob[1] / 2
  expect_error(run(wrong), "must agree")
  wrong <- d
  wrong$variables$id[1] <- 999
  expect_error(run(wrong), "Every fitted observation")
  wrong <- d
  wrong$variables$id[1] <- 2
  expect_error(run(wrong), "unique")
  expect_error(run(subset(d, psu %% 5 != 0)), "all sampled PSUs")
  multi <- survey::svydesign(~psu + id, strata = ~strata, weights = ~w, data = d$variables)
  expect_error(run(multi), "one-stage")
  repd <- survey::as.svrepdesign(d, type = "JKn")
  expect_error(run(repd), "replicate-weight")
  cal <- survey::calibrate(d, ~x, population = c(sum(weights(d)), sum(weights(d) * d$variables$x)))
  expect_error(run(cal), "Calibrated")
  withr::local_options(survey.lonely.psu = "adjust")
  expect_error(run(), "survey.lonely.psu")
})

test_that("survey certainty PSUs and removal policies are explicit", {
  withr::local_options(survey.lonely.psu = "fail", survey.adjust.domain.lonely = FALSE)
  raw <- survey_fixture()$data
  raw$strata[raw$psu == 1] <- 5
  raw$pop <- ifelse(raw$strata == 5, 1, 30)
  d <- survey::svydesign(~psu, strata = ~strata, weights = ~w, fpc = ~pop, data = raw)
  a <- survey::svyglm(y ~ x, d, influence = TRUE)
  s <- suest(a, a, survey_design = d, observation_id = "id")
  expect_equal(unname(vcov(s)), survey_manual_joint(list(a, a), d), tolerance = 1e-10)
  d <- survey::svydesign(~psu, strata = ~strata, weights = ~w, data = raw)
  for (policy in c("remove", "certainty")) {
    options(survey.lonely.psu = policy)
    a <- survey::svyglm(y ~ x, d)
    s <- suest(a, a, survey_design = d, observation_id = "id")
    expect_equal(unname(vcov(s)[s$index[[1]], s$index[[1]]]), unname(vcov(a)), tolerance = 1e-10)
  }
  options(survey.lonely.psu = "fail")
  expect_error(suest(a, a, survey_design = d, observation_id = "id"), "only one PSU")
})

test_that("survey metadata changes and nonstandard fits cannot silently pass", {
  withr::local_options(survey.lonely.psu = "fail", survey.adjust.domain.lonely = FALSE)
  d <- survey_fixture(TRUE)$design
  a <- survey::svyglm(y ~ x, d)
  for (field in c("cluster", "strata")) {
    wrong <- a
    wrong$survey.design[[field]][1, 1] <- 999
    expect_error(suest(a, wrong, survey_design = d, observation_id = "id"), "must agree")
  }
  wrong <- a
  wrong$survey.design$fpc$popsize[, 1] <- wrong$survey.design$fpc$popsize[, 1] * 2
  expect_error(suest(a, wrong, survey_design = d, observation_id = "id"), "must agree")
  wrong <- a
  wrong$cov.unscaled <- wrong$cov.unscaled * 2
  expect_error(suest(a, wrong, survey_design = d, observation_id = "id"), "native design-based covariance")
  wrong <- a
  wrong$coefficients[2] <- NA_real_
  expect_error(suest(a, wrong, survey_design = d, observation_id = "id"), "nonaliased")
  extra <- survey::svyglm(y ~ x, d, weights = w)
  expect_error(suest(a, extra, survey_design = d, observation_id = "id"), "Additional")
  off <- survey::svyglm(y ~ x, d, offset = z)
  expect_error(suest(a, off, survey_design = d, observation_id = "id"), "offsets in the formula")
  withr::local_options(survey.adjust.domain.lonely = TRUE)
  expect_error(suest(a, a, survey_design = d, observation_id = "id"), "survey.adjust.domain.lonely")
})

test_that("survey nested PSU labels, factors, and weight rescaling are stable", {
  withr::local_options(survey.lonely.psu = "fail", survey.adjust.domain.lonely = FALSE)
  raw <- survey_fixture()$data
  raw$psu <- (raw$psu - 1) %% 5 + 1
  raw$f <- factor(raw$unit %% 3)
  d <- survey::svydesign(~psu, strata = ~strata, weights = ~w, nest = TRUE, data = raw)
  a <- survey::svyglm(y ~ x * f, d)
  b <- survey::svyglm(y2 ~ x + f, d)
  s <- suest(a, b, survey_design = d, observation_id = c("strata", "id"))
  nd <- raw[1:6, ]
  expect_equal(get_predict.suest_model(s, nd)$estimate,
    c(as.numeric(predict(a, nd, se.fit = FALSE)), as.numeric(predict(b, nd, se.fit = FALSE))), tolerance = 1e-10)
  raw$w <- raw$w * 17
  scaled <- survey::svydesign(~psu, strata = ~strata, weights = ~w, nest = TRUE, data = raw)
  a2 <- survey::svyglm(y ~ x * f, scaled, rescale = FALSE)
  b2 <- survey::svyglm(y2 ~ x + f, scaled, rescale = FALSE)
  s2 <- suest(a2, b2, survey_design = scaled, observation_id = c("strata", "id"))
  expect_equal(unname(coef(s2)), unname(coef(s)), tolerance = 1e-10)
  expect_equal(unname(vcov(s2)), unname(vcov(s)), tolerance = 1e-10)
})

test_that("survey independent observations and intercept-only models work", {
  withr::local_options(survey.lonely.psu = "fail", survey.adjust.domain.lonely = FALSE)
  raw <- survey_fixture()$data
  d <- survey::svydesign(~1, weights = ~w, data = raw)
  a <- survey::svyglm(y ~ 1, d)
  b <- survey::svyglm(y2 ~ 1, d)
  s <- suest(a, b, survey_design = d, observation_id = "id")
  means <- survey::svymean(~y + y2, d)
  expect_equal(unname(coef(s)), unname(coef(means)), tolerance = 1e-10)
  expect_equal(unname(vcov(s)), unname(vcov(means)), tolerance = 1e-10)
})

survey_binary_fixture <- function(fpc = FALSE) {
  raw <- survey_fixture()$data
  set.seed(8042)
  raw$f <- factor(raw$unit %% 3)
  raw$off <- (raw$unit - 3.5) / 8
  raw$yb <- rbinom(nrow(raw), 1, plogis(-.35 + .75 * raw$x - .4 * raw$z + raw$off))
  raw$yb2 <- rbinom(nrow(raw), 1, pnorm(.15 - .5 * raw$x + .55 * raw$z - raw$off / 2))
  design <- if (fpc) survey::svydesign(~psu, strata = ~strata, weights = ~w, fpc = ~pop, data = raw) else
    survey::svydesign(~psu, strata = ~strata, weights = ~w, data = raw)
  list(data = raw, design = design)
}

survey_binary_influence <- function(model) {
  X <- stats::model.matrix(model)
  offset <- stats::model.offset(model$model)
  if (is.null(offset)) offset <- rep(0, nrow(X))
  .suest_survey_binary_influence(model, X, model$y, model$prior.weights,
    stats::coef(model), offset)
}

test_that("survey logit/probit influences reproduce survey's native linearization", {
  withr::local_options(survey.lonely.psu = "fail", survey.adjust.domain.lonely = FALSE)
  d <- survey_binary_fixture()$design
  for (link in c("logit", "probit")) {
    a <- survey::svyglm(yb ~ x + z + offset(off), d,
      family = quasibinomial(link), influence = TRUE,
      control = list(epsilon = 1e-12, maxit = 100))
    b <- survey::svyglm(yb ~ x + z + offset(off), d,
      family = quasibinomial(link), influence = TRUE, rescale = FALSE,
      control = list(epsilon = 1e-12, maxit = 100))
    expect_equal(unname(survey_binary_influence(a)), unname(attr(a, "influence")), tolerance = 2e-10)
    expect_equal(unname(survey_binary_influence(b)), unname(attr(b, "influence")), tolerance = 2e-10)
    expect_lt(max(abs(unname(coef(a)) - unname(coef(b)))), 5e-8)
    expect_lt(max(abs(unname(attr(a, "influence")) -
      unname(attr(b, "influence")))), 1e-7)
  }
})

test_that("survey binary covariance reproduces native blocks and independent PSU algebra", {
  withr::local_options(survey.lonely.psu = "fail", survey.adjust.domain.lonely = FALSE)
  for (fpc in c(FALSE, TRUE)) {
    d <- survey_binary_fixture(fpc)$design
    a <- survey::svyglm(yb ~ x, d, family = quasibinomial("logit"), influence = TRUE)
    b <- survey::svyglm(yb2 ~ x + z, d, family = quasibinomial("probit"), influence = TRUE)
    s <- suest(a, b, model_names = c("A", "B"), survey_design = d, observation_id = "id")
    expect_equal(unname(vcov(s)), survey_manual_joint(list(a, b), d), tolerance = 2e-10)
    expect_equal(unname(vcov(s)[s$index[[1]], s$index[[1]]]), unname(vcov(a)), tolerance = 2e-10)
    expect_equal(unname(vcov(s)[s$index[[2]], s$index[[2]]]), unname(vcov(b)), tolerance = 2e-10)
    expect_equal(unname(coef(s)), c(unname(coef(a)), unname(coef(b))))
    expect_equal(unname(s$model_types), c("logit", "probit"))
    expect_equal(s$comparison_scale, "predicted probabilities")
    expect_true(s$mixed_models)
  }
})

test_that("survey binary partial overlap and disjoint domains align on the full design", {
  withr::local_options(survey.lonely.psu = "fail", survey.adjust.domain.lonely = FALSE)
  d <- survey_binary_fixture(TRUE)$design
  subsets <- list(list(subset(d, unit <= 4), subset(d, unit >= 3), FALSE),
    list(subset(d, unit <= 3), subset(d, unit > 3), FALSE),
    list(subset(d, strata <= 2), subset(d, strata > 2), TRUE))
  for (sub in subsets) {
    a <- survey::svyglm(yb ~ x, sub[[1]],
      family = quasibinomial("logit"), influence = TRUE)
    b <- survey::svyglm(yb2 ~ x + z, sub[[2]],
      family = quasibinomial("probit"), influence = TRUE)
    s <- suest(a, b, survey_design = d, observation_id = "id")
    expect_equal(unname(vcov(s)), survey_manual_joint(list(a, b), d), tolerance = 2e-10)
    expect_equal(s$nobs_overlap, length(intersect(a$survey.design$variables$id, b$survey.design$variables$id)))
    cross <- vcov(s)[s$index[[1]], s$index[[2]]]
    if (sub[[3]]) expect_equal(unname(cross), matrix(0, 2, 3), tolerance = 1e-12) else
      expect_gt(max(abs(cross)), 1e-7)
  }
})

test_that("survey binary factors, formula offsets, predictions, and weight rescaling work", {
  withr::local_options(survey.lonely.psu = "fail", survey.adjust.domain.lonely = FALSE)
  fixture <- survey_binary_fixture()
  raw <- fixture$data
  d <- fixture$design
  a <- survey::svyglm(yb ~ x * f + offset(off), d,
    family = quasibinomial("logit"), influence = TRUE,
    control = list(epsilon = 1e-12, maxit = 100))
  b <- survey::svyglm(yb2 ~ x + f + offset(off), d,
    family = quasibinomial("logit"), influence = TRUE,
    control = list(epsilon = 1e-12, maxit = 100))
  s <- suest(a, b, survey_design = d, observation_id = "id")
  nd <- raw[1:6, ]
  analytic <- function(model, data, response = TRUE) {
    terms <- stats::delete.response(stats::terms(model))
    mf <- stats::model.frame(terms, data, na.action = stats::na.pass, xlev = model$xlevels)
    X <- stats::model.matrix(terms, mf, contrasts.arg = model$contrasts)
    offset <- stats::model.offset(mf)
    if (is.null(offset)) offset <- rep(0, nrow(X))
    eta <- as.numeric(X %*% coef(model)) + offset
    if (response) model$family$linkinv(eta) else eta
  }
  expected_response <- c(analytic(a, nd), analytic(b, nd))
  expect_equal(get_predict.suest_model(s, nd)$estimate, expected_response, tolerance = 1e-10)
  expect_equal(marginaleffects::predictions(s, newdata = nd)$estimate, expected_response, tolerance = 1e-8)
  expect_equal(get_predict.suest_model(s, nd, type = "link")$estimate,
    c(analytic(a, nd, FALSE), analytic(b, nd, FALSE)), tolerance = 1e-10)
  perturbed <- coef(s)
  perturbed[1] <- perturbed[1] + 1
  s_perturbed <- set_coef.suest_model(s, perturbed)
  expect_equal(get_predict.suest_model(s_perturbed, nd, type = "link")$estimate -
    get_predict.suest_model(s, nd, type = "link")$estimate, c(rep(1, 6), rep(0, 6)), tolerance = 1e-10)
  expect_equal(unname(vcov(s)), survey_manual_joint(list(a, b), d), tolerance = 2e-10)

  raw$w <- raw$w * 17
  scaled <- survey::svydesign(~psu, strata = ~strata, weights = ~w, data = raw)
  a2 <- survey::svyglm(yb ~ x * f + offset(off), scaled,
    family = quasibinomial("logit"), influence = TRUE, rescale = FALSE,
    control = list(epsilon = 1e-12, maxit = 100))
  b2 <- survey::svyglm(yb2 ~ x + f + offset(off), scaled,
    family = quasibinomial("logit"), influence = TRUE, rescale = FALSE,
    control = list(epsilon = 1e-12, maxit = 100))
  s2 <- suest(a2, b2, survey_design = scaled, observation_id = "id")
  expect_lt(max(abs(unname(coef(s2)) - unname(coef(s)))), 1e-8)
  expect_lt(max(abs(unname(vcov(s2)) - unname(vcov(s)))), 1e-7)
})

test_that("survey binary restrictions remain narrow", {
  withr::local_options(survey.lonely.psu = "fail", survey.adjust.domain.lonely = FALSE)
  fixture <- survey_binary_fixture()
  raw <- fixture$data
  d <- fixture$design
  logit <- survey::svyglm(yb ~ x, d, family = quasibinomial("logit"))
  expect_no_error(suest(logit, logit, survey_design = d, observation_id = "id"))
  wrong <- logit
  wrong$cov.unscaled <- wrong$cov.unscaled * 1.01
  expect_error(suest(logit, wrong, survey_design = d, observation_id = "id"), "native design-based covariance")
  cloglog <- survey::svyglm(yb ~ x, d, family = quasibinomial("cloglog"))
  expect_error(suest(logit, cloglog, survey_design = d, observation_id = "id"), "quasibinomial logit/probit")
  binomial <- suppressWarnings(survey::svyglm(yb ~ x, d, family = binomial("logit")))
  expect_error(suest(logit, binomial, survey_design = d, observation_id = "id"), "quasibinomial logit/probit")
  linear <- survey::svyglm(y ~ x, d)
  expect_error(suest(linear, logit, survey_design = d, observation_id = "id"), "does not combine Gaussian and binary")
  raw$frac <- (raw$yb + .25) / 1.5
  frac_design <- survey::svydesign(~psu, strata = ~strata, weights = ~w, data = raw)
  frac <- survey::svyglm(frac ~ x, frac_design, family = quasibinomial("logit"))
  expect_error(suest(frac, frac, survey_design = frac_design, observation_id = "id"), "numeric 0/1 response")
})
