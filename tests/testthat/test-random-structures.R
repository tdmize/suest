rs_terms_slopes <- function(d, columns = c("intercept", "x1", "x2"), diagonal = FALSE)
  list(rs_term("school", d$school, columns, diagonal))

rs_terms_nested <- function(d, name = "school.class", columns = c("intercept", "x1"))
  list(rs_term("school", d$school, columns),
    rs_term(name, interaction(d$school, d$class), "intercept"))

rs_components <- function(fit) {
  if (inherits(fit, "merMod")) {
    w <- suest:::.suest_wrap_lme4(fit)
    a <- suest:::.suest_lme4_adapter(w)
    engine <- if (inherits(fit, "lmerMod")) "lme4::lmer" else "lme4::glmer"
    c(suest:::.suest_model_components(w, a$type, engine), type = a$type)
  } else if (inherits(fit, "clmm")) {
    c(suest:::.suest_clmm_components(fit), type = suest:::.suest_model_adapter(fit)$type)
  } else {
    a <- suest:::.suest_model_adapter(fit)
    c(suest:::.suest_model_components(fit, a$type, a$engine), type = a$type)
  }
}

test_that("lmer with several slopes or a nested level matches the exact Gaussian likelihood", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("numDeriv")
  d <- random_structure_data()
  cases <- list(
    slopes = list(fit = lme4::lmer(yc ~ x1 + x2 + (1 + x1 + x2 | school), d, REML = FALSE),
      terms = rs_terms_slopes(d)),
    nested = list(fit = lme4::lmer(yc ~ x1 + x2 + (1 + x1 | school) + (1 | school:class), d,
      REML = FALSE), terms = rs_terms_nested(d)))
  for (case in names(cases)) {
    fit <- cases[[case]]$fit
    component <- rs_components(fit)
    expect_identical(component$type, "glmm_gaussian_re")
    theta <- component$parameters
    spec <- rs_spec(d, ~ x1 + x2, d$yc, "gaussian", cases[[case]]$terms)
    total <- function(p) sum(rs_school_loglik(p, spec))
    expect_equal(total(theta), as.numeric(logLik(fit)), tolerance = 1e-10, label = case)
    schools <- c("1", "17", "40")
    expect_equal(unname(rowsum(component$score, d$school)[schools, ]),
      unname(rs_numeric_scores(theta, spec, schools)), tolerance = 1e-6, label = case)
    # Stata mixed's layout: lmer's fixed-effect covariance, and the variance
    # block of the inverse observed information.
    bread <- component$bread/nrow(d)
    fixed <- 1:3
    expect_equal(unname(bread[fixed, fixed]), unname(as.matrix(vcov(fit))), tolerance = 1e-6)
    expect_true(all(bread[fixed, -fixed] == 0))
    inverse <- solve(-numDeriv::hessian(total, theta, method.args = list(r = 2)))
    expect_equal(unname(bread[-fixed, -fixed]), unname(inverse[-fixed, -fixed]),
      tolerance = 1e-4, label = case)
  }
})

test_that("Laplace fits with several random effects match an independent Laplace likelihood", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  d <- random_structure_data()
  tight <- lme4::glmerControl(tolPwrss = 1e-12)
  cases <- list(
    glmer_slopes = list(fit = lme4::glmer(yb ~ x1 + x2 + (1 + x1 + x2 | school), d, binomial,
      control = tight), family = "logit", y = d$yb, terms = rs_terms_slopes(d)),
    glmer_nested = list(fit = lme4::glmer(count ~ x1 + x2 + (1 + x1 | school) + (1 | school:class),
      d, poisson, control = tight), family = "poisson", y = d$count, terms = rs_terms_nested(d)),
    glmer_three_level = list(fit = lme4::glmer(yb ~ x1 + x2 + (1 | school/class), d, binomial,
      control = tight), family = "logit", y = d$yb,
      terms = rs_terms_nested(d, "class.school", "intercept")),
    glmer_independent = list(fit = lme4::glmer(count ~ x1 + x2 + (1 + x1 || school), d, poisson,
      control = tight), family = "poisson", y = d$count,
      terms = list(rs_term("school", d$school, "intercept"), rs_term("school", d$school, "x1"))),
    tmb_slope_only = list(fit = glmmTMB::glmmTMB(count ~ x1 + (0 + x1 | school) + (1 | school:class),
      d, family = poisson), family = "poisson", y = d$count,
      terms = list(rs_term("school", d$school, "x1"),
        rs_term("school.class", interaction(d$school, d$class), "intercept"))),
    tmb_slopes = list(fit = glmmTMB::glmmTMB(count ~ x1 + x2 + (1 + x1 + x2 | school), d,
      family = poisson), family = "poisson", y = d$count, terms = rs_terms_slopes(d)),
    tmb_diag = list(fit = glmmTMB::glmmTMB(yb ~ x1 + x2 + diag(1 + x1 + x2 | school), d,
      family = binomial), family = "logit", y = d$yb, terms = rs_terms_slopes(d, diagonal = TRUE)),
    tmb_nested = list(fit = glmmTMB::glmmTMB(yb ~ x1 + x2 + (1 + x1 | school) + (1 | school:class),
      d, family = binomial), family = "logit", y = d$yb, terms = rs_terms_nested(d)))
  for (case in names(cases)) {
    x <- cases[[case]]
    component <- rs_components(x$fit)
    expect_match(component$type, "^glmm_(logit|poisson)_re$", label = case)
    theta <- component$parameters
    spec <- rs_spec(d, if (case == "tmb_slope_only") ~ x1 else ~ x1 + x2, x$y, x$family, x$terms)
    expect_equal(sum(rs_school_loglik(theta, spec)), as.numeric(logLik(x$fit)),
      tolerance = 1e-7, label = case)
    schools <- c("2", "23", "39")
    expect_equal(unname(rowsum(component$score, d$school)[schools, ]),
      unname(rs_numeric_scores(theta, spec, schools)), tolerance = 1e-6, label = case)
  }
})

test_that("glmmTMB systems equal the sandwich built from glmmTMB's own group scores", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  skip_if_not_installed("numDeriv")
  d <- random_structure_data()
  models <- list(
    glmmTMB::glmmTMB(count ~ x1 + x2 + (1 + x1 + x2 | school), d, family = poisson),
    glmmTMB::glmmTMB(yb ~ x1 + x2 + (1 + x1 | school) + (1 | school:class), d, family = binomial),
    glmmTMB::glmmTMB(yb ~ x2 + diag(1 + x1 + x2 | school), d, family = binomial))
  influence <- function(model) {
    J <- rs_glmmtmb_jacobian(model)
    U <- sandwich::estfun(model, full = TRUE) %*% as.matrix(vcov(model, full = TRUE)) %*% t(J)
    colnames(U) <- c(names(glmmTMB::fixef(model)$cond), rownames(J)[-seq_along(glmmTMB::fixef(model)$cond)])
    U
  }
  fit <- suest(models[[1]], models[[2]], models[[3]], model_names = c("A", "B", "C"))
  expect_identical(unname(fit$model_types), rep(c("glmm_poisson_re", "glmm_logit_re"), c(1, 2)))
  parts <- lapply(models, influence)
  expected <- crossprod(do.call(cbind, parts))*40/39
  # Match suest's parameter order (top level first) to glmmTMB's by name.
  offset <- c(0, cumsum(vapply(parts, ncol, 0)))
  index <- unlist(lapply(seq_along(models), function(i) {
    own <- sub("^[ABC]::", "", names(coef(fit))[fit$index[[i]]])
    offset[i] + match(own, colnames(parts[[i]]))
  }))
  expect_false(anyNA(index))
  # suest's information is numerical; glmmTMB's is exact at its own optimum.
  expect_equal(unname(vcov(fit)), unname(expected[index, index]), tolerance = 1e-4)
  expect_identical(fit$n_clusters, 40L)
})

test_that("general structures give glmmTMB's own group scores for every family", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  skip_if_not_installed("numDeriv")
  d <- random_structure_data()
  fit <- function(...) suppressWarnings(glmmTMB::glmmTMB(..., data = d))
  models <- list(
    logit = fit(yb ~ x1 + x2 + (1 + x1 + x2 | school), family = binomial),
    probit = fit(yb ~ x1 + (1 + x1 | school) + (1 | school:class), family = binomial("probit")),
    cloglog = fit(yb ~ x1 + x2 + diag(1 + x1 + x2 | school), family = binomial("cloglog")),
    poisson = fit(count ~ x1 + (1 + x1 | school) + (1 | school:class), family = poisson),
    nbinom2 = fit(nb ~ x1 + (1 + x1 | school) + (1 | school:class), family = glmmTMB::nbinom2),
    gaussian = fit(yc ~ x1 + x2 + (1 + x1 + x2 | school)),
    gamma = fit(positive ~ x1 + (1 + x1 | school) + (1 | school:class), family = Gamma("log")))
  for (family in names(models)) {
    model <- models[[family]]
    component <- rs_components(model)
    expect_identical(component$type, paste0("glmm_", family, "_re"))
    J <- rs_glmmtmb_jacobian(model)
    covariance <- grep("^log_sd_|^atanh_rho_", rownames(J), value = TRUE)
    expect_equal(unname(component$parameters[covariance]),
      unname(rs_glmmtmb_covariance(model, model$fit$par[names(model$fit$par) == "theta"])),
      label = family)
    # Scores transform with the inverse Jacobian: s_suest = J^-T s_glmmTMB.
    theirs <- sandwich::estfun(model, full = TRUE) %*% solve(J)
    ours <- rowsum(component$score, d$school)
    position <- c(seq_len(ncol(J) - length(covariance)),
      match(covariance, names(component$parameters)))
    expect_equal(unname(ours[rownames(theirs), position]), unname(theirs), tolerance = 1e-6,
      label = family)
  }
})

test_that("clmm random slopes and nested levels match an independent Laplace likelihood", {
  skip_if_not_installed("ordinal")
  d <- random_structure_data()
  cases <- list(
    slopes = list(fit = ordinal::clmm(yo ~ x1 + x2 + (1 + x1 + x2 | school), data = d),
      terms = rs_terms_slopes(d)),
    nested = list(fit = ordinal::clmm(yo ~ x1 + x2 + (1 + x1 | school) + (1 | school:class),
      data = d), terms = rs_terms_nested(d)))
  for (case in names(cases)) {
    fit <- cases[[case]]$fit
    component <- rs_components(fit)
    expect_identical(component$type, "re_ologit")
    theta <- component$parameters
    # The threshold model is P(Y <= j) = F(threshold_j - x b - z u).
    spec <- rs_spec(d, ~ 0 + x1 + x2, as.integer(d$yo), "ologit", cases[[case]]$terms,
      thresholds = c("a|b", "b|c", "c|d"))
    fixed_names <- c("a|b", "b|c", "c|d", "x1", "x2")
    # clmm stops its inner mode search at a gradient of 1e-4.
    expect_equal(sum(rs_school_loglik(theta, spec)), as.numeric(logLik(fit)),
      tolerance = 1e-6, label = case)
    schools <- c("5", "31")
    expect_equal(unname(rowsum(component$score, d$school)[schools, ]),
      unname(rs_numeric_scores(theta, spec, schools)), tolerance = 1e-6, label = case)
    # The fixed-effect block of the inverse information does not depend on
    # how the covariance is parameterized.
    bread <- component$bread/nrow(d)
    expect_equal(unname(bread[fixed_names, fixed_names]),
      unname(as.matrix(vcov(fit))[fixed_names, fixed_names]), tolerance = 1e-3, label = case)
  }
})

test_that("glmer and glmmTMB give the same system for the same general structure", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  d <- random_structure_data()
  partner <- lm(yc ~ x1, d)
  # The two packages' optimizers stop at slightly different points, which
  # differ more across platforms for weakly identified variances.
  cases <- list(list(count ~ x1 + x2 + (1 + x1 + x2 | school), poisson),
    list(yb ~ x1 + x2 + (1 + x1 | school) + (1 | school:class), binomial))
  for (case in cases) {
    a <- suest(lme4::glmer(case[[1]], d, case[[2]], control = lme4::glmerControl(tolPwrss = 1e-12)),
      partner, model_names = c("m", "p"))
    b <- suest(glmmTMB::glmmTMB(case[[1]], d, family = case[[2]]), partner, model_names = c("m", "p"))
    expect_identical(names(coef(a)), names(coef(b)))
    expect_equal(coef(a), coef(b), tolerance = 1e-3)
    expect_equal(vcov(a), vcov(b), tolerance = 5e-3)
  }
  linear <- list(lme4::lmer(yc ~ x1 + x2 + (1 + x1 | school) + (1 | school:class), d, REML = FALSE),
    glmmTMB::glmmTMB(yc ~ x1 + x2 + (1 + x1 | school) + (1 | school:class), d))
  a <- suest(linear[[1]], glm(yb ~ x1, binomial, d), model_names = c("m", "p"))
  b <- suest(linear[[2]], glm(yb ~ x1, binomial, d), model_names = c("m", "p"))
  expect_equal(vcov(a), vcov(b), tolerance = 1e-3)
})

test_that("population-averaged predictions integrate over every random effect", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  skip_if_not_installed("ordinal")
  skip_if_not_installed("marginaleffects")
  d <- random_structure_data()
  rows <- d[c(1, 200, 500), ]
  rows$x1 <- c(-2, 0.5, 3)
  variance <- function(fit, rows) {
    vc <- if (inherits(fit, "clmm")) ordinal::VarCorr(fit) else if (inherits(fit, "glmmTMB"))
      glmmTMB::VarCorr(fit)$cond else lme4::VarCorr(fit)
    Reduce(`+`, lapply(vc, function(S) {
      S <- as.matrix(S)
      Z <- sapply(colnames(S), function(col) if (col == "(Intercept)") rep(1, nrow(rows)) else rows[[col]])
      Z <- matrix(Z, nrow(rows))
      rowSums((Z %*% S)*Z)
    }))
  }
  pa <- function(fit, partner = lm(yc ~ x1, d)) {
    system <- suest(fit, partner, model_names = c("m", "p"))
    p <- marginaleffects::predictions(system, newdata = rows, vcov = FALSE)
    list(system = system, estimate = p$estimate[startsWith(as.character(p$group), "m")])
  }
  logit <- glmmTMB::glmmTMB(yb ~ x1 + x2 + (1 + x1 | school) + (1 | school:class), d, family = binomial)
  eta <- unname(drop(model.matrix(~ x1 + x2, rows) %*% glmmTMB::fixef(logit)$cond))
  v <- variance(logit, rows)
  expected <- vapply(seq_along(eta), function(i) integrate(function(u)
    plogis(eta[i] + sqrt(v[i])*u)*dnorm(u), -Inf, Inf, rel.tol = 1e-12)$value, 0)
  expect_equal(pa(logit)$estimate, expected, tolerance = 1e-8)

  probit <- glmmTMB::glmmTMB(yb ~ x1 + x2 + (1 + x1 + x2 | school), d, family = binomial("probit"))
  eta <- unname(drop(model.matrix(~ x1 + x2, rows) %*% glmmTMB::fixef(probit)$cond))
  expect_equal(pa(probit)$estimate, pnorm(eta/sqrt(1 + variance(probit, rows))), tolerance = 1e-8)

  count <- lme4::glmer(count ~ x1 + x2 + (1 + x1 | school) + (1 | school:class), d, poisson)
  eta <- unname(drop(model.matrix(~ x1 + x2, rows) %*% lme4::fixef(count)))
  expect_equal(pa(count)$estimate, exp(eta + variance(count, rows)/2), tolerance = 1e-8)

  ordered <- ordinal::clmm(yo ~ x1 + x2 + (1 + x1 | school) + (1 | school:class), data = d)
  eta <- unname(drop(as.matrix(rows[c("x1", "x2")]) %*% ordered$beta))
  v <- variance(ordered, rows)
  cumulative <- sapply(ordered$alpha, function(t) vapply(seq_along(eta), function(i)
    integrate(function(u) plogis(t - eta[i] - sqrt(v[i])*u)*dnorm(u), -Inf, Inf,
      rel.tol = 1e-12)$value, 0))
  expected <- cbind(cumulative, 1) - cbind(0, cumulative)
  result <- pa(ordered)
  expect_equal(result$estimate, as.vector(expected), tolerance = 1e-8)
  expect_identical(result$system$n_clusters, 40L)
})

test_that("unsupported random-effects structures are refused with a reason", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  d <- random_structure_data()
  d$region <- factor(as.integer(d$school) %% 4)
  d$crossed <- factor(rep(1:5, length.out = nrow(d)))
  partner <- lm(yc ~ x1, d)
  crossed <- suppressMessages(lme4::lmer(yc ~ x1 + (1 | school) + (1 | crossed), d, REML = FALSE))
  expect_error(suest(crossed, partner), "crossed random effects are not supported")
  d$area <- factor(as.integer(d$school) %% 10)
  four <- suppressMessages(lme4::lmer(yc ~ x1 + (1 | area/school/class), d, REML = FALSE))
  expect_error(suest(four, partner), "at most one grouping variable nested")
  transformed <- lme4::lmer(yc ~ x1 + (1 + I(x1^2) | school), d, REML = FALSE)
  expect_error(suest(transformed, partner), "untransformed numeric data columns")
  d$f <- factor(d$x2 > 0)
  factor_slope <- suppressMessages(lme4::lmer(yc ~ x1 + (1 + f | school), d, REML = FALSE))
  expect_error(suest(factor_slope, partner), "untransformed numeric data columns")
  ar1 <- suppressWarnings(glmmTMB::glmmTMB(yc ~ x1 + ar1(class + 0 | school), d))
  expect_error(suest(ar1, partner))
  nested <- lme4::lmer(yc ~ x1 + x2 + (1 + x1 | school) + (1 | school:class), d, REML = FALSE)
  d$school_class <- interaction(d$school, d$class)
  expect_error(suest(nested, partner, cluster = "school_class"))
  expect_s3_class(suest(nested, partner, cluster = "region"), "suest_model")
})

test_that("nearly singular random-effect covariances are refused", {
  skip_if_not_installed("lme4")
  d <- random_structure_data()
  set.seed(9)
  # The x1 slope is an exact multiple of the intercept: correlation 1.
  u <- rnorm(40)[d$school]
  d$yd <- 1 + d$x1 + u*(1 + 0.5*d$x1) + rnorm(40, sd = 0.4)[d$school]*d$x2 + rnorm(nrow(d))
  fit <- suppressMessages(suppressWarnings(lme4::lmer(yd ~ x1 + (1 + x1 + x2 | school), d,
    REML = FALSE)))
  vc <- as.matrix(lme4::VarCorr(fit)$school)
  skip_if(min(eigen(cov2cor(vc), only.values = TRUE)$values) > 1e-4)
  expect_error(suest(fit, lm(yc ~ x1, d)), "singular")
})
