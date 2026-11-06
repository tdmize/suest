svy_models_fixture <- function() {
  skip_if_not_installed("survey")
  skip_if_not_installed("svyVGAM")
  skip_if_not_installed("numDeriv")
  set.seed(5130)
  n <- 480
  d <- data.frame(id = seq_len(n), strata = rep(1:4, each = n / 4),
    psu = rep(1:32, each = n / 32))
  d$x <- rnorm(n) + rep(rnorm(32), each = n / 32)
  d$z <- rnorm(n)
  d$g <- factor(sample(c("u", "v", "w"), n, TRUE))
  d$w <- 1 + (d$id %% 7) / 3
  d$o <- cut(.8 * d$x - .5 * d$z + rlogis(n), c(-Inf, -1, .5, 1.5, Inf),
    labels = c("a", "b", "c", "d"), ordered_result = TRUE)
  d$m <- factor(ifelse(d$x > .8, "r", sample(c("p", "q", "r"), n, TRUE)))
  d$cnt <- rnbinom(n, size = 1.5, mu = exp(.5 + .4 * d$x - .3 * d$z))
  d$t <- rweibull(n, 1.4, exp(1 + .3 * d$x))
  d$ev <- rbinom(n, 1, .8)
  d$inst <- d$x + rnorm(n)
  d$endo <- .6 * d$inst + .3 * d$z + rnorm(n)
  d$yiv <- 1 + .5 * d$endo - .4 * d$z + rnorm(n)
  list(data = d, design = survey::svydesign(~psu, strata = ~strata, weights = ~w, data = d))
}

svy_tight <- function() VGAM::vglm.control(epsilon = 1e-13, maxit = 200)

# Independent influence: numerical scores and observed information of the
# weighted log likelihood, aligned to the full design by ID.
svy_numeric_influence <- function(loglik, theta, w, ids, design) {
  J <- numDeriv::jacobian(loglik, theta)
  H <- numDeriv::hessian(function(b) sum(w * loglik(b)), theta)
  out <- matrix(0, nrow(design$variables), length(theta))
  out[match(ids, design$variables$id), ] <- (w * J) %*% solve(-H)
  out
}

# Explicit PSU totals and stratum centering, with no survey variance call.
svy_manual_vcov <- function(z, design) {
  V <- matrix(0, ncol(z), ncol(z))
  for (h in unique(design$variables$strata)) {
    keep <- design$variables$strata == h
    totals <- rowsum(z[keep, , drop = FALSE], design$variables$psu[keep])
    G <- nrow(totals)
    totals <- sweep(totals, 2, colMeans(totals))
    V <- V + crossprod(totals) * G / (G - 1)
  }
  V
}

svy_ordered_loglik <- function(X, y, k) function(b) {
  eta <- drop(X %*% b[seq_len(ncol(X))])
  cut <- c(-Inf, b[ncol(X) + seq_len(k)], Inf)
  log(plogis(cut[y + 1] - eta) - plogis(cut[y] - eta))
}

test_that("survey ordered, multinomial, negative-binomial, and gologit blocks equal observed-information sandwiches", {
  withr::local_options(survey.lonely.psu = "fail", survey.adjust.domain.lonely = FALSE)
  f <- svy_models_fixture()
  des <- f$design
  d <- des$variables
  ol <- survey::svyolr(o ~ x + z, des, control = list(reltol = 1e-15, maxit = 1000))
  ml <- svyVGAM::svy_vglm(m ~ x + z, VGAM::multinomial(refLevel = 1), design = des, control = svy_tight())
  nb <- svyVGAM::svy_vglm(cnt ~ x + z, VGAM::negbinomial(), design = des, control = svy_tight())
  go <- svyVGAM::svy_vglm(o ~ x + z, VGAM::cumulative(parallel = FALSE, reverse = TRUE),
    design = des, control = svy_tight())
  s <- suest(ol, ml, nb, go, survey_design = des, observation_id = "id")
  expect_equal(unname(s$model_types), c("ologit", "multinom", "negbin", "gologit"))

  X <- cbind(d$x, d$z)
  X1 <- cbind(1, d$x, d$z)
  y <- as.integer(d$o)
  ym <- as.integer(d$m)
  multinomial <- function(b) {
    B <- cbind(0, matrix(b, 3, 2, byrow = TRUE))
    eta <- X1 %*% B
    eta[cbind(seq_along(ym), ym)] - log(rowSums(exp(eta)))
  }
  negbin <- function(b) dnbinom(d$cnt, size = exp(b[2]), mu = exp(drop(X1 %*% b[c(1, 3, 4)])), log = TRUE)
  gologit <- function(b) {
    above <- cbind(1, plogis(X1 %*% matrix(b, 3, 3, byrow = TRUE)), 0)
    log(above[cbind(seq_along(y), y)] - above[cbind(seq_along(y), y + 1)])
  }
  z <- cbind(
    svy_numeric_influence(svy_ordered_loglik(X, y, 3), coef(ol), d$w, d$id, des),
    svy_numeric_influence(multinomial, coef(ml), d$w, d$id, des),
    svy_numeric_influence(negbin, coef(nb), d$w, d$id, des),
    svy_numeric_influence(gologit, coef(go), d$w, d$id, des))
  expect_equal(unname(vcov(s)), svy_manual_vcov(z, des), tolerance = 1e-5)
  expect_equal(unname(coef(s)), unname(c(coef(ol), coef(ml), coef(nb), coef(go))))
  # The fits' own covariances use approximate or expected information.
  expect_lt(max(s$survey$native_vcov_max_abs_error[1:2]), 1e-7)
  expect_gt(s$survey$native_vcov_max_abs_error[[3]], 1e-5)
})

test_that("survey survreg and instrumental-variable blocks equal the native fits", {
  skip_if_not_installed("AER")
  withr::local_options(survey.lonely.psu = "fail", survey.adjust.domain.lonely = FALSE)
  f <- svy_models_fixture()
  des <- f$design
  sr <- survey::svysurvreg(survival::Surv(t, ev) ~ x + z, des, dist = "weibull")
  tb <- survey::svysurvreg(survival::Surv(pmax(yiv, 0), yiv > 0, type = "left") ~ x + z, des,
    dist = "gaussian")
  iv <- survey::svyivreg(yiv ~ endo + z | inst + z, des)
  lin <- survey::svyglm(yiv ~ endo + z, des)
  s <- suest(sr, tb, iv, lin, survey_design = des, observation_id = "id")
  for (i in 1:4)
    expect_equal(unname(vcov(s)[s$index[[i]], s$index[[i]]]), unname(vcov(list(sr, tb, iv, lin)[[i]])),
      tolerance = 1e-10)
  expect_true(all(c("Log(scale)") %in% s$local_names[[1]]))
  # Row names are the IDs in this fixture. svysurvreg() replaces $var with the
  # design-based covariance after computing its influences.
  dfb <- function(m) {
    m$var <- m$inv.info
    out <- matrix(0, nrow(des$variables), length(coef(m)) + 1L)
    out[as.numeric(rownames(m$model)), ] <- residuals(m, "dfbeta", weighted = TRUE)
    out
  }
  ivz <- matrix(0, nrow(des$variables), 3)
  ivz[as.numeric(rownames(iv$model)), ] <- (model.matrix(iv, component = "projected") *
    (weights(iv) * residuals(iv))) %*% iv$invinf
  z <- cbind(dfb(sr), dfb(tb), ivz, as.matrix(attr(survey::svyglm(yiv ~ endo + z, des, influence = TRUE), "influence")))
  expect_equal(unname(vcov(s)), unname(svy_manual_vcov(z, des)), tolerance = 1e-8)
})

test_that("survey models align domains, missing data, and the common design", {
  withr::local_options(survey.lonely.psu = "fail", survey.adjust.domain.lonely = FALSE)
  f <- svy_models_fixture()
  des <- f$design
  des$variables$z[des$variables$id %% 11 == 0] <- NA
  a <- survey::svyolr(o ~ x, subset(des, strata <= 3), control = list(reltol = 1e-15, maxit = 1000))
  b <- survey::svyolr(o ~ x + z, des, control = list(reltol = 1e-15, maxit = 1000))
  m <- svyVGAM::svy_vglm(m ~ x + z, VGAM::multinomial(refLevel = 1), design = subset(des, strata >= 2),
    control = svy_tight())
  s <- suest(a, b, m, survey_design = des, observation_id = "id")
  expect_equal(unname(s$nobs_models), c(nrow(a$survey.design), nrow(b$survey.design), nrow(m$design)))
  expect_equal(s$nobs_overlap, length(intersect(a$survey.design$variables$id,
    intersect(b$survey.design$variables$id, m$design$variables$id))))
  expect_lt(max(abs(vcov(s)[s$index[[1]], s$index[[1]]] - vcov(a))), 1e-7)
  expect_lt(max(abs(vcov(s)[s$index[[3]], s$index[[3]]] - vcov(m))), 1e-7)
  expect_gt(max(abs(vcov(s)[s$index[[1]], s$index[[3]]])), 1e-5)
  expect_equal(nrow(suest_newdata(s)), sum(s$nobs_models))
})

test_that("survey model predictions follow perturbed coefficients", {
  skip_if_not_installed("AER")
  withr::local_options(survey.lonely.psu = "fail", survey.adjust.domain.lonely = FALSE)
  f <- svy_models_fixture()
  des <- f$design
  ol <- survey::svyolr(o ~ x + g, des)
  ml <- svyVGAM::svy_vglm(m ~ x + g, VGAM::multinomial(refLevel = 1), design = des)
  nb <- svyVGAM::svy_vglm(cnt ~ x, VGAM::negbinomial(), design = des)
  sr <- survey::svysurvreg(survival::Surv(t, ev) ~ x, des, dist = "weibull")
  iv <- survey::svyivreg(yiv ~ endo + z | inst + z, des)
  s <- suest(ol, ml, nb, sr, iv, survey_design = des, observation_id = "id")
  nd <- des$variables[1:4, ]
  p <- get_predict.suest_model(s, nd)
  eta <- drop(model.matrix(~ x + g, nd)[, -1] %*% coef(ol, FALSE))
  expect_equal(p$estimate[p$group == "ol::a"], unname(plogis(ol$zeta[1] - eta)), tolerance = 1e-10)
  B <- cbind(0, matrix(coef(ml), 4, 2, byrow = TRUE))
  E <- exp(model.matrix(~ x + g, nd) %*% B)
  expect_equal(p$estimate[p$group == "ml::q"], unname(E[, 2] / rowSums(E)), tolerance = 1e-10)
  expect_equal(p$estimate[p$group == "nb"], unname(exp(coef(nb)[1] + coef(nb)[3] * nd$x)), tolerance = 1e-10)
  expect_equal(p$estimate[p$group == "sr"], unname(predict(sr, nd)), tolerance = 1e-10)
  expect_equal(p$estimate[p$group == "iv"], unname(drop(cbind(1, nd$endo, nd$z) %*% coef(iv))), tolerance = 1e-10)
  b <- coef(s)
  b["nb::(Intercept):1"] <- b["nb::(Intercept):1"] + 1
  b["ol::x"] <- b["ol::x"] + 1
  p2 <- get_predict.suest_model(set_coef.suest_model(s, b), nd)
  expect_equal(p2$estimate[p2$group == "nb"], exp(1) * p$estimate[p$group == "nb"], tolerance = 1e-10)
  expect_equal(p2$estimate[p2$group == "ol::a"], unname(plogis(ol$zeta[1] - eta - nd$x)), tolerance = 1e-10)
  effects <- marginaleffects::avg_comparisons(s, variables = "x", newdata = suest_newdata(s),
    wts = ".suest_weight")
  expect_equal(nrow(effects), 4 + 3 + 3)
  expect_true(all(is.finite(effects$std.error[effects$group != "iv"])))
  groups <- marginaleffects::avg_comparisons(s, variables = "g", newdata = suest_newdata(s))
  expect_true(all(c("ol::d", "ml::r") %in% groups$group))
})

test_that("unsupported survey model forms are refused", {
  withr::local_options(survey.lonely.psu = "fail", survey.adjust.domain.lonely = FALSE)
  f <- svy_models_fixture()
  des <- f$design
  ok <- survey::svyolr(o ~ x, des)
  expect_error(suest(ok, survey::svyolr(o ~ x, des, method = "cloglog"), survey_design = des,
    observation_id = "id"), "logistic or probit")
  pois <- svyVGAM::svy_vglm(cnt ~ x, VGAM::poissonff(), design = des)
  expect_error(suest(ok, pois, survey_design = des, observation_id = "id"), "cumulative\\(\\), multinomial")
  nb <- svyVGAM::svy_vglm(cnt ~ x, VGAM::negbinomial(zero = NULL), design = des)
  expect_error(suest(ok, nb, survey_design = des, observation_id = "id"), "intercept-only size")
  expect_error(suest(ok, MASS::polr(o ~ x, f$data), survey_design = des, observation_id = "id"),
    "only survey models")
})
