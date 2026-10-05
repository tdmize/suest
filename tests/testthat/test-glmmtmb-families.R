glmm_family_data <- function(groups = 50L, size = 6L) {
  set.seed(1006)
  id <- rep(seq_len(groups), each = size)
  x <- rnorm(groups*size)
  z <- rnorm(groups*size)
  u <- rep(rnorm(groups, sd = 0.7), each = size)
  data.frame(
    id = factor(id), time = rep(seq_len(size), groups), x = x, z = z,
    yb = rbinom(groups*size, 1L, pnorm(-0.2 + 0.5*x - 0.3*z + u)),
    yc = 1 + 0.5*x - 0.3*z + u + rnorm(groups*size),
    yg = rgamma(groups*size, shape = 3, rate = 3/exp(0.3 + 0.2*x - 0.1*z + 0.5*u))
  )
}

glmm_family_fit <- function(family, data, formula = NULL) {
  switch(family,
    probit = glmmTMB::glmmTMB(yb ~ x + z + (1 | id), data = data,
      family = stats::binomial("probit")),
    cloglog = glmmTMB::glmmTMB(yb ~ x + z + (1 | id), data = data,
      family = stats::binomial("cloglog")),
    gaussian = glmmTMB::glmmTMB(yc ~ x + z + (1 | id), data = data,
      family = stats::gaussian()),
    gamma = glmmTMB::glmmTMB(yg ~ x + z + (1 | id), data = data,
      family = stats::Gamma("log")))
}

# Log density of one observation and its first two derivatives in eta.
glmm_family_density <- function(family, y, eta, dispersion) {
  if (family == "probit") {
    lambda1 <- exp(stats::dnorm(eta, log = TRUE) - stats::pnorm(eta, log.p = TRUE))
    lambda0 <- exp(stats::dnorm(eta, log = TRUE) - stats::pnorm(-eta, log.p = TRUE))
    list(l = ifelse(y == 1, stats::pnorm(eta, log.p = TRUE), stats::pnorm(-eta, log.p = TRUE)),
      d1 = ifelse(y == 1, lambda1, -lambda0),
      d2 = ifelse(y == 1, -lambda1*(eta + lambda1), -lambda0*(lambda0 - eta)))
  } else if (family == "cloglog") {
    t <- exp(eta)
    list(l = ifelse(y == 1, log(-expm1(-t)), -t),
      d1 = ifelse(y == 1, t/expm1(t), -t),
      d2 = ifelse(y == 1, t*(expm1(t) - t*exp(t))/expm1(t)^2, -t))
  } else if (family == "gaussian") {
    s <- exp(dispersion)
    list(l = stats::dnorm(y, eta, s, log = TRUE), d1 = (y - eta)/s^2,
      d2 = rep(-1/s^2, length(y)))
  } else {
    k <- exp(dispersion)
    list(l = stats::dgamma(y, shape = k, rate = k/exp(eta), log = TRUE),
      d1 = -k + k*y*exp(-eta), d2 = -k*y*exp(-eta))
  }
}

# Independent Laplace log likelihood of each group's random intercept.
glmm_family_loglik <- function(parameters, model, family, panels = NULL) {
  frame <- suest:::.suest_model_frame(model, "glmmTMB::glmmTMB")
  panel <- as.character(suest:::.suest_panel_id(model))
  if (is.null(panels)) panels <- unique(panel)
  X <- stats::model.matrix(stats::delete.response(stats::terms(model)), frame)
  beta <- parameters[colnames(X)]
  eta <- as.numeric(X %*% beta)
  y <- as.numeric(stats::model.response(frame))
  dispersion <- if (length(parameters) > ncol(X) + 1L) parameters[ncol(X) + 1L] else NA
  sigma2 <- exp(2*unname(parameters["log_sigma"]))
  sum(vapply(panels, function(value) {
    rows <- panel == value
    b <- 0
    for (iteration in 1:100) {
      parts <- glmm_family_density(family, y[rows], eta[rows] + b, dispersion)
      step <- (sum(parts$d1) - b/sigma2)/(sum(parts$d2) - 1/sigma2)
      b <- b - step
      if (abs(step) < 1e-13) break
    }
    parts <- glmm_family_density(family, y[rows], eta[rows] + b, dispersion)
    sum(parts$l) - b^2/(2*sigma2) - 0.5*log(sigma2) -
      0.5*log(1/sigma2 - sum(parts$d2))
  }, numeric(1)))
}

test_that("glmmTMB probit, cloglog, Gaussian, and Gamma scores match their Laplace likelihoods", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  dat <- glmm_family_data()
  for (family in c("probit", "cloglog", "gaussian", "gamma")) {
    model <- glmm_family_fit(family, dat)
    type <- paste0("glmm_", family, "_ri")
    expect_identical(suest:::.suest_model_adapter(model)$type, type)
    component <- suest:::.suest_model_components(model, type, "glmmTMB::glmmTMB")
    parameters <- component$parameters
    expect_equal(glmm_family_loglik(parameters, model, family),
      as.numeric(stats::logLik(model)), tolerance = 1e-8, label = family)

    native <- as.matrix(stats::vcov(model, full = TRUE))
    if (family == "gaussian") {
      # Stata mixed's layout: fixed effects (X'V^-1 X)^-1, variance parameters
      # from the full inverse, no cross block.
      fixed <- seq_len(length(glmmTMB::fixef(model)$cond))
      expected <- matrix(0, nrow(native), ncol(native))
      expected[fixed, fixed] <- solve(solve(native)[fixed, fixed])
      expected[-fixed, -fixed] <- native[-fixed, -fixed]
      native <- expected
    }
    expect_equal(unname(component$bread/nrow(component$score)), unname(native),
      tolerance = 1e-10)

    panel <- as.character(suest:::.suest_panel_id(model))
    first <- unique(panel)[1L]
    step <- .Machine$double.eps^(1/3)*pmax(1, abs(parameters))
    numerical <- vapply(seq_along(parameters), function(j) {
      plus <- minus <- parameters
      plus[j] <- plus[j] + step[j]
      minus[j] <- minus[j] - step[j]
      (glmm_family_loglik(plus, model, family, first) -
        glmm_family_loglik(minus, model, family, first))/(2*step[j])
    }, numeric(1))
    expect_equal(colSums(component$score[panel == first, , drop = FALSE]),
      numerical, tolerance = 2e-6, ignore_attr = TRUE, label = family)
  }
})

test_that("new glmmTMB families name their dispersion parameters", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  dat <- glmm_family_data()
  gaussian <- glmm_family_fit("gaussian", dat)
  gamma <- glmm_family_fit("gamma", dat)
  fit <- suest(gaussian, gamma, model_names = c("mixed", "gamma"),
    observation_id = c("id", "time"))
  expect_identical(names(coef(fit)), c(
    "mixed::(Intercept)", "mixed::x", "mixed::z", "mixed::log_sigma_e", "mixed::log_sigma",
    "gamma::(Intercept)", "gamma::x", "gamma::z", "gamma::log_shape", "gamma::log_sigma"))
  expect_equal(unname(coef(fit)["mixed::log_sigma_e"]), log(stats::sigma(gaussian)),
    tolerance = 1e-10)
  expect_equal(unname(coef(fit)["gamma::log_shape"]), -2*log(stats::sigma(gamma)),
    tolerance = 1e-10)
})

test_that("population-averaged predictions integrate over the random intercept", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  dat <- glmm_family_data()
  rows <- dat[c(1, 50, 200), ]
  for (family in c("probit", "cloglog", "gaussian", "gamma")) {
    model <- glmm_family_fit(family, dat)
    fit <- suest(model, model, model_names = c("a", "b"))
    predicted <- marginaleffects::predictions(fit, newdata = rows)$estimate[1:3]
    beta <- glmmTMB::fixef(model)$cond
    eta <- drop(cbind(1, rows$x, rows$z) %*% beta)
    sigma <- exp(model$fit$par[["theta"]])
    mean_function <- switch(family,
      probit = stats::pnorm, cloglog = function(v) -expm1(-exp(v)),
      gaussian = identity, gamma = exp)
    expected <- vapply(eta, function(e) stats::integrate(function(u)
      mean_function(e + sigma*u)*stats::dnorm(u), -30, 30, rel.tol = 1e-12)$value,
      numeric(1))
    expect_equal(predicted, expected, tolerance = 1e-9, label = family)
  }
})

test_that("new families combine with ordinary models and support random slopes", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  dat <- glmm_family_data()
  probit <- glmm_family_fit("probit", dat)
  ordinary <- glm(yb ~ x + z, binomial("probit"), dat)
  fit <- suest(ordinary, probit, model_names = c("probit", "meprobit"),
    observation_id = c("id", "time"))
  expect_identical(fit$comparison_scale, "predicted probabilities")
  effect <- marginaleffects::avg_comparisons(fit, variables = "x",
    hypothesis = difference ~ revpairwise)
  expect_true(is.finite(effect$std.error) && effect$std.error > 0)

  slope <- glmmTMB::glmmTMB(yc ~ x + z + (1 + x | id), data = dat,
    family = stats::gaussian())
  expect_identical(suest:::.suest_model_adapter(slope)$type, "glmm_gaussian_rs")
  sfit <- suest(slope, glmm_family_fit("gaussian", dat), model_names = c("slope", "intercept"),
    observation_id = c("id", "time"))
  expect_identical(names(coef(sfit))[1:7], c("slope::(Intercept)", "slope::x", "slope::z",
    "slope::log_sigma_e", "slope::log_sd_intercept", "slope::log_sd_slope", "slope::atanh_rho"))
  rows <- dat[1:3, ]
  expect_equal(marginaleffects::predictions(sfit, newdata = rows)$estimate[1:3],
    drop(cbind(1, rows$x, rows$z) %*% glmmTMB::fixef(slope)$cond), tolerance = 1e-10)
})

test_that("random-slope probit and Gamma predictions integrate over the slope variance", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  dat <- glmm_family_data(groups = 80L)
  set.seed(1007)
  u0 <- rep(rnorm(80, sd = 0.7), each = 6)
  u1 <- rep(rnorm(80, sd = 0.5), each = 6)
  dat$yb <- rbinom(nrow(dat), 1L, pnorm(-0.2 + 0.5*dat$x - 0.3*dat$z + u0 + u1*dat$x))
  dat$yg <- rgamma(nrow(dat), shape = 3,
    rate = 3/exp(0.3 + 0.2*dat$x - 0.1*dat$z + u0 + u1*dat$x))
  rows <- dat[c(2, 90, 250), ]
  for (family in c("probit", "gamma")) {
    model <- switch(family,
      probit = glmmTMB::glmmTMB(yb ~ x + z + (1 + x | id), data = dat,
        family = stats::binomial("probit")),
      gamma = glmmTMB::glmmTMB(yg ~ x + z + (1 + x | id), data = dat,
        family = stats::Gamma("log")))
    fit <- suest(model, model, model_names = c("a", "b"))
    p <- coef(fit)
    eta <- drop(cbind(1, rows$x, rows$z) %*% p[c("a::(Intercept)", "a::x", "a::z")])
    sd0 <- exp(p[["a::log_sd_intercept"]]); sd1 <- exp(p[["a::log_sd_slope"]])
    rho <- tanh(p[["a::atanh_rho"]])
    sd <- sqrt((sd0 + rho*sd1*rows$x)^2 + (sd1*rows$x)^2*(1 - rho^2))
    mean_function <- if (family == "probit") stats::pnorm else exp
    expected <- vapply(seq_along(eta), function(i) stats::integrate(function(u)
      mean_function(eta[i] + sd[i]*u)*stats::dnorm(u), -30, 30, rel.tol = 1e-12)$value,
      numeric(1))
    expect_equal(marginaleffects::predictions(fit, newdata = rows)$estimate[1:3],
      expected, tolerance = 1e-9, label = family)
    expect_lt(max(abs(colSums(suest:::.suest_model_components(model,
      paste0("glmm_", family, "_rs"), "glmmTMB::glmmTMB")$score))), 1e-3)
  }
})
