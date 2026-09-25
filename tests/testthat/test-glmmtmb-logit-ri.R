glmmtmb_logit_ri_data <- function(groups = 65L, periods = 6L) {
  set.seed(916)
  id <- rep(seq_len(groups), each = periods)
  x <- rnorm(groups*periods)
  z <- rep(rnorm(groups), each = periods)
  woman <- rep(rep(0:1, length.out = groups), each = periods)
  u1 <- rep(rnorm(groups, sd = 0.8), each = periods)
  u2 <- rep(rnorm(groups, sd = 0.65), each = periods)
  data.frame(
    id = factor(id), time = rep(seq_len(periods), groups),
    higher = factor(rep(rep(seq_len(ceiling(groups/5)), each = 5),
      each = periods, length.out = groups*periods)),
    x = x, z = z, woman = factor(woman),
    y1 = rbinom(groups*periods, 1,
      plogis(-0.35 + 0.7*x - 0.25*z + 0.2*woman + u1)),
    y2 = rbinom(groups*periods, 1,
      plogis(-0.15 + 0.5*x - 0.15*z - 0.1*woman + u2))
  )
}

glmmtmb_logit_ri_fit <- function(formula, data, ...) {
  glmmTMB::glmmTMB(
    formula, data = data, family = stats::binomial("logit"), ...
  )
}

glmmtmb_logit_ri_panel_loglik <- function(parameters, model, panels = NULL) {
  frame <- suest:::.suest_model_frame(model, "glmmTMB::glmmTMB")
  panel <- as.character(suest:::.suest_panel_id(model))
  if (is.null(panels)) panels <- unique(panel)
  terms <- stats::delete.response(stats::terms(model))
  X <- stats::model.matrix(
    terms, frame, contrasts.arg = model$modelInfo$contrasts$cond)
  beta_names <- setdiff(names(parameters), "log_sigma")
  eta <- as.numeric(X[, beta_names, drop = FALSE] %*% parameters[beta_names])
  y <- as.numeric(stats::model.response(frame))
  sigma2 <- exp(2*unname(parameters["log_sigma"]))

  sum(vapply(panels, function(value) {
    rows <- panel == value
    mode_score <- function(b)
      sum(y[rows] - stats::plogis(eta[rows] + b)) - b/sigma2
    bhat <- stats::uniroot(mode_score, c(-30, 30), tol = 1e-12)$root
    probability <- stats::plogis(eta[rows] + bhat)
    curvature <- sum(probability*(1 - probability)) + 1/sigma2
    sum(stats::dbinom(y[rows], 1, probability, log = TRUE)) -
      bhat^2/(2*sigma2) - unname(parameters["log_sigma"]) -
      0.5*log(curvature)
  }, numeric(1)))
}

test_that("glmmTMB random-intercept logit scores match its Laplace likelihood", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  dat <- glmmtmb_logit_ri_data()
  model <- glmmtmb_logit_ri_fit(y1 ~ x + z + (1 | id), dat)
  adapter <- suest:::.suest_model_adapter(model)
  component <- suest:::.suest_glmmtmb_logit_ri_components(model)
  parameters <- component$parameters
  panel <- as.character(suest:::.suest_panel_id(model))

  expect_identical(adapter$type, "glmm_logit_ri")
  expect_identical(adapter$engine, "glmmTMB::glmmTMB")
  expect_equal(
    glmmtmb_logit_ri_panel_loglik(parameters, model),
    as.numeric(stats::logLik(model)), tolerance = 1e-8)
  expect_lt(max(abs(colSums(component$score))), 5e-4)

  native <- as.matrix(stats::vcov(model, full = TRUE))
  dimnames(native) <- list(names(parameters), names(parameters))
  expect_equal(component$bread/nrow(component$score), native,
    tolerance = 1e-12)

  first_panel <- unique(panel)[1L]
  rows <- which(panel == first_panel)
  step <- .Machine$double.eps^(1/3)*pmax(1, abs(parameters))
  numerical_score <- vapply(seq_along(parameters), function(j) {
    plus <- minus <- parameters
    plus[j] <- plus[j] + step[j]
    minus[j] <- minus[j] - step[j]
    (glmmtmb_logit_ri_panel_loglik(plus, model, first_panel) -
       glmmtmb_logit_ri_panel_loglik(minus, model, first_panel))/(2*step[j])
  }, numeric(1))
  expect_equal(colSums(component$score[rows, , drop = FALSE]),
    numerical_score, tolerance = 2e-6, ignore_attr = TRUE)
  expect_equal(sum(rowSums(abs(component$score)) > 0),
    length(unique(panel)))
})

test_that("glmmTMB logit systems support panel alignment and higher clusters", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  dat <- glmmtmb_logit_ri_data(groups = 70L)
  first <- glmmtmb_logit_ri_fit(y1 ~ x + z + (1 | id), dat)
  second <- glmmtmb_logit_ri_fit(y2 ~ x + z + (1 | id), dat)
  fit <- suest(first, second, model_names = c("First", "Second"),
    observation_id = c("id", "time"))
  higher <- suest(first, second, model_names = c("First", "Second"),
    observation_id = c("id", "time"), cluster = "higher")

  expect_identical(unname(fit$model_types), rep("glmm_logit_ri", 2L))
  expect_true(all(grepl("log_sigma$", names(coef(fit))[c(4L, 8L)])))
  expect_true(all(is.finite(vcov(fit))))
  expect_true(all(is.finite(vcov(higher))))
  expect_equal(vcov(fit), t(vcov(fit)), tolerance = 1e-12)
  expect_gt(max(abs(offdiag_vcov(fit))), 0)
  expect_gt(max(abs(vcov(fit) - vcov(higher))), 0)

  partial_first <- glmmtmb_logit_ri_fit(
    y1 ~ x + z + (1 | id), dat[as.integer(dat$id) <= 60, ])
  partial_second <- glmmtmb_logit_ri_fit(
    y2 ~ x + z + (1 | id), dat[as.integer(dat$id) >= 11, ])
  partial <- suest(partial_first, partial_second,
    observation_id = c("id", "time"))
  expect_true(all(is.finite(vcov(partial))))
  expect_gt(max(abs(offdiag_vcov(partial))), 0)

  low <- glmmtmb_logit_ri_fit(
    y1 ~ x + z + (1 | id), dat[as.integer(dat$id) <= 35, ])
  high <- glmmtmb_logit_ri_fit(
    y2 ~ x + z + (1 | id), dat[as.integer(dat$id) > 35, ])
  disjoint <- suest(low, high, observation_id = c("id", "time"))
  expect_equal(unname(offdiag_vcov(disjoint)), matrix(0, 4L, 4L),
    tolerance = 1e-12)
})

test_that("glmmTMB integrated predictions work with marginaleffects", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  dat <- glmmtmb_logit_ri_data()
  first <- glmmtmb_logit_ri_fit(y1 ~ x*woman + z + (1 | id), dat)
  second <- glmmtmb_logit_ri_fit(y2 ~ x*woman + z + (1 | id), dat)
  fit <- suest(first, second, model_names = c("First", "Second"),
    observation_id = c("id", "time"))
  predictions <- marginaleffects::avg_predictions(fit, newdata = dat)
  slopes <- marginaleffects::avg_slopes(
    fit, variables = "x", newdata = dat)

  expected <- vapply(list(first, second), function(model) {
    mean(suest:::.suest_predict_values(
      model, dat, "response", "glmmTMB::glmmTMB", "glmm_logit_ri"))
  }, numeric(1))
  expect_equal(predictions$estimate, expected, tolerance = 1e-8)
  expect_true(all(is.finite(predictions$std.error)))
  expect_true(all(is.finite(slopes$estimate)))
  expect_true(all(is.finite(slopes$std.error)))

  parameters <- suest:::.suest_extract_parameters(
    first, "glmm_logit_ri", "glmmTMB::glmmTMB")
  eta <- as.numeric(stats::model.matrix(
    stats::delete.response(stats::terms(first)), dat[1:3, ]) %*%
      parameters[setdiff(names(parameters), "log_sigma")])
  sigma <- exp(unname(parameters["log_sigma"]))
  numerical_integral <- vapply(eta, function(value) stats::integrate(
    function(u) stats::plogis(value + sigma*u)*stats::dnorm(u),
    -Inf, Inf, rel.tol = 1e-12
  )$value, numeric(1))
  quadrature_prediction <- suest:::.suest_predict_values(
    first, dat[1:3, ], "response", "glmmTMB::glmmTMB", "glmm_logit_ri")
  expect_equal(quadrature_prediction, numerical_integral, tolerance = 1e-9)

  changed <- coef(fit)
  changed["First::log_sigma"] <- changed["First::log_sigma"] + log(1.5)
  altered <- marginaleffects::set_coef(fit, changed)
  original_prediction <- marginaleffects::predictions(
    fit, newdata = dat[1:5, ], vcov = FALSE)
  altered_prediction <- marginaleffects::predictions(
    altered, newdata = dat[1:5, ], vcov = FALSE)
  first_rows <- original_prediction$group == "First"
  expect_gt(max(abs(
    original_prediction$estimate[first_rows] -
      altered_prediction$estimate[first_rows]
  )), 1e-5)
})

test_that("glmmTMB validation is deliberately narrow", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  dat <- glmmtmb_logit_ri_data(groups = 45L)
  dat$constant_offset <- 0.1
  probit <- glmmTMB::glmmTMB(
    y1 ~ x + z + (1 | id), data = dat,
    family = stats::binomial("probit"))
  poisson <- glmmTMB::glmmTMB(
    y1 ~ x + z + (1 | id), data = dat, family = stats::poisson())
  random_slope <- suppressWarnings(glmmTMB::glmmTMB(
    y1 ~ x + z + (1 + x + z | id), data = dat,
    family = stats::binomial("logit")))
  weighted <- glmmtmb_logit_ri_fit(
    y1 ~ x + z + (1 | id), dat, weights = rep(2, nrow(dat)))
  offset <- glmmtmb_logit_ri_fit(
    y1 ~ x + z + offset(constant_offset) + (1 | id), dat)

  expect_error(suest(probit, probit), "binomial-logit", fixed = TRUE)
  expect_s3_class(suest(poisson, poisson), "suest_model")
  expect_error(suest(random_slope, random_slope),
    "one conditional random intercept", fixed = TRUE)
  expect_error(suest(weighted, weighted), "Weighted glmmTMB", fixed = TRUE)
  expect_error(suest(offset, offset), "Offsets are not supported", fixed = TRUE)
})

test_that("glmmTMB logit reproduces the returned Stata benchmark", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  skip_if_not_installed("marginaleffects")
  fixture <- readRDS(test_path("fixtures", "glmmtmb-logit-ri-stata.rds"))
  dat <- fixture$data
  dat$id <- factor(dat$id)
  fit_model <- function(formula, data = dat) glmmtmb_logit_ri_fit(
    formula, data)
  extract_model <- function(model) {
    parameters <- suest:::.suest_extract_parameters(
      model, "glmm_logit_ri", "glmmTMB::glmmTMB")
    covariance <- as.matrix(vcov(model, full = TRUE))
    dimnames(covariance) <- list(names(parameters), names(parameters))
    list(coefficients = parameters, vcov = covariance,
      logLik = as.numeric(logLik(model)))
  }

  y1 <- fit_model(y1 ~ x + z + (1 | id))
  y2 <- fit_model(y2 ~ x + z + (1 | id))
  for (i in seq_along(list(y1, y2))) {
    actual <- extract_model(list(y1, y2)[[i]])
    expected <- fixture$balanced[[paste0("model", i)]]
    expect_lt(max(abs(actual$coefficients - expected$coefficients)), 2e-5)
    expect_lt(max(abs(actual$vcov - expected$vcov)), 1e-6)
    expect_lt(abs(actual$logLik - expected$logLik), 1e-6)
  }

  balanced <- suest(y1, y2, model_names = c("Y1", "Y2"),
    observation_id = c("id", "time"))
  predictions <- marginaleffects::avg_predictions(balanced, newdata = dat)
  slopes <- marginaleffects::avg_slopes(
    balanced, variables = "x", newdata = dat)
  expect_lt(max(abs(predictions$estimate -
    fixture$balanced$prediction_mean)), 4e-6)
  expect_lt(max(abs(slopes$estimate -
    fixture$balanced$average_slope_x)), 4e-6)

  partial <- suest(
    fit_model(y1_partial ~ x + z + (1 | id)),
    fit_model(y2_partial ~ x + z + (1 | id)),
    model_names = c("Y1", "Y2"), observation_id = c("id", "time"),
    cluster = "higher")
  left <- dat[as.integer(dat$id) <= 40, ]
  right <- dat[as.integer(dat$id) > 40, ]
  disjoint <- suest(
    fit_model(y1_left ~ x + z + (1 | id), left),
    fit_model(y2_right ~ x + z + (1 | id), right),
    model_names = c("Y1", "Y2"), observation_id = c("id", "time"))
  systems <- list(
    balanced = balanced, partial_higher = partial, disjoint = disjoint)
  covariance_tolerance <- c(
    balanced = 7e-4, partial_higher = 1.2e-3, disjoint = 1.1e-3)

  for (case in names(systems)) {
    actual <- systems[[case]]
    expected <- fixture[[case]]$gsem_laplace
    expect_lt(max(abs(coef(actual) - expected$coefficients)), 3e-5)
    expect_lt(max(abs(vcov(actual) - expected$vcov)),
      covariance_tolerance[case])
    expect_identical(actual$nobs_union, expected$nobs_union)
    expect_identical(actual$nobs_overlap, expected$nobs_overlap)
    expect_identical(actual$n_clusters, expected$n_clusters)
  }
  expect_equal(unname(vcov(disjoint)[1:4, 5:8]), matrix(0, 4L, 4L),
    tolerance = 1e-12)
  expect_lt(max(abs(fixture$disjoint$gsem_laplace$vcov[1:4, 5:8])),
    2.1e-4)
  expect_equal(fixture$disjoint$system$vcov[1:4, 5:8],
    matrix(0, 4L, 4L), tolerance = 1e-12, ignore_attr = TRUE)
})
