panel_poisson_re_test_data <- function(groups = 70L, periods = 6L) {
  set.seed(41903)
  id <- rep(seq_len(groups), each = periods)
  time <- rep(seq_len(periods), groups)
  higher <- ceiling(id / 5)
  x <- rnorm(length(id))
  z <- rnorm(length(id))
  woman <- factor(rbinom(length(id), 1, 0.48), labels = c("Men", "Women"))
  common <- rgamma(groups, shape = 2.2, rate = 2.2)[id]
  second <- rgamma(groups, shape = 3.0, rate = 3.0)[id]
  mu1 <- exp(0.20 + 0.35*x - 0.20*z + 0.15*(woman == "Women"))*common
  mu2 <- exp(-0.10 + 0.25*x + 0.20*z - 0.10*(woman == "Women"))*
    sqrt(common)*second
  data.frame(
    id, time, higher, x, z, woman,
    y1 = rpois(length(id), mu1), y2 = rpois(length(id), mu2)
  )
}

panel_poisson_re_start <- function(formula, data) {
  c(coef(glm(formula, data = data, family = poisson("log"))), alpha = 0.5)
}

panel_poisson_re_fit <- function(formula, data, other = "sd") {
  maxLik <- maxLik::maxLik
  suppressWarnings(pglm::pglm(
    formula, data = data, family = poisson("log"), model = "random",
    effect = "individual", other = other, index = c("id", "time"),
    method = "bfgs", print.level = 0, reltol = 1e-12,
    start = panel_poisson_re_start(formula, data)
  ))
}

panel_poisson_re_loglik <- function(parameters, model, panels = NULL) {
  fit_environment <- environment(model$objectiveFn)
  X <- get("X", fit_environment, inherits = FALSE)
  y <- as.numeric(get("y", fit_environment, inherits = FALSE))
  id <- as.character(get("id", fit_environment, inherits = FALSE))
  beta_names <- setdiff(names(parameters), "alpha")
  mean <- exp(as.numeric(X[, beta_names, drop = FALSE] %*%
    parameters[beta_names]))
  dispersion <- 1/unname(parameters["alpha"])
  if (is.null(panels)) panels <- unique(id)

  sum(vapply(panels, function(panel) {
    rows <- id == panel
    total_mean <- sum(mean[rows])
    total_y <- sum(y[rows])
    sum(y[rows]*log(mean[rows]) - lgamma(y[rows] + 1)) +
      dispersion*log(dispersion) -
      (dispersion + total_y)*log(total_mean + dispersion) +
      lgamma(dispersion + total_y) - lgamma(dispersion)
  }, numeric(1)))
}

test_that("pglm gamma random-effects Poisson scores match its likelihood", {
  skip_if_not_installed("pglm")
  dat <- panel_poisson_re_test_data()
  model <- panel_poisson_re_fit(y1 ~ x + z, dat)
  adapter <- suest:::.suest_model_adapter(model)
  component <- suest:::.suest_pglm_poisson_re_components(model)
  parameters <- component$parameters
  panel <- as.character(suest:::.suest_panel_id(model))

  expect_identical(adapter$type, "panel_poisson_re")
  expect_named(parameters, c("(Intercept)", "x", "z", "alpha"))
  expect_gt(parameters["alpha"], 0)
  expect_equal(panel_poisson_re_loglik(parameters, model),
    as.numeric(logLik(model)), tolerance = 1e-9)
  expect_lt(max(abs(colSums(component$score))), 2e-5)
  expect_equal(component$bread/nrow(component$score), vcov(model),
    tolerance = 1e-10, ignore_attr = TRUE)

  first_panel <- unique(panel)[1L]
  rows <- which(panel == first_panel)
  step <- 1e-5*pmax(1, abs(parameters))
  numerical_score <- vapply(seq_along(parameters), function(j) {
    plus <- minus <- parameters
    plus[j] <- plus[j] + step[j]
    minus[j] <- minus[j] - step[j]
    (panel_poisson_re_loglik(plus, model, first_panel) -
       panel_poisson_re_loglik(minus, model, first_panel))/(2*step[j])
  }, numeric(1))
  expect_equal(colSums(component$score[rows, , drop = FALSE]),
    numerical_score, tolerance = 2e-6, ignore_attr = TRUE)
  expect_equal(sum(rowSums(abs(component$score)) > 0), length(unique(panel)))
})

test_that("gamma random-effects Poisson supports panel system alignment", {
  skip_if_not_installed("pglm")
  dat <- panel_poisson_re_test_data(groups = 80L)
  dat$y1[dat$time == 6 & dat$id %% 4 == 0] <- NA
  dat$y2[dat$time == 1 & dat$id %% 5 == 0] <- NA
  first <- panel_poisson_re_fit(y1 ~ x + z, dat)
  second <- panel_poisson_re_fit(y2 ~ x + z, dat)
  higher_cluster <- lapply(list(first, second), function(model) {
    frame <- suest:::.suest_model_frame(model, "pglm::pglm")
    dat$higher[match(rownames(frame), rownames(dat))]
  })
  partial <- suest(first, second, model_names = c("First", "Second"),
    observation_id = c("id", "time"), cluster = higher_cluster)

  expect_equal(partial$n_clusters, length(unique(dat$higher)))
  expect_lt(partial$nobs_overlap, partial$nobs_union)
  expect_true(all(c("First::alpha", "Second::alpha") %in% names(coef(partial))))
  expect_true(all(coef(partial)[c("First::alpha", "Second::alpha")] > 0))
  expect_true(all(is.finite(vcov(partial))))
  expect_gt(max(abs(offdiag_vcov(partial))), 0)

  low <- dat[dat$id <= 40, ]
  high <- dat[dat$id > 40, ]
  disjoint <- suest(
    panel_poisson_re_fit(y1 ~ x + z, low),
    panel_poisson_re_fit(y2 ~ x + z, high),
    model_names = c("Low", "High"), observation_id = c("id", "time")
  )
  expect_equal(disjoint$nobs_overlap, 0L)
  expect_equal(unname(offdiag_vcov(disjoint)), matrix(0, 4L, 4L),
    tolerance = 1e-12)
})

test_that("gamma random-effects Poisson works with marginaleffects", {
  skip_if_not_installed("pglm")
  skip_if_not_installed("marginaleffects")
  dat <- panel_poisson_re_test_data()
  first <- panel_poisson_re_fit(y1 ~ x*woman + z, dat)
  second <- panel_poisson_re_fit(y2 ~ x*woman + z, dat)
  fit <- suest(first, second, model_names = c("First", "Second"),
    observation_id = c("id", "time"))
  predictions <- marginaleffects::avg_predictions(fit, newdata = dat)
  slopes <- marginaleffects::avg_slopes(fit, variables = "x", newdata = dat)
  expected <- vapply(list(first, second), function(model) mean(
    suest:::.suest_predict_values(
      model, dat, "response", "pglm::pglm", "panel_poisson_re"
    )
  ), numeric(1))

  expect_equal(predictions$estimate, expected, tolerance = 1e-8)
  expect_true(all(is.finite(predictions$std.error)))
  expect_true(all(is.finite(slopes$estimate)))
  expect_true(all(is.finite(slopes$std.error)))

  changed <- coef(fit)
  changed["First::alpha"] <- changed["First::alpha"]*1.5
  altered <- marginaleffects::set_coef(fit, changed)
  original_prediction <- marginaleffects::predictions(
    fit, newdata = dat[1:5, ], vcov = FALSE)
  altered_prediction <- marginaleffects::predictions(
    altered, newdata = dat[1:5, ], vcov = FALSE)
  expect_equal(original_prediction$estimate, altered_prediction$estimate,
    tolerance = 1e-12)
})

test_that("gamma random-effects Poisson validation is deliberately narrow", {
  skip_if_not_installed("pglm")
  dat <- panel_poisson_re_test_data(groups = 45L)
  standard <- panel_poisson_re_fit(y1 ~ x + z, dat)
  inverse <- panel_poisson_re_fit(y1 ~ x + z, dat, other = "inv")
  expect_s3_class(suest(standard, standard), "suest_model")
  expect_error(suest(inverse, inverse), "other = \"sd\"", fixed = TRUE)

  noninteger <- dat
  noninteger$y1 <- noninteger$y1 + 0.25
  invalid <- panel_poisson_re_fit(y1 ~ x + z, noninteger)
  expect_error(suest(invalid, invalid), "nonnegative integer counts")
})

test_that("panel Poisson reproduces the returned Stata component benchmark", {
  skip_if_not_installed("pglm")
  skip_if_not_installed("marginaleffects")
  fixture <- readRDS(test_path("fixtures", "panel-poisson-re-stata.rds"))
  dat <- fixture$data
  fit_model <- function(formula) panel_poisson_re_fit(formula, dat)

  transform_model <- function(model) {
    parameters <- coef(model)
    names(parameters)[length(parameters)] <- "alpha"
    covariance <- vcov(model)
    dimnames(covariance) <- list(names(parameters), names(parameters))
    jacobian <- diag(length(parameters))
    jacobian[nrow(jacobian), ncol(jacobian)] <- 1/parameters["alpha"]
    parameters["alpha"] <- log(parameters["alpha"])
    names(parameters)[length(parameters)] <- "lnalpha"
    covariance <- jacobian %*% covariance %*% t(jacobian)
    dimnames(covariance) <- list(names(parameters), names(parameters))
    list(coefficients = parameters, vcov = covariance,
      logLik = as.numeric(logLik(model)))
  }
  transform_system <- function(fit) {
    parameters <- coef(fit)
    alpha <- grep("::alpha$", names(parameters))
    jacobian <- diag(length(parameters))
    jacobian[cbind(alpha, alpha)] <- 1/parameters[alpha]
    parameters[alpha] <- log(parameters[alpha])
    names(parameters)[alpha] <- sub("alpha$", "lnalpha", names(parameters)[alpha])
    covariance <- jacobian %*% vcov(fit) %*% t(jacobian)
    dimnames(covariance) <- list(names(parameters), names(parameters))
    list(coefficients = parameters, vcov = covariance,
      nobs_union = fit$nobs_union, nobs_overlap = fit$nobs_overlap,
      n_clusters = fit$n_clusters)
  }

  first <- fit_model(y1 ~ x + z)
  second <- fit_model(y2 ~ x + z)
  r_models <- lapply(list(first, second), transform_model)
  stata_models <- fixture$balanced[c("model1", "model2")]
  for (i in seq_along(r_models)) {
    expect_equal(r_models[[i]]$coefficients,
      stata_models[[i]]$coefficients, tolerance = 2e-6)
    expect_lt(max(abs(r_models[[i]]$vcov -
      stata_models[[i]]$vcov)), 1e-8)
    expect_equal(r_models[[i]]$logLik,
      stata_models[[i]]$logLik, tolerance = 1e-6)
  }

  balanced <- suest(first, second, model_names = c("Y1", "Y2"),
    observation_id = c("id", "time"))
  predictions <- marginaleffects::avg_predictions(balanced, newdata = dat)
  slopes <- marginaleffects::avg_slopes(
    balanced, variables = "x", newdata = dat)
  expect_equal(predictions$estimate, fixture$balanced$prediction_mean,
    tolerance = 2e-7, ignore_attr = TRUE)
  expect_equal(slopes$estimate, fixture$balanced$average_slope_x,
    tolerance = 1e-7, ignore_attr = TRUE)

  partial_models <- list(
    fit_model(y1_partial ~ x + z), fit_model(y2_partial ~ x + z))
  higher_cluster <- lapply(partial_models, function(model) {
    frame <- suest:::.suest_model_frame(model, "pglm::pglm")
    dat$higher[match(rownames(frame), rownames(dat))]
  })
  partial <- suest(
    partial_models[[1]], partial_models[[2]],
    model_names = c("Y1", "Y2"), observation_id = c("id", "time"),
    cluster = higher_cluster)
  disjoint <- suest(
    fit_model(y1_left ~ x + z), fit_model(y2_right ~ x + z),
    model_names = c("Y1", "Y2"), observation_id = c("id", "time"))
  r_systems <- list(
    balanced = transform_system(balanced),
    partial_higher = transform_system(partial),
    disjoint = transform_system(disjoint))

  for (case in names(r_systems)) {
    expect_equal(r_systems[[case]]$coefficients,
      fixture[[case]]$system$coefficients, tolerance = 2e-6)
    expect_identical(r_systems[[case]]$nobs_union,
      fixture[[case]]$system$nobs_union)
    expect_identical(r_systems[[case]]$nobs_overlap,
      fixture[[case]]$system$nobs_overlap)
    expect_identical(r_systems[[case]]$n_clusters,
      fixture[[case]]$system$n_clusters)
  }

  # suest2 1.0.0 deliberately repeats the full ancillary cluster score on
  # every observation, inflating lnalpha rows/columns. Its higher-cluster route
  # also integrates at the higher cluster instead of the retained panel.
  fixed <- !grepl("lnalpha", names(r_systems$balanced$coefficients))
  expect_lt(max(abs(r_systems$balanced$vcov[fixed, fixed] -
    fixture$balanced$system$vcov[fixed, fixed])), 4e-4)
  expect_gt(max(abs(r_systems$balanced$vcov -
    fixture$balanced$system$vcov)), 1)
  expect_gt(max(abs(r_systems$partial_higher$vcov[fixed, fixed] -
    fixture$partial_higher$system$vcov[fixed, fixed])), 0.005)
  expect_lt(max(abs(r_systems$disjoint$vcov[fixed, fixed] -
    fixture$disjoint$system$vcov[fixed, fixed])), 7e-4)
  expect_equal(fixture$disjoint$system$vcov[1:4, 5:8],
    matrix(0, 4, 4), tolerance = 1e-12, ignore_attr = TRUE)
})
