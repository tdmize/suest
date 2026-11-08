flexsurv_test_data <- function(n = 400L) {
  set.seed(4101)
  x <- rnorm(n)
  z <- rbinom(n, 1L, 0.5)
  t <- rweibull(n, 1.3, exp(1 + 0.3*x - 0.2*z))
  censor <- runif(n, 0, 8)
  entry <- ifelse(runif(n) < 0.3, runif(n, 0, 0.5), 0)
  keep <- pmin(t, censor) > entry
  data.frame(id = seq_len(n), time = pmin(t, censor), dead = as.integer(t <= censor),
    entry = entry, x = x, z = z)[keep, ]
}

flexsurv_test_models <- function(d) list(
  gompertz = flexsurv::flexsurvreg(survival::Surv(time, dead) ~ x + z, data = d,
    dist = "gompertz"),
  gompertz_anc = flexsurv::flexsurvreg(survival::Surv(time, dead) ~ x + z, data = d,
    dist = "gompertz", anc = list(shape = ~ z)),
  gengamma = flexsurv::flexsurvreg(survival::Surv(time, dead) ~ x + z, data = d,
    dist = "gengamma"),
  gengamma_anc = flexsurv::flexsurvreg(survival::Surv(time, dead) ~ x + z, data = d,
    dist = "gengamma", anc = list(sigma = ~ z, Q = ~ x)),
  gompertz_entry = flexsurv::flexsurvreg(survival::Surv(entry, time, dead) ~ x + z,
    data = d, dist = "gompertz"))

# flexsurv's own log likelihood at theta: refit with every parameter fixed.
flexsurv_loglik_at <- function(model, theta, d) {
  inits <- theta
  for (j in seq_along(model$dlist$pars))
    inits[model$basepars[j]] <- model$dlist$inv.transforms[[j]](theta[model$basepars[j]])
  call <- model$call
  call$inits <- unname(inits)
  call$fixedpars <- TRUE
  call$data <- d
  eval(call)$loglik
}

test_that("Gompertz and generalized gamma scores and information are exact", {
  skip_if_not_installed("flexsurv")
  skip_if_not_installed("numDeriv")
  d <- flexsurv_test_data()
  for (name in names(flexsurv_test_models(d))) {
    model <- flexsurv_test_models(d)[[name]]
    adapter <- suest:::.suest_model_adapter(model)
    expect_identical(adapter$type, "flexsurv", label = name)
    component <- suest:::.suest_model_components(model, "flexsurv", adapter$engine)
    theta <- component$parameters
    expect_identical(names(theta), names(model$coefficients), label = name)
    X <- suest:::.suest_flexsurv_design(model)
    loglik <- function(b) suest:::.suest_flexsurv_loglik(model, b, X, model$data$Y)
    expect_equal(unname(loglik(theta)), unname(model$logliki), tolerance = 1e-10, label = name)
    moved <- theta + 0.03
    expect_equal(sum(loglik(moved)), flexsurv_loglik_at(model, moved, d), tolerance = 1e-10,
      label = name)
    expect_equal(unname(component$score), unname(numDeriv::jacobian(loglik, theta)),
      tolerance = 1e-7, label = name)
    information <- -numDeriv::hessian(function(b) sum(loglik(b)), theta)
    expect_equal(unname(component$bread), unname(nrow(X[[1L]])*solve(information)),
      tolerance = 1e-5, label = name)
    expect_equal(unname(component$bread/nrow(X[[1L]])), unname(model$cov), tolerance = 1e-3,
      label = name)
  }
})

test_that("Gompertz scores are exact near zero shape", {
  skip_if_not_installed("flexsurv")
  skip_if_not_installed("numDeriv")
  d <- flexsurv_test_data()
  model <- flexsurv_test_models(d)$gompertz_entry
  X <- suest:::.suest_flexsurv_design(model)
  for (shape in c(-0.4, -1e-5, 0, 1e-6, 0.3)) {
    theta <- model$coefficients
    theta["shape"] <- shape
    analytic <- suest:::.suest_gompertz_scores(model, theta, X, model$data$Y)
    numeric <- numDeriv::jacobian(function(b)
      suest:::.suest_flexsurv_loglik(model, b, X, model$data$Y), theta)
    expect_equal(unname(analytic), numeric, tolerance = 1e-7, label = format(shape))
  }
})

test_that("flexsurv models combine with survreg and predict median survival times", {
  skip_if_not_installed("flexsurv")
  skip_if_not_installed("marginaleffects")
  d <- flexsurv_test_data()
  models <- flexsurv_test_models(d)
  weibull <- survival::survreg(survival::Surv(time, dead) ~ x + z, data = d, dist = "weibull")
  fit <- suest(models$gompertz, models$gengamma_anc, weibull,
    model_names = c("Gompertz", "GenGamma", "Weibull"), observation_id = "id")
  expect_identical(unname(fit$model_types), c("flexsurv", "flexsurv", "survreg"))
  expect_equal(fit$nobs_union, nrow(d))
  for (name in c("gompertz", "gengamma_anc")) {
    model <- models[[name]]
    expected <- summary(model, newdata = d[1:5, ], type = "median", tidy = TRUE)$est
    expect_equal(suest:::.suest_predict_values(model, d[1:5, ], "response",
      "flexsurv::flexsurvreg", "flexsurv"), expected, tolerance = 1e-10, label = name)
  }
  expect_equal(suest:::.suest_predict_values(models$gompertz, d[1:5, ], "link",
    "flexsurv::flexsurvreg", "flexsurv"),
    unname(log(predict(models$gompertz, newdata = d[1:5, ], type = "link")$.pred_link)),
    tolerance = 1e-10)

  averages <- marginaleffects::avg_predictions(fit, newdata = d)
  expect_equal(averages$estimate[1:2], c(
    mean(summary(models$gompertz, newdata = d, type = "median", tidy = TRUE)$est),
    mean(summary(models$gengamma_anc, newdata = d, type = "median", tidy = TRUE)$est)),
    tolerance = 1e-8)
  expect_true(all(is.finite(averages$std.error)))

  component <- suest:::.suest_model_components(models$gompertz, "flexsurv",
    "flexsurv::flexsurvreg")
  n <- nrow(d)
  index <- fit$index[["Gompertz"]]
  expected <- component$bread %*% crossprod(component$score) %*% component$bread/n^2*n/(n - 1)
  expect_equal(unname(vcov(fit)[index, index]), unname(expected), tolerance = 1e-10)
})

test_that("unsupported flexsurv models are refused with a reason", {
  skip_if_not_installed("flexsurv")
  d <- flexsurv_test_data()
  weibull <- flexsurv::flexsurvreg(survival::Surv(time, dead) ~ x, data = d, dist = "weibull")
  expect_error(suest(weibull, weibull), "survival::survreg")
  weighted <- flexsurv::flexsurvreg(survival::Surv(time, dead) ~ x, data = d,
    dist = "gompertz", weights = rep(2, nrow(d)))
  expect_error(suest(weighted, weighted), "weights")
  d$upper <- ifelse(d$dead == 1, d$time, Inf)
  interval <- flexsurv::flexsurvreg(survival::Surv(time, upper, type = "interval2") ~ x,
    data = d, dist = "gompertz")
  expect_error(suest(interval, interval), "right-censored")
})

test_that("the generalized gamma density stays accurate when Q is near 0", {
  skip_if_not_installed("flexsurv")
  density <- suest:::.suest_gengamma_log_density
  # Exact values (400-bit arithmetic) at t = 2, mu = 1, sigma = 0.8.
  expect_equal(density(2, 1, 0.8, c(-1e-4, 3e-3, 0.05)),
    c(-1.4625045513138395, -1.4624761523921921, -1.4622439281217321), tolerance = 1e-14)
  expect_equal(density(c(0.05, 2, 30), 1, 0.8, 0), dlnorm(c(0.05, 2, 30), 1, 0.8, log = TRUE),
    tolerance = 1e-14)
  expect_equal(density(c(0.05, 2, 30), 1, 0.8, 1e-9), dlnorm(c(0.05, 2, 30), 1, 0.8, log = TRUE),
    tolerance = 1e-8)
  grid <- expand.grid(t = c(0.05, 1, 2, 30), Q = c(-2, -0.5, -0.11, -0.099, -0.03, 0.03, 0.099,
    0.101, 0.4, 1, 2.5))
  expect_equal(density(grid$t, 1, 0.8, grid$Q), flexsurv::dgengamma(grid$t, 1, 0.8, grid$Q,
    log = TRUE), tolerance = 1e-11)
})

test_that("generalized gamma scores are smooth when Q depends on a covariate and crosses 0", {
  skip_if_not_installed("flexsurv")
  skip_if_not_installed("numDeriv")
  set.seed(4102)
  n <- 300L
  d <- data.frame(id = seq_len(n), x = rnorm(n), z = rnorm(n))
  t <- flexsurv::rgengamma(n, 1 + 0.4*d$x, exp(-0.3 + 0.2*d$z), 0.3 + 0.3*d$x)
  censor <- runif(n, 0.5, 15)
  d$time <- pmin(t, censor); d$dead <- as.integer(t <= censor)
  model <- flexsurv::flexsurvreg(survival::Surv(time, dead) ~ x + z, data = d,
    dist = "gengamma", anc = list(sigma = ~ z, Q = ~ x))
  # Newton steps from flexsurv's estimates converge to a zero score only when
  # the scores are free of rounding noise.
  for (step in 1:3) {
    component <- suest:::.suest_model_components(model, "flexsurv", "flexsurv::flexsurvreg")
    total <- colSums(component$score)
    model$suest_parameters <- component$parameters + drop(component$bread %*% total)/n
  }
  expect_lt(max(abs(total)), 1e-6)
})
