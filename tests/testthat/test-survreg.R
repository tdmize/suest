survreg_test_data <- function(n = 700L) {
  set.seed(2951)
  x <- rnorm(n)
  z <- rnorm(n)
  latent <- 0.4 + 0.7*x - 0.35*z + rnorm(n, sd = 0.9)
  time <- exp(0.3 + 0.25*x - 0.15*z + rnorm(n, sd = 0.7))

  data.frame(
    id = seq_len(n),
    x = x,
    z = z,
    y = pmax(0, latent),
    event = latent > 0,
    lower = floor(2*latent)/2,
    upper = floor(2*latent)/2 + 0.5,
    time = time,
    stratum = factor(rep(c("A", "B"), length.out = n))
  )
}

test_that("Gaussian censored regressions include the ancillary scale", {
  skip_if_not_installed("survival")
  dat <- survreg_test_data()
  base <- survival::survreg(
    survival::Surv(y, event, type = "left") ~ x,
    data = dat,
    dist = "gaussian",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  adjusted <- survival::survreg(
    survival::Surv(y, event, type = "left") ~ x + z,
    data = dat,
    dist = "gaussian",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  fit <- suest(base, adjusted, model_names = c("Base", "Adjusted"))

  expect_equal(unname(fit$model_types), rep("survreg", 2L))
  expect_equal(
    unname(fit$model_engines),
    rep("survival::survreg", 2L)
  )
  expect_true(all(c(
    "Base::Log(scale)",
    "Adjusted::Log(scale)"
  ) %in% names(coef(fit))))
  expect_equal(
    unname(vcov(fit)[fit$index[[1L]], fit$index[[1L]], drop = FALSE]),
    unname(robust_vcov_direct(base)),
    tolerance = 1e-8
  )
  expect_equal(
    unname(vcov(fit)[fit$index[[2L]], fit$index[[2L]], drop = FALSE]),
    unname(robust_vcov_direct(adjusted)),
    tolerance = 1e-8
  )

  effects <- marginaleffects::avg_comparisons(
    fit,
    variables = "x",
    newdata = dat
  )
  expect_equal(nrow(effects), 2L)
  expect_true(all(is.finite(effects$estimate)))
  expect_true(all(is.finite(effects$std.error)))
})

test_that("Gaussian interval regression works with marginaleffects", {
  skip_if_not_installed("survival")
  dat <- survreg_test_data()
  base <- survival::survreg(
    survival::Surv(lower, upper, type = "interval2") ~ x,
    data = dat,
    dist = "gaussian",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  adjusted <- survival::survreg(
    survival::Surv(lower, upper, type = "interval2") ~ x + z,
    data = dat,
    dist = "gaussian",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  fit <- suest(base, adjusted)
  effects <- marginaleffects::avg_comparisons(
    fit,
    variables = "x",
    newdata = dat
  )

  expect_equal(unname(fit$model_types), rep("survreg", 2L))
  expect_equal(
    unname(vcov(fit)[fit$index[[1L]], fit$index[[1L]], drop = FALSE]),
    unname(robust_vcov_direct(base)),
    tolerance = 1e-8
  )
  expect_equal(nrow(effects), 2L)
  expect_true(all(is.finite(effects$std.error)))
})

test_that("parametric survival models work with marginaleffects", {
  skip_if_not_installed("survival")
  dat <- survreg_test_data()
  weibull <- survival::survreg(
    survival::Surv(time) ~ x,
    data = dat,
    dist = "weibull",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  lognormal <- survival::survreg(
    survival::Surv(time) ~ x + z,
    data = dat,
    dist = "lognormal",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  fit <- suest(weibull, lognormal)
  predictions <- marginaleffects::avg_predictions(fit, newdata = dat)
  effects <- marginaleffects::avg_comparisons(
    fit,
    variables = "x",
    newdata = dat
  )

  expect_equal(nrow(predictions), 2L)
  expect_equal(nrow(effects), 2L)
  expect_true(all(is.finite(predictions$std.error)))
  expect_true(all(is.finite(effects$std.error)))
})

test_that("fixed and estimated survreg scales can be combined", {
  skip_if_not_installed("survival")
  dat <- survreg_test_data()
  fixed <- survival::survreg(
    survival::Surv(y, event, type = "left") ~ x,
    data = dat,
    dist = "gaussian",
    scale = 0.9,
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  estimated <- survival::survreg(
    survival::Surv(y, event, type = "left") ~ x + z,
    data = dat,
    dist = "gaussian",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  fit <- suest(fixed, estimated)

  expect_false(any(grepl("^fixed::Log\\(scale\\)$", names(coef(fit)))))
  expect_true("estimated::Log(scale)" %in% names(coef(fit)))
  expect_true(all(is.finite(vcov(fit))))
})

test_that("stratum-specific survreg scales are rejected", {
  skip_if_not_installed("survival")
  dat <- survreg_test_data()
  strata <- survival::strata
  stratified <- survival::survreg(
    survival::Surv(time) ~ x + strata(stratum),
    data = dat,
    dist = "weibull",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  ordinary <- survival::survreg(
    survival::Surv(time) ~ x,
    data = dat,
    dist = "weibull",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )

  expect_error(
    suest(stratified, ordinary),
    "stratum-specific scales"
  )
})
