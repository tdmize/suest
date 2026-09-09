censreg_test_data <- function(n = 900L) {
  set.seed(6814)
  x <- rnorm(n)
  z <- rnorm(n)
  grp <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  eta <- 0.35 + 0.75*x - 0.3*z + c(A = 0, B = 0.25, C = -0.2)[grp]
  latent <- eta + rnorm(n, sd = 0.85)

  data.frame(
    id = seq_len(n),
    x = x,
    z = z,
    grp = grp,
    latent = latent,
    left_y = pmax(0, latent),
    right_y = pmin(1, latent)
  )
}

test_that("direct Tobit models preserve robust covariance and scale", {
  skip_if_not_installed("censReg")
  dat <- censreg_test_data()
  base <- censReg::censReg(left_y ~ x + grp, left = 0, data = dat)
  adjusted <- censReg::censReg(
    left_y ~ x + z + grp,
    left = 0,
    data = dat
  )
  fit <- suest(base, adjusted, model_names = c("Base", "Adjusted"))

  expect_equal(unname(fit$model_types), rep("censreg", 2L))
  expect_equal(unname(fit$model_engines), rep("censReg::censReg", 2L))
  expect_true(all(c("Base::logSigma", "Adjusted::logSigma") %in%
    names(coef(fit))))
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

test_that("censReg Tobit equations agree with Gaussian survreg", {
  skip_if_not_installed("censReg")
  skip_if_not_installed("survival")
  dat <- censreg_test_data()
  direct <- censReg::censReg(left_y ~ x + z, left = 0, data = dat)
  survival <- survival::survreg(
    survival::Surv(left_y, latent > 0, type = "left") ~ x + z,
    data = dat,
    dist = "gaussian",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )

  expect_equal(
    unname(stats::coef(direct)),
    unname(c(stats::coef(survival), log(survival$scale))),
    tolerance = 1e-7
  )
  expect_equal(
    unname(sandwich::estfun(direct)),
    unname(sandwich::estfun(survival)),
    tolerance = 1e-6
  )
  expect_equal(
    unname(sandwich::bread(direct)),
    unname(sandwich::bread(survival)),
    tolerance = 1e-6
  )

  fit <- suest(direct, survival)
  expect_equal(unname(fit$model_types), c("censreg", "survreg"))
  expect_true(any(abs(offdiag_vcov(fit)) > 1e-12))
  expect_true(all(is.finite(vcov(fit))))
})

test_that("right-censored Tobit models support marginaleffects", {
  skip_if_not_installed("censReg")
  dat <- censreg_test_data()
  base <- censReg::censReg(
    right_y ~ x,
    left = -Inf,
    right = 1,
    data = dat
  )
  adjusted <- censReg::censReg(
    right_y ~ x + z,
    left = -Inf,
    right = 1,
    data = dat
  )
  fit <- suest(base, adjusted)

  effects <- marginaleffects::avg_comparisons(
    fit,
    variables = "x",
    newdata = dat
  )
  expect_equal(nrow(effects), 2L)
  expect_equal(
    effects$estimate,
    unname(c(stats::coef(base)["x"], stats::coef(adjusted)["x"])),
    tolerance = 1e-5
  )
  expect_true(all(is.finite(effects$std.error)))
})
