test_that("unsupported model diagnostics identify the offending model", {
  fit <- lm(mpg ~ wt, data = mtcars)
  unknown <- structure(list(), class = "unknown_regression")
  expect_error(
    suest(fit, unknown, model_names = c("Linear", "Unknown")),
    "Each model must be a supported.*Unknown.*unknown_regression"
  )
})

test_that("unsupported glmmTMB family diagnostics give the family and link", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  set.seed(293)
  d <- data.frame(id = factor(rep(seq_len(18), each = 6)), x = rnorm(108))
  d$y <- exp(.3 + .2*d$x + rep(rnorm(18, sd = .3), each = 6) +
    rnorm(108, sd = .15))
  fit <- glmmTMB::glmmTMB(y ~ x + (1 | id), data = d,
    family = stats::Gamma(link = "log"))
  expect_error(suest(fit, fit), "family 'Gamma'.*link 'log'")
})
