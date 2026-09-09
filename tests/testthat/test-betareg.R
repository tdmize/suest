betareg_test_data <- function(n = 900L) {
  set.seed(4812)
  x <- rnorm(n)
  z <- rnorm(n)
  precision <- exp(3.4 - 0.25*z)
  mean <- plogis(-0.3 + 0.55*x - 0.2*z)

  data.frame(
    id = seq_len(n),
    x = x,
    z = z,
    y = rbeta(n, mean*precision, (1 - mean)*precision)
  )
}

test_that("beta-regression mean links work with marginaleffects", {
  skip_if_not_installed("betareg")
  dat <- betareg_test_data()

  for (link in c("logit", "probit", "cloglog", "loglog")) {
    base <- betareg::betareg(
      y ~ x,
      data = dat,
      link = link,
      model = TRUE,
      x = TRUE,
      y = TRUE
    )
    adjusted <- betareg::betareg(
      y ~ x + z,
      data = dat,
      link = link,
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

    expect_equal(unname(fit$model_types), rep("betareg", 2L))
    expect_equal(
      unname(fit$model_engines),
      rep("betareg::betareg", 2L)
    )
    expect_equal(nrow(effects), 2L)
    expect_true(all(is.finite(effects$estimate)))
    expect_true(all(is.finite(effects$std.error)))
  }
})

test_that("beta-regression precision parameters enter the joint covariance", {
  skip_if_not_installed("betareg")
  dat <- betareg_test_data()
  constant_precision <- betareg::betareg(
    y ~ x + z,
    data = dat,
    link = "logit",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  variable_precision <- betareg::betareg(
    y ~ x + z | z,
    data = dat,
    link = "probit",
    link.phi = "log",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  fit <- suest(
    constant_precision,
    variable_precision,
    model_names = c("Constant", "Variable")
  )

  expect_true("Constant::(phi)" %in% names(coef(fit)))
  expect_true("Variable::(phi)_z" %in% names(coef(fit)))
  expect_equal(
    unname(vcov(fit)[fit$index[[1L]], fit$index[[1L]], drop = FALSE]),
    unname(robust_vcov_direct(constant_precision)),
    tolerance = 1e-8
  )
  expect_equal(
    unname(vcov(fit)[fit$index[[2L]], fit$index[[2L]], drop = FALSE]),
    unname(robust_vcov_direct(variable_precision)),
    tolerance = 1e-8
  )
})

test_that("beta regression can be combined with other scalar models", {
  skip_if_not_installed("betareg")
  dat <- betareg_test_data()
  beta <- betareg::betareg(
    y ~ x + z,
    data = dat,
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  fractional <- glm(y ~ x + z, family = quasibinomial(), data = dat)
  fit <- suest(beta, fractional, model_names = c("Beta", "Fractional"))
  effects <- marginaleffects::avg_comparisons(
    fit,
    variables = "x",
    newdata = dat
  )

  expect_equal(nrow(effects), 2L)
  expect_true(all(is.finite(effects$std.error)))
})
