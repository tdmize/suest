glm_family_data <- function(n = 900L) {
  set.seed(8135)
  x <- rnorm(n)
  z <- rnorm(n)
  eta <- -0.2 + 0.45*x - 0.25*z
  probability <- 1 - exp(-exp(eta))
  fractional <- plogis(eta + rnorm(n, sd = 0.8))
  data.frame(
    x = x,
    z = z,
    y_gaussian = 1 + eta + rnorm(n),
    y_gamma = rgamma(n, shape = 3, scale = exp(0.4 + 0.15*eta)/3),
    y_binary = rbinom(n, 1, probability),
    y_fractional = fractional
  )
}

loglog_link <- function() {
  structure(
    list(
      linkfun = function(mu) -log(-log(mu)),
      linkinv = function(eta) exp(-exp(-eta)),
      mu.eta = function(eta) exp(-eta - exp(-eta)),
      valideta = function(eta) TRUE,
      name = "loglog"
    ),
    class = "link-glm"
  )
}

test_that("additional GLM families work with marginaleffects", {
  dat <- glm_family_data()
  gaussian <- glm(
    y_gaussian ~ x + z,
    family = gaussian("identity"),
    data = dat
  )
  gamma <- glm(y_gamma ~ x + z, family = Gamma("log"), data = dat)
  cloglog <- glm(
    y_binary ~ x + z,
    family = binomial("cloglog"),
    data = dat
  )
  fractional <- glm(
    y_fractional ~ x + z,
    family = quasibinomial("logit"),
    data = dat
  )
  fit <- suest(
    gaussian,
    gamma,
    cloglog,
    fractional,
    model_names = c("Gaussian", "Gamma", "Cloglog", "Fractional")
  )
  effects <- marginaleffects::avg_comparisons(
    fit,
    variables = "x",
    newdata = dat
  )

  expect_equal(
    unname(fit$model_types),
    c("glm", "glm", "cloglog", "fraclogit")
  )
  expect_equal(nrow(effects), 4L)
  expect_true(all(is.finite(fit$vcov)))
  expect_true(all(is.finite(effects$estimate)))
  expect_true(all(is.finite(effects$std.error)))
})

test_that("same-family additional GLMs preserve robust covariance", {
  dat <- glm_family_data()
  specifications <- list(
    gaussian = list(outcome = "y_gaussian", family = gaussian("identity")),
    gamma = list(outcome = "y_gamma", family = Gamma("log")),
    cloglog = list(outcome = "y_binary", family = binomial("cloglog")),
    fraclogit = list(
      outcome = "y_fractional",
      family = quasibinomial("logit")
    ),
    fracprobit = list(
      outcome = "y_fractional",
      family = quasibinomial("probit")
    ),
    fraccloglog = list(
      outcome = "y_fractional",
      family = quasibinomial("cloglog")
    ),
    fracloglog = list(
      outcome = "y_fractional",
      family = quasibinomial(loglog_link())
    )
  )

  for (specification in specifications) {
    base <- glm(
      reformulate("x", response = specification$outcome),
      family = specification$family,
      data = dat
    )
    adjusted <- glm(
      reformulate(c("x", "z"), response = specification$outcome),
      family = specification$family,
      data = dat
    )
    fit <- suest(base, adjusted)

    expect_equal(
      unname(fit$vcov[fit$index[[1L]], fit$index[[1L]], drop = FALSE]),
      unname(robust_vcov_direct(base)),
      tolerance = 1e-8
    )
    expect_equal(
      unname(fit$vcov[fit$index[[2L]], fit$index[[2L]], drop = FALSE]),
      unname(robust_vcov_direct(adjusted)),
      tolerance = 1e-8
    )
  }
})

test_that("fractional log-log works with marginaleffects", {
  dat <- glm_family_data()
  base <- glm(
    y_fractional ~ x,
    family = quasibinomial(loglog_link()),
    data = dat
  )
  adjusted <- glm(
    y_fractional ~ x + z,
    family = quasibinomial(loglog_link()),
    data = dat
  )
  fit <- suest(base, adjusted)
  effects <- marginaleffects::avg_comparisons(
    fit,
    variables = "x",
    newdata = dat
  )

  expect_equal(unname(fit$model_types), rep("fracloglog", 2L))
  expect_equal(nrow(effects), 2L)
  expect_true(all(is.finite(effects$std.error)))
})

test_that("fractional outcomes must remain in the unit interval", {
  dat <- glm_family_data()
  model1 <- glm(y_fractional ~ x, family = quasibinomial(), data = dat)
  model2 <- glm(y_fractional ~ x + z, family = quasibinomial(), data = dat)
  model1$model[[1L]][1L] <- 2

  expect_error(
    suest(model1, model2),
    "numeric values from zero to one"
  )
})
