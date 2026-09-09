hetprob_test_data <- function(n = 1000L) {
  set.seed(7406)
  x <- rnorm(n)
  z <- rnorm(n)
  w <- rnorm(n)
  mean <- 0.2 + 0.7*x - 0.25*w
  scale <- exp(0.35*z)

  data.frame(
    id = seq_len(n),
    y = rbinom(n, 1, pnorm(mean/scale)),
    x = x,
    z = z,
    w = w
  )
}

test_that("heteroskedastic probit preserves the full robust covariance", {
  skip_if_not_installed("Rchoice")
  dat <- hetprob_test_data()
  base <- Rchoice::hetprob(y ~ x | z, data = dat, link = "probit")
  adjusted <- Rchoice::hetprob(
    y ~ x + w | z + w,
    data = dat,
    link = "probit"
  )
  fit <- suest(base, adjusted, model_names = c("Base", "Adjusted"))

  expect_equal(unname(fit$model_types), rep("hetprobit", 2L))
  expect_equal(unname(fit$model_engines), rep("Rchoice::hetprob", 2L))
  expect_true(any(grepl("het\\.", names(coef(fit)))))
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
  expect_true(all(is.finite(effects$std.error)))
})

test_that("heteroskedastic binary scores equal likelihood derivatives", {
  skip_if_not_installed("Rchoice")
  dat <- hetprob_test_data(600L)
  model <- Rchoice::hetprob(y ~ x + w | z, data = dat, link = "probit")
  parameters <- stats::coef(model)
  numerical <- matrix(NA_real_, nrow(dat), length(parameters))
  step <- 1e-6

  loglikelihood <- function(candidate) {
    updated <- suest:::.suest_set_parameters(
      model,
      candidate,
      "hetprobit",
      "Rchoice::hetprob"
    )
    probability <- suest:::.suest_predict_values(
      updated,
      dat,
      "response",
      "Rchoice::hetprob",
      "hetprobit"
    )
    dat$y*log(probability) + (1 - dat$y)*log1p(-probability)
  }

  for (j in seq_along(parameters)) {
    plus <- minus <- parameters
    plus[j] <- plus[j] + step
    minus[j] <- minus[j] - step
    numerical[, j] <- (loglikelihood(plus) - loglikelihood(minus))/(2*step)
  }

  expect_equal(
    unname(sandwich::estfun(model)),
    unname(numerical),
    tolerance = 1e-5
  )
})

test_that("the adapter repairs Rchoice's omitted-link prediction bug", {
  skip_if_not_installed("Rchoice")
  dat <- hetprob_test_data()
  default <- Rchoice::hetprob(y ~ x | z, data = dat)
  explicit <- Rchoice::hetprob(y ~ x + w | z, data = dat, link = "probit")

  expect_error(
    stats::predict(default, newdata = dat, type = "pr"),
    "length 1 vector"
  )

  fit <- suest(default, explicit)
  effects <- marginaleffects::avg_comparisons(
    fit,
    variables = "x",
    newdata = dat
  )
  expect_equal(unname(fit$model_types), rep("hetprobit", 2L))
  expect_equal(nrow(effects), 2L)
  expect_true(all(is.finite(effects$std.error)))
})

test_that("heteroskedastic logit is available as an R extension", {
  skip_if_not_installed("Rchoice")
  dat <- hetprob_test_data()
  base <- Rchoice::hetprob(y ~ x | z, data = dat, link = "logit")
  adjusted <- Rchoice::hetprob(y ~ x + w | z, data = dat, link = "logit")
  fit <- suest(base, adjusted)
  effects <- marginaleffects::avg_comparisons(
    fit,
    variables = "x",
    newdata = dat
  )

  expect_equal(unname(fit$model_types), rep("hetlogit", 2L))
  expect_true(all(is.finite(effects$std.error)))
})
