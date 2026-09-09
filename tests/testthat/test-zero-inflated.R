zeroinfl_test_data <- function(n = 1000L) {
  set.seed(7642)
  x <- rnorm(n)
  z <- rnorm(n)
  structural_zero <- rbinom(n, 1, plogis(-1 + 0.45*z))
  count <- rnbinom(n, mu = exp(0.4 + 0.35*x - 0.2*z), size = 2.8)

  data.frame(
    id = seq_len(n),
    x = x,
    z = z,
    y = ifelse(structural_zero == 1, 0, count)
  )
}

test_that("zero-inflated Poisson preserves its robust covariance", {
  skip_if_not_installed("pscl")
  dat <- zeroinfl_test_data()
  base <- pscl::zeroinfl(
    y ~ x | z,
    data = dat,
    dist = "poisson",
    link = "logit",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  adjusted <- pscl::zeroinfl(
    y ~ x + z | x + z,
    data = dat,
    dist = "poisson",
    link = "probit",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  fit <- suest(base, adjusted)

  expect_equal(unname(fit$model_types), rep("zip", 2L))
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
})

test_that("ZINB includes dispersion and matches numerical likelihood scores", {
  skip_if_not_installed("pscl")
  dat <- zeroinfl_test_data()
  model <- pscl::zeroinfl(
    y ~ x + z | z,
    data = dat,
    dist = "negbin",
    link = "logit",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  components <- .suest_zinb_components(model)
  parameters <- components$parameters

  expect_true("ln_theta" %in% names(parameters))
  expect_equal(
    components$score[, names(stats::coef(model)), drop = FALSE],
    sandwich::estfun(model),
    tolerance = 1e-10
  )
  expect_equal(
    unname(components$bread/nrow(components$score)),
    unname(solve(-model$optim$hessian)),
    tolerance = 1e-10
  )

  observation_loglik <- function(par, row) {
    count_names <- names(model$coefficients$count)
    zero_names <- names(model$coefficients$zero)
    count_index <- seq_along(count_names)
    zero_index <- length(count_names) + seq_along(zero_names)
    theta_index <- length(par)
    mu <- exp(sum(model$x$count[row, ]*par[count_index]))
    zero_eta <- sum(model$x$zero[row, ]*par[zero_index])
    probability_zero <- model$linkinv(zero_eta)
    theta <- exp(par[theta_index])
    y <- model$y[row]
    likelihood <- if (y == 0) {
      probability_zero +
        (1 - probability_zero)*stats::dnbinom(0, size = theta, mu = mu)
    } else {
      (1 - probability_zero)*stats::dnbinom(y, size = theta, mu = mu)
    }
    model$weights[row]*log(likelihood)
  }

  step <- 1e-6
  rows <- c(which(model$y == 0)[1L], which(model$y > 0)[1L])
  for (row in rows) {
    numerical <- vapply(seq_along(parameters), function(j) {
      upper <- lower <- parameters
      upper[j] <- upper[j] + step
      lower[j] <- lower[j] - step
      (observation_loglik(upper, row) - observation_loglik(lower, row)) /
        (2*step)
    }, numeric(1))

    expect_equal(
      unname(components$score[row, ]),
      unname(numerical),
      tolerance = 1e-6
    )
  }
})

test_that("ZIP and ZINB work together with marginaleffects", {
  skip_if_not_installed("pscl")
  dat <- zeroinfl_test_data()
  zip <- pscl::zeroinfl(
    y ~ x + z | z,
    data = dat,
    dist = "poisson",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  zinb <- pscl::zeroinfl(
    y ~ x + z | z,
    data = dat,
    dist = "negbin",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  fit <- suest(zip, zinb, model_names = c("ZIP", "ZINB"))
  effects <- marginaleffects::avg_comparisons(
    fit,
    variables = "x",
    newdata = dat
  )

  expect_true("ZINB::ln_theta" %in% names(coef(fit)))
  expect_equal(nrow(effects), 2L)
  expect_true(all(is.finite(effects$estimate)))
  expect_true(all(is.finite(effects$std.error)))
})
