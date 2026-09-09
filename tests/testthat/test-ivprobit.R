ivprobit_test_data <- function(n = 220L) {
  set.seed(7108)
  x <- rnorm(n)
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  errors <- mvtnorm::rmvnorm(
    n,
    sigma = matrix(c(1, 0.35, 0.35, 1), nrow = 2L)
  )
  endogenous <- 0.25*x + 0.70*z1 + 0.20*z2 + errors[, 2L]
  y <- as.integer(
    0.15 + 0.50*x + 0.65*endogenous + errors[, 1L] > 0
  )
  data.frame(y, x, z1, z2, endogenous)
}

test_that("Rchoice IV probit systems retain every likelihood equation", {
  skip_if_not_installed("Rchoice")
  skip_if_not_installed("mvtnorm")
  dat <- ivprobit_test_data()
  first <- Rchoice::ivpml(
    y ~ x + endogenous | x + z1,
    data = dat, messages = FALSE, printLevel = 0
  )
  second <- Rchoice::ivpml(
    y ~ x + endogenous | x + z1 + z2,
    data = dat, messages = FALSE, printLevel = 0
  )
  fit <- suest(first, second, model_names = c("Base", "Adjusted"))

  expect_equal(fit$nobs_overlap, nrow(dat))
  expect_true(all(c(
    "Base::lnsigma", "Base::athrho",
    "Adjusted::lnsigma", "Adjusted::athrho"
  ) %in% names(coef(fit))))
  expect_true(all(is.finite(vcov(fit))))

  component <- suest:::.suest_ivprobit_components(first)
  expect_equal(component$score, sandwich::estfun(first), ignore_attr = TRUE)
  expect_equal(component$bread, sandwich::bread(first), ignore_attr = TRUE)
  expect_identical(tail(names(component$parameters), 2L), c(
    "lnsigma", "athrho"
  ))

  predictions <- marginaleffects::predictions(fit, newdata = dat[1:3, ])
  structural_probability <- function(model, data) {
    X <- model.matrix(~x + endogenous, data)
    unname(pnorm(drop(X %*% coef(model)[seq_len(ncol(X))])))
  }
  expected <- c(structural_probability(first, dat[1:3, ]),
                structural_probability(second, dat[1:3, ]))
  expect_equal(predictions$estimate, expected, tolerance = 1e-8)
  expected_se <- unlist(lapply(c("Base", "Adjusted"), function(label) {
    X <- model.matrix(~x + endogenous, dat[1:3, ])
    index <- which(startsWith(names(coef(fit)), paste0(label, "::")))
    index <- index[seq_len(ncol(X))]
    eta <- drop(X %*% coef(fit)[index])
    dnorm(eta) * sqrt(rowSums((X %*% vcov(fit)[index, index]) * X))
  }), use.names = FALSE)
  expect_equal(predictions$std.error, expected_se, tolerance = 1e-6)

  # Instruments affect estimation, not the structural prediction at fixed X.
  changed <- dat[1:3, ]
  changed$z1 <- changed$z1 + 10
  changed$z2 <- changed$z2 - 10
  expect_equal(
    marginaleffects::predictions(fit, newdata = changed)$estimate,
    expected, tolerance = 1e-8
  )

  slopes <- marginaleffects::avg_slopes(
    fit, variables = "x", newdata = dat[1:20, ]
  )
  expect_true(all(is.finite(slopes$estimate)))
  expect_true(all(is.finite(slopes$std.error)))
  analytic_slopes <- vapply(list(first, second), function(model) {
    X <- model.matrix(~x + endogenous, dat[1:20, ])
    beta <- coef(model)[seq_len(ncol(X))]
    mean(dnorm(drop(X %*% beta))) * unname(beta[2L])
  }, numeric(1))
  expect_equal(slopes$estimate, analytic_slopes, tolerance = 1e-7)

  # Full-system delta method: includes covariance between the two fits.
  comparisons <- marginaleffects::avg_comparisons(fit, variables = list(x = 1),
    newdata = dat[1:20, ])
  expect_identical(as.character(comparisons$group), c("Base", "Adjusted"))
  delta <- marginaleffects::hypotheses(comparisons,
    hypothesis = matrix(c(1, -1), ncol = 1))
  low <- high <- dat[1:20, ]
  high$x <- high$x + 1
  X0 <- model.matrix(~x + endogenous, low)
  X1 <- model.matrix(~x + endogenous, high)
  gradient <- numeric(length(coef(fit)))
  expected_delta <- 0
  for (j in seq_along(c("Base", "Adjusted"))) {
    index <- which(startsWith(names(coef(fit)),
      paste0(c("Base", "Adjusted")[j], "::")))[seq_len(ncol(X0))]
    beta <- coef(fit)[index]
    eta0 <- drop(X0 %*% beta)
    eta1 <- drop(X1 %*% beta)
    sign <- c(1, -1)[j]
    expected_delta <- expected_delta + sign * mean(pnorm(eta1) - pnorm(eta0))
    gradient[index] <- sign * colMeans(dnorm(eta1) * X1 - dnorm(eta0) * X0)
  }
  expected_delta_se <- sqrt(drop(gradient %*% vcov(fit) %*% gradient))
  expect_equal(as.numeric(delta$estimate), expected_delta, tolerance = 1e-8)
  expect_equal(as.numeric(delta$std.error), expected_delta_se, tolerance = 1e-6)
})

test_that("Rchoice IV probit exposes the structural index on the link scale", {
  skip_if_not_installed("Rchoice")
  skip_if_not_installed("mvtnorm")
  dat <- ivprobit_test_data(n = 160L)
  model <- Rchoice::ivpml(
    y ~ x + endogenous | x + z1,
    data = dat, messages = FALSE, printLevel = 0
  )
  fit <- suest(model, model, model_names = c("First", "Second"))
  prediction <- suest:::get_predict.suest_model(
    fit, dat[1:4, ], type = "link"
  )

  expect_equal(
    prediction$estimate,
    rep(predict(model, newdata = dat[1:4, ], type = "xb"), 2L),
    tolerance = 1e-10
  )
})

test_that("overidentified IV-probit bread matches the independent joint likelihood", {
  skip_if_not_installed("Rchoice")
  skip_if_not_installed("mvtnorm")
  dat <- ivprobit_test_data(n = 600L)
  model <- Rchoice::ivpml(y ~ x + endogenous | x + z1 + z2,
    data = dat, messages = FALSE, printLevel = 0)
  component <- suest:::.suest_ivprobit_components(model)
  X <- model.matrix(~x + endogenous, dat)
  Z <- model.matrix(~x + z1 + z2, dat)
  likelihood <- function(b) {
    sigma <- exp(b[8])
    rho <- tanh(b[9])
    u <- (dat$endogenous - drop(Z %*% b[4:7]))/sigma
    index <- (drop(X %*% b[1:3]) + rho*u)/sqrt(1 - rho^2)
    dnorm(u, log = TRUE) - log(sigma) + pnorm((2*dat$y - 1)*index, log.p = TRUE)
  }
  b <- coef(model)
  expect_equal(sum(likelihood(b)), as.numeric(logLik(model)), tolerance = 1e-10)
  for (step in c(1e-4, 2e-4)) {
    H <- optimHess(b, function(x) -sum(likelihood(x)),
      control = list(ndeps = rep(step, length(b))))
    expect_equal(unname(component$bread/nrow(dat)), unname(solve(H)), tolerance = 1e-6)
  }
  scores <- vapply(seq_along(b), function(j) {
    plus <- minus <- b
    plus[j] <- plus[j] + 1e-5
    minus[j] <- minus[j] - 1e-5
    (likelihood(plus) - likelihood(minus))/2e-5
  }, numeric(nrow(dat)))
  expect_equal(unname(component$score), unname(scores), tolerance = 1e-7)
})
