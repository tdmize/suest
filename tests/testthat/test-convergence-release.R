# Removing an admission guard must fail the corresponding low-iteration test.
release_convergence_data <- function() {
  set.seed(129)
  d <- data.frame(x = rnorm(300), z = rnorm(300))
  d$y <- ordered(sample(1:3, 300, TRUE))
  d$n <- rnbinom(300, mu = exp(.4 + .4*d$x), size = .7)
  d
}

test_that("nonconverged multinom fits are rejected", {
  skip_if_not_installed("nnet")
  d <- release_convergence_data()
  bad <- nnet::multinom(y ~ x + z, data = d, trace = FALSE, maxit = 1)
  expect_equal(bad$convergence, 1)
  expect_error(suest(bad, bad), "did not converge")
  good <- nnet::multinom(y ~ x + z, data = d, trace = FALSE, maxit = 200)
  expect_s3_class(suest(good, good), "suest_model")
})

test_that("nonconverged polr fits are rejected", {
  d <- release_convergence_data()
  bad <- MASS::polr(y ~ x + z, data = d, Hess = TRUE, control = list(maxit = 1))
  expect_equal(bad$convergence, 1)
  expect_error(suest(bad, bad), "did not converge")
  good <- MASS::polr(y ~ x + z, data = d, Hess = TRUE)
  expect_s3_class(suest(good, good), "suest_model")
})

test_that("failed clm convergence diagnostics are rejected", {
  skip_if_not_installed("ordinal")
  d <- release_convergence_data()
  bad <- suppressWarnings(ordinal::clm(y ~ x + z, data = d,
    control = ordinal::clm.control(maxIter = 1)))
  expect_lt(bad$convergence$code, 0)
  expect_error(suest(bad, bad), "did not converge")
  good <- ordinal::clm(y ~ x + z, data = d)
  expect_s3_class(suest(good, good), "suest_model")
})

test_that("negative binomial alternation failures are rejected", {
  d <- release_convergence_data()
  bad <- suppressWarnings(MASS::glm.nb(n ~ x + z, data = d,
    control = glm.control(maxit = 2)))
  expect_match(bad$th.warn, "limit reached")
  expect_error(suest(bad, bad), "did not converge")
  good <- MASS::glm.nb(n ~ x + z, data = d)
  expect_s3_class(suest(good, good), "suest_model")
  good$converged <- FALSE
  expect_error(suest(good, good), "did not converge")
})
