test_that("ordinary GLMs must converge before combining", {
  set.seed(20260928)
  dat <- data.frame(x = runif(400, -1, 1))
  dat$y <- rbinom(nrow(dat), 1, plogis(-.2 + .8 * dat$x))
  good <- glm(y ~ x, data = dat, family = binomial())
  expect_true(good$converged)
  expect_s3_class(suest(good, good), "suest_model")

  expect_warning(
    bad <- glm(y ~ x, data = dat, family = binomial(),
      control = glm.control(maxit = 1)),
    "did not converge"
  )
  expect_false(bad$converged)
  expect_error(suest(good, bad), "converg")
})

test_that("glm2 fits must converge before combining", {
  skip_if_not_installed("glm2")
  set.seed(20260928)
  dat <- data.frame(x = runif(400, -1, 1))
  dat$y <- rbinom(nrow(dat), 1, plogis(-.2 + .8 * dat$x))
  good <- glm2::glm2(y ~ x, data = dat, family = binomial())
  expect_true(good$converged)
  expect_s3_class(suest(good, good), "suest_model")
  expect_warning(
    bad <- glm2::glm2(y ~ x, data = dat, family = binomial(),
      control = glm.control(maxit = 1)),
    "did not converge"
  )
  expect_false(bad$converged)
  expect_error(suest(good, bad), "converg")
})

test_that("multinomial weight decay cannot use ordinary ML scores", {
  skip_if_not_installed("nnet")
  set.seed(20260928)
  dat <- data.frame(x = runif(400, -1, 1))
  dat$y <- factor(sample(c("A", "B", "C"), nrow(dat), replace = TRUE))
  ordinary <- nnet::multinom(y ~ x, data = dat, decay = 0,
    trace = FALSE, Hess = TRUE, maxit = 500)
  penalized <- nnet::multinom(y ~ x, data = dat, decay = 2,
    trace = FALSE, Hess = TRUE, maxit = 500)
  expect_equal(ordinary$convergence, 0)
  expect_equal(penalized$convergence, 0)
  expect_s3_class(suest(ordinary, ordinary), "suest_model")
  expect_error(suest(ordinary, penalized), "[Pp]enal|decay")
})
