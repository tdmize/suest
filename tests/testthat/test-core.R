test_that("same, partial, and disjoint samples are aligned", {
  set.seed(376)
  dat <- data.frame(
    id = seq_len(600),
    year = rep(c(1986, 2016), each = 300),
    x = rnorm(600),
    z = rnorm(600)
  )
  dat$y <- rbinom(
    600,
    1,
    plogis(-0.2 + 0.4 * dat$x - 0.2 * dat$z)
  )

  same1 <- glm(y ~ x, family = binomial(), data = dat)
  same2 <- glm(y ~ x + z, family = binomial(), data = dat)
  partial1 <- glm(
    y ~ x,
    family = binomial(),
    data = dat,
    subset = id <= 450
  )
  partial2 <- glm(
    y ~ x + z,
    family = binomial(),
    data = dat,
    subset = id >= 151
  )
  disjoint1 <- glm(
    y ~ x + z,
    family = binomial(),
    data = dat,
    subset = year == 1986
  )
  disjoint2 <- glm(
    y ~ x + z,
    family = binomial(),
    data = dat,
    subset = year == 2016
  )

  same <- suest(same1, same2)
  partial <- suest(partial1, partial2)
  disjoint <- suest(disjoint1, disjoint2)

  expect_equal(same$nobs_overlap, 600)
  expect_gt(max(abs(offdiag_vcov(same))), 0)

  expect_equal(partial$nobs_overlap, 300)
  expect_equal(partial$nobs_union, 600)
  expect_gt(max(abs(offdiag_vcov(partial))), 0)

  expect_equal(disjoint$nobs_overlap, 0)
  expect_true(all(offdiag_vcov(disjoint) == 0))
})

test_that("model-specific covariance blocks are preserved", {
  dat <- mtcars
  dat$am <- factor(dat$am)

  model1 <- glm(
    am ~ wt + hp,
    family = binomial(),
    data = dat,
    subset = cyl == 4
  )
  model2 <- glm(
    am ~ wt + hp,
    family = binomial(),
    data = dat,
    subset = cyl != 4
  )

  fit <- suest(model1, model2)
  V <- vcov(fit)

  expect_equal(
    unname(V[fit$index[[1]], fit$index[[1]]]),
    unname(robust_vcov_direct(model1)),
    tolerance = 1e-10
  )
  expect_equal(
    unname(V[fit$index[[2]], fit$index[[2]]]),
    unname(robust_vcov_direct(model2)),
    tolerance = 1e-10
  )
})

test_that("unsupported model combinations are rejected", {
  dat <- mtcars
  dat$am <- factor(dat$am)

  logit <- glm(am ~ wt, family = binomial(), data = dat)
  poisson <- glm(cyl ~ wt, family = poisson(), data = dat)

  expect_error(
    suest(logit, poisson),
    "not supported"
  )
})
