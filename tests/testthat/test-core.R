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

test_that("observation IDs align models fitted from different data objects", {
  set.seed(317)
  dat <- data.frame(
    id = seq_len(600),
    x = rnorm(600),
    z = rnorm(600)
  )
  dat$y <- rbinom(600, 1, plogis(-0.2 + 0.5 * dat$x - 0.3 * dat$z))
  data1 <- dat[dat$id <= 450, ]
  data2 <- dat[dat$id >= 151, ]

  separate1 <- glm(y ~ x, family = binomial(), data = data1)
  separate2 <- glm(y ~ x + z, family = binomial(), data = data2)
  common1 <- glm(
    y ~ x,
    family = binomial(),
    data = dat,
    subset = id <= 450
  )
  common2 <- glm(
    y ~ x + z,
    family = binomial(),
    data = dat,
    subset = id >= 151
  )

  default <- suest(separate1, separate2)
  reference <- suest(common1, common2)
  aligned <- suest(
    separate1,
    separate2,
    model_names = c("A", "B"),
    observation_id = "id"
  )
  supplied <- suest(
    separate1,
    separate2,
    model_names = c("A", "B"),
    observation_id = list(A = data1$id, B = data2$id)
  )

  expect_equal(default$nobs_overlap, 0)
  expect_equal(aligned$nobs_overlap, 300)
  expect_equal(aligned$nobs_union, 600)
  expect_equal(unname(aligned$vcov), unname(reference$vcov), tolerance = 1e-12)
  expect_equal(supplied$vcov, aligned$vcov, tolerance = 1e-12)
})

test_that("composite observation IDs are supported", {
  set.seed(727)
  dat <- data.frame(
    person = rep(seq_len(150), each = 3),
    wave = rep(seq_len(3), 150),
    x = rnorm(450),
    z = rnorm(450)
  )
  dat$y <- 1 + 0.4 * dat$x - 0.2 * dat$z + rnorm(450)
  data1 <- dat[dat$wave <= 2, ]
  data2 <- dat[dat$wave >= 2, ]

  model1 <- lm(y ~ x, data = data1)
  model2 <- lm(y ~ x + z, data = data2)
  fit <- suest(
    model1,
    model2,
    observation_id = c("person", "wave")
  )

  expect_equal(fit$nobs_overlap, 150)
  expect_equal(fit$nobs_union, 450)
})

test_that("linear ancillary parameters match the Stata suest2 reference", {
  set.seed(3371)
  n <- 1200L
  x <- rnorm(n)
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  disturbance <- rnorm(n)
  endogenous <- 0.7*z1 + 0.25*z2 + 0.4*disturbance + rnorm(n)
  y <- 0.45 + 0.65*x + 1.15*endogenous + disturbance
  hetero_probability <- pnorm((0.2 + 0.7*x - 0.25*z2)/exp(0.35*z1))
  benchmark <- data.frame(
    id = seq_len(n),
    y = y,
    x = x,
    endogenous = endogenous,
    z1 = z1,
    z2 = z2,
    y_het = rbinom(n, 1, hetero_probability)
  )
  base <- lm(y ~ x, data = benchmark)
  adjusted <- lm(y ~ x + z1, data = benchmark)
  fit <- suest(base, adjusted, model_names = c("Base", "Adjusted"))

  expected_b <- c(
    `Base::(Intercept)` = 0.484446605779495,
    `Base::x` = 0.658672818262305,
    `Base::lnvar` = 1.4496969842131,
    `Adjusted::(Intercept)` = 0.462877426932272,
    `Adjusted::x` = 0.676895799208494,
    `Adjusted::z1` = 0.753540029894244,
    `Adjusted::lnvar` = 1.31654110370512
  )
  expect_equal(coef(fit), expected_b, tolerance = 5e-9)

  # Lower triangle copied from the Stata log in Stata's display order:
  # base x, intercept, lnvar; adjusted x, z1, intercept, lnvar.
  lower <- c(
    0.00318639311789977,
    0.000161005378427106, 0.00354748979959301,
    -0.0000113592979373545, -0.000190647476596527,
    0.0017472882690647,
    0.00279860956106366, 0.000126512170871698,
    0.0000605202239118408, 0.00284454928696272,
    -0.000201652333688137, -0.0000336649682425002,
    0.00118854869842638, -0.000104353607675721,
    0.00350441767970183,
    0.000133098352266485, 0.00310367750371945,
    -0.000214930468180607, 0.0000945706618101833,
    -0.000129151337297262, 0.00310746088770924,
    0.0000745775261584513, -0.000200286599734616,
    0.00150995500189096, 0.00011096235472822,
    0.0000203790323224869, -0.000195395710492785,
    0.00171002106463944
  )
  expected_v <- matrix(0, nrow = 7L, ncol = 7L)
  cursor <- 1L
  for (i in seq_len(7L)) {
    for (j in seq_len(i)) {
      expected_v[i, j] <- lower[cursor]
      expected_v[j, i] <- lower[cursor]
      cursor <- cursor + 1L
    }
  }
  r_order <- c(2L, 1L, 3L, 6L, 4L, 5L, 7L)
  expected_v <- expected_v[r_order, r_order]
  expect_lt(max(abs(unname(vcov(fit)) - expected_v)), 3e-11)

  effects <- marginaleffects::avg_slopes(
    fit,
    variables = "x",
    newdata = benchmark
  )
  expect_equal(nrow(effects), 2L)
  expect_true(all(is.finite(effects$std.error)))
})

test_that("invalid observation IDs are rejected", {
  dat <- mtcars
  model1 <- lm(mpg ~ wt, data = dat)
  model2 <- lm(mpg ~ wt + hp, data = dat)

  expect_error(
    suest(model1, model2, observation_id = list(seq_len(nrow(dat)))),
    "one ID vector"
  )
  expect_error(
    suest(
      model1,
      model2,
      observation_id = list(rep(1, nrow(dat)), seq_len(nrow(dat)))
    ),
    "unique within model"
  )
  expect_error(
    suest(
      model1,
      model2,
      observation_id = list(c(NA, seq_len(nrow(dat) - 1L)), seq_len(nrow(dat)))
    ),
    "cannot be missing"
  )
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

test_that("scalar and categorical response structures can be mixed", {
  dat <- mtcars
  dat$am <- factor(dat$am)
  dat$cyl_ord <- ordered(dat$cyl)

  logit <- glm(am ~ wt, family = binomial(), data = dat)
  ordered <- MASS::polr(
    cyl_ord ~ wt,
    data = dat,
    method = "logistic",
    Hess = TRUE
  )

  fit <- suest(logit, ordered, model_names = c("Logit", "Ordered"))
  effects <- marginaleffects::avg_comparisons(
    fit,
    variables = "wt",
    newdata = dat
  )

  expect_equal(nrow(effects), 4L)
  expect_true(all(is.finite(effects$estimate)))
  expect_true(all(is.finite(effects$std.error)))
})
