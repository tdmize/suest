ivreg_test_data <- function(n = 1200L) {
  set.seed(3371)
  x <- rnorm(n)
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  disturbance <- rnorm(n)
  endogenous <- 0.7*z1 + 0.25*z2 + 0.4*disturbance + rnorm(n)
  y <- 0.45 + 0.65*x + 1.15*endogenous + disturbance

  data.frame(
    id = seq_len(n),
    panel = factor(rep(seq_len(n/4L), each = 4L)),
    y = y,
    x = x,
    endogenous = endogenous,
    z1 = z1,
    z2 = z2
  )
}

test_that("fixest 2SLS scores and bread equal the IV equations", {
  skip_if_not_installed("fixest")
  dat <- ivreg_test_data()
  model <- fixest::feols(y ~ x | endogenous ~ z1 + z2, data = dat)

  X <- cbind(1, dat$endogenous, dat$x)
  Z <- cbind(1, dat$z1, dat$z2, dat$x)
  Xhat <- qr.fitted(qr(Z), X)
  beta <- solve(crossprod(Xhat, X), crossprod(Xhat, dat$y))
  residual <- as.numeric(dat$y - X %*% beta)
  score <- Xhat*residual
  bread <- nrow(dat)*solve(crossprod(Xhat, X))

  expect_equal(
    unname(stats::coef(model)),
    as.numeric(beta),
    tolerance = 1e-10
  )
  expect_equal(
    unname(sandwich::estfun(model)),
    unname(score),
    tolerance = 1e-10
  )
  expect_equal(
    unname(sandwich::bread(model)),
    unname(bread),
    tolerance = 1e-10
  )
})

test_that("fixest 2SLS preserves covariance and marginaleffects", {
  skip_if_not_installed("fixest")
  dat <- ivreg_test_data()
  base <- fixest::feols(y ~ x | endogenous ~ z1, data = dat)
  adjusted <- fixest::feols(
    y ~ x | endogenous ~ z1 + z2,
    data = dat
  )
  fit <- suest(base, adjusted, model_names = c("Base", "Adjusted"))

  expect_equal(unname(fit$model_types), rep("ivreg", 2L))
  expect_equal(unname(fit$model_engines), rep("fixest::feols", 2L))
  expect_equal(
    unname(vcov(fit)[fit$index[[1L]], fit$index[[1L]], drop = FALSE]),
    unname(robust_vcov_direct(base, correction = FALSE)),
    tolerance = 1e-8
  )
  expect_equal(
    unname(vcov(fit)[fit$index[[2L]], fit$index[[2L]], drop = FALSE]),
    unname(robust_vcov_direct(adjusted, correction = FALSE)),
    tolerance = 1e-8
  )

  effects <- marginaleffects::avg_comparisons(
    fit,
    variables = "x",
    newdata = dat
  )
  expect_equal(nrow(effects), 2L)
  expect_equal(
    effects$estimate,
    unname(c(stats::coef(base)["x"], stats::coef(adjusted)["x"])),
    tolerance = 1e-5
  )
  expect_true(all(is.finite(effects$std.error)))
})

test_that("fixest 2SLS matches the Stata suest2 reference", {
  skip_if_not_installed("fixest")
  dat <- ivreg_test_data()
  base <- fixest::feols(y ~ x | endogenous ~ z1, data = dat)
  adjusted <- fixest::feols(
    y ~ x | endogenous ~ z1 + z2,
    data = dat
  )
  fit <- suest(base, adjusted, model_names = c("Base", "Adjusted"))

  expected_b <- c(
    `Base::(Intercept)` = 0.448580187579429,
    `Base::fit_endogenous` = 1.06691818679677,
    `Base::x` = 0.659181399605071,
    `Adjusted::(Intercept)` = 0.447632033544247,
    `Adjusted::fit_endogenous` = 1.0951229207732,
    `Adjusted::x` = 0.659194844311193
  )
  expect_equal(coef(fit), expected_b, tolerance = 1e-8)

  V <- vcov(fit)
  first <- fit$index[[1L]]
  second <- fit$index[[2L]]
  names(first) <- fit$local_names[[1L]]
  names(second) <- fit$local_names[[2L]]
  expect_equal(
    c(
      V[first["fit_endogenous"], first["fit_endogenous"]],
      V[first["x"], first["x"]],
      V[second["fit_endogenous"], second["fit_endogenous"]],
      V[second["x"], second["x"]],
      V[first["fit_endogenous"], second["fit_endogenous"]],
      V[first["x"], second["x"]],
      V[first["(Intercept)"], second["(Intercept)"]]
    ),
    c(
      0.00205632246455327,
      0.00082722527208704,
      0.00169041773123478,
      0.000804515226974383,
      0.00171497102466982,
      0.000815361439295064,
      0.000902651213387818
    ),
    tolerance = 2e-8
  )
})

test_that("fixest 2SLS samples align through original observation IDs", {
  skip_if_not_installed("fixest")
  dat <- ivreg_test_data()
  left_data <- dat[dat$id <= 900L, ]
  right_data <- dat[dat$id > 300L, ]
  left <- fixest::feols(
    y ~ x | endogenous ~ z1 + z2,
    data = left_data
  )
  right <- fixest::feols(
    y ~ x | endogenous ~ z1 + z2,
    data = right_data
  )
  fit <- suest(left, right, observation_id = "id")

  expect_equal(fit$nobs_overlap, 600L)
  expect_true(any(abs(offdiag_vcov(fit)) > 1e-12))
  expect_true(all(is.finite(vcov(fit))))
})

test_that("unsupported fixest specifications are rejected", {
  skip_if_not_installed("fixest")
  dat <- ivreg_test_data()
  ordinary <- fixest::feols(y ~ x + endogenous, data = dat)
  fixed_effect <- fixest::feols(
    y ~ x | panel | endogenous ~ z1 + z2,
    data = dat
  )
  valid <- fixest::feols(y ~ x | endogenous ~ z1 + z2, data = dat)

  expect_error(suest(ordinary, valid), "supported")
  expect_error(suest(fixed_effect, valid), "absorbed fixed effects")
})
