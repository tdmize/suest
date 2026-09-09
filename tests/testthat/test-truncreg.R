truncreg_test_data <- function(n = 2400L) {
  set.seed(5264)
  x <- rnorm(n)
  z <- rnorm(n)
  data.frame(
    id = seq_len(n),
    x = x,
    z = z,
    y = 0.4 + 0.7*x - 0.3*z + rnorm(n)
  )
}

test_that("left-truncated regressions preserve robust covariance", {
  skip_if_not_installed("truncreg")
  full <- truncreg_test_data()
  dat <- full[full$y > 0, ]
  base <- truncreg::truncreg(
    y ~ x,
    data = dat,
    point = 0,
    direction = "left",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  adjusted <- truncreg::truncreg(
    y ~ x + z,
    data = dat,
    point = 0,
    direction = "left",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  fit <- suest(base, adjusted)

  expect_equal(unname(fit$model_types), rep("truncreg", 2L))
  expect_true(any(grepl("::sigma$", names(coef(fit)))))
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

test_that("left- and right-truncated samples align by observation ID", {
  skip_if_not_installed("truncreg")
  full <- truncreg_test_data()
  left_data <- full[full$y > 0, ]
  right_data <- full[full$y < 1, ]
  left <- truncreg::truncreg(
    y ~ x + z,
    data = left_data,
    point = 0,
    direction = "left",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  right <- truncreg::truncreg(
    y ~ x + z,
    data = right_data,
    point = 1,
    direction = "right",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  fit <- suest(left, right, observation_id = "id")

  expect_equal(
    fit$nobs_overlap,
    length(intersect(left_data$id, right_data$id))
  )
  expect_true(any(abs(offdiag_vcov(fit)) > 1e-12))
  expect_true(all(is.finite(vcov(fit))))
})
