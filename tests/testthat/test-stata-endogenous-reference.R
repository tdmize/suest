# Returned Stata 19.5 benchmark, including the cross-model covariance blocks.
test_that("IV probit matches the full returned Stata coefficient and covariance system", {
  skip_if_not_installed("Rchoice")
  fixture <- readRDS(test_path("fixtures", "endogenous-stata.rds"))
  dat <- fixture$data
  first <- Rchoice::ivpml(binary_iv ~ x + endogenous | x + z1, data = dat,
    messages = FALSE, printLevel = 0, reltol = 1e-12)
  second <- Rchoice::ivpml(binary_iv ~ x + endogenous | x + z1 + z2, data = dat,
    messages = FALSE, printLevel = 0, reltol = 1e-12)
  fit <- suest(first, second)
  order <- c(3, 2, 1, 5, 6, 4, 8, 7, 11, 10, 9, 13, 14, 15, 12, 17, 16)
  ref <- fixture$reference$ivprobit
  expect_lt(max(abs(coef(fit)[order] - ref$b)), 1e-6)
  expect_lt(max(abs(vcov(fit)[order, order] - ref$V)), 2e-8)
  expect_lt(max(abs(diag(vcov(fit)[order, order])/diag(ref$V) - 1)), 2e-6)

  # A corrupted upstream Hessian must not enter the recomputed sandwich.
  second$hessian <- second$hessian * 0.5
  expect_equal(vcov(suest(first, second)), vcov(fit), tolerance = 1e-12)
})

