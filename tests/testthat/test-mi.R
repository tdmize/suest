mi_test_fits <- function(m = 5L) {
  set.seed(1907)
  n <- 180L
  x <- rnorm(n)
  z <- rnorm(n)
  y <- 0.4 + 0.7*x - 0.3*z + rnorm(n)
  lapply(seq_len(m), function(j) {
    dat <- data.frame(y, x = x + rnorm(n, sd = 0.03*j), z)
    suest(lm(y ~ x, dat), lm(y ~ x + z, dat),
      model_names = c("Base", "Adjusted"))
  })
}

test_that("suest_mi reproduces the complete Rubin covariance calculation", {
  fits <- mi_test_fits()
  pooled <- suest_mi(fits)
  estimates <- vapply(fits, coef, numeric(length(coef(fits[[1L]]))))
  rownames(estimates) <- names(coef(fits[[1L]]))
  within <- Reduce(`+`, lapply(fits, vcov))/length(fits)
  between <- stats::cov(t(estimates))
  total <- within + (1 + 1/length(fits))*between

  expect_s3_class(pooled, "suest_mi")
  expect_equal(coef(pooled), rowMeans(estimates), tolerance = 1e-14)
  expect_equal(pooled$within_vcov, within, tolerance = 1e-14)
  expect_equal(pooled$between_vcov, between, tolerance = 1e-14)
  expect_equal(vcov(pooled), total, tolerance = 1e-14)
  expect_gt(abs(pooled$within_vcov[1L, 4L]), 0)
  expect_gt(abs(pooled$between_vcov[1L, 4L]), 0)
  expect_identical(nobs(pooled), fits[[1L]]$nobs_union)
  expect_identical(pooled$m, 5L)

  added <- (1 + 1/pooled$m)*diag(between)
  expected_df <- (pooled$m - 1)*(1 + diag(within)/added)^2
  expect_equal(pooled$degrees_freedom, expected_df, tolerance = 1e-12)
  expect_equal(pooled$fraction_missing_information,
    added/diag(total), tolerance = 1e-14)
})

test_that("suest_mi summary uses Rubin degrees of freedom", {
  pooled <- suest_mi(mi_test_fits())
  out <- summary(pooled, conf.level = 0.90)
  expected_se <- sqrt(diag(vcov(pooled)))
  expected_statistic <- coef(pooled)/expected_se

  expect_identical(out$term, names(coef(pooled)))
  expect_equal(out$std.error, unname(expected_se), tolerance = 1e-14)
  expect_equal(out$statistic, unname(expected_statistic), tolerance = 1e-14)
  expect_equal(out$p.value,
    unname(2*pt(-abs(expected_statistic), pooled$degrees_freedom)),
    tolerance = 1e-14)
  expect_error(summary(pooled, conf.level = 1), "strictly between")
  expect_output(print(pooled), "Imputations: 5", fixed = TRUE)
})

test_that("suest_mi validates complete system compatibility", {
  fits <- mi_test_fits(3L)
  expect_error(suest_mi(fits[1L]), "at least two")
  expect_error(suest_mi(list(fits[[1L]], lm(mpg ~ wt, mtcars))),
    "Every element")

  changed <- fits
  changed[[2L]]$model_names[1L] <- "Different"
  expect_error(suest_mi(changed), "model_names differs")
  changed <- fits
  changed[[2L]]$coefficients <- changed[[2L]]$coefficients[-1L]
  expect_error(suest_mi(changed), "coefficient names or ordering")
  changed <- fits
  changed[[2L]]$vcov[1L, 1L] <- NA_real_
  expect_error(suest_mi(changed), "nonfinite covariance")
})

test_that("suest_mi accepts current mira wrappers and getfit-style lists", {
  fits <- mi_test_fits(3L)
  expected <- suest_mi(fits)
  wrapper <- structure(list(
    call = quote(with(imp, suest(lm(y ~ x), lm(y ~ x + z)))),
    call1 = quote(mice(data, m = 3)),
    nmis = c(y = 0L, x = 12L, z = 0L),
    analyses = fits
  ), class = "mira")
  from_wrapper <- suest_mi(wrapper)

  expect_equal(coef(from_wrapper), coef(expected), tolerance = 1e-14)
  expect_equal(vcov(from_wrapper), vcov(expected), tolerance = 1e-14)
  expect_identical(from_wrapper$input_source, "mice::mira")
  expect_identical(from_wrapper$mice_call, wrapper$call)
  expect_identical(from_wrapper$mice_call1, wrapper$call1)
  expect_output(print(from_wrapper), "Input: mice::mira", fixed = TRUE)

  getfit_style <- fits
  class(getfit_style) <- c("mira", "list")
  from_getfit <- suest_mi(getfit_style)
  expect_equal(coef(from_getfit), coef(expected), tolerance = 1e-14)
  expect_equal(vcov(from_getfit), vcov(expected), tolerance = 1e-14)

  malformed <- wrapper
  malformed$analyses <- fits[1L]
  expect_error(suest_mi(malformed), "at least two")
})

test_that("suest_mi accepts a real mice with workflow", {
  skip_if_not_installed("mice")
  # Avoid predictive mean matching here so the integration test also works
  # with older lightweight mice releases used on legacy R installations.
  imp <- mice::mice(mice::nhanes, m = 2L, maxit = 1L,
    method = c("", "norm", "norm", "norm"),
    printFlag = FALSE, seed = 1908)
  analyses <- with(imp, suest(
    lm(bmi ~ age + hyp),
    lm(bmi ~ age + hyp + chl),
    model_names = c("Base", "Adjusted")
  ))
  pooled <- suest_mi(analyses)
  pooled_list <- suest_mi(analyses$analyses)

  expect_s3_class(analyses, "mira")
  expect_s3_class(pooled, "suest_mi")
  expect_identical(pooled$m, 2L)
  expect_identical(pooled$input_source, "mice::mira")
  expect_equal(coef(pooled), coef(pooled_list), tolerance = 1e-14)
  expect_equal(vcov(pooled), vcov(pooled_list), tolerance = 1e-14)
  expect_equal(vcov(suest_mi(mice::getfit(analyses))), vcov(pooled),
    tolerance = 1e-14)
})

test_that("suest_mi accepts a real mitools imputationList workflow", {
  skip_if_not_installed("mitools")
  set.seed(1909)
  n <- 160L
  base <- data.frame(x = rnorm(n), z = rnorm(n))
  base$y <- 0.5 + 0.8*base$x - 0.25*base$z + rnorm(n)
  datasets <- lapply(seq_len(4L), function(j) transform(base,
    x = x + rnorm(n, sd = 0.02*j)))
  imputations <- mitools::imputationList(datasets)
  analyses <- with(imputations, suest(
    lm(y ~ x), lm(y ~ x + z), model_names = c("Base", "Adjusted")
  ))
  pooled <- suest_mi(analyses)
  reference <- mitools::MIcombine(analyses)

  expect_identical(class(analyses), "list")
  expect_false(is.null(attr(analyses, "call")))
  expect_identical(pooled$input_source, "mitools::imputationList")
  expect_identical(pooled$mitools_call, attr(analyses, "call"))
  expect_equal(coef(pooled), coef(reference), tolerance = 1e-14)
  expect_equal(vcov(pooled), vcov(reference), tolerance = 1e-14)
  expect_gt(abs(pooled$within_vcov[1L, 4L]), 0)
  expect_gt(abs(pooled$between_vcov[1L, 4L]), 0)
  expect_output(print(pooled), "Input: mitools::imputationList", fixed = TRUE)

  legacy <- unclass(analyses)
  class(legacy) <- "imputationResultList"
  legacy$call <- attr(analyses, "call")
  from_legacy <- suest_mi(legacy)
  expect_equal(coef(from_legacy), coef(pooled), tolerance = 1e-14)
  expect_equal(vcov(from_legacy), vcov(pooled), tolerance = 1e-14)
})

test_that("coefficient hypotheses work but MI predictions fail explicitly", {
  pooled <- suest_mi(mi_test_fits())
  hypothesis <- suppressWarnings(
    marginaleffects::hypotheses(pooled, hypothesis = "b1 = b2")
  )
  gradient <- c(1, -1, rep(0, length(coef(pooled)) - 2L))
  expect_equal(as.numeric(hypothesis$estimate),
    unname(coef(pooled)[1L] - coef(pooled)[2L]), tolerance = 1e-12)
  expect_equal(hypothesis$std.error,
    sqrt(drop(gradient %*% vcov(pooled) %*% gradient)), tolerance = 1e-8)
  expect_error(marginaleffects::predictions(pooled, newdata = mtcars[1L, ]),
    "not yet implemented")
})
