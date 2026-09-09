.biprobit_test_cache <- new.env(parent = emptyenv())

biprobit_test_fit <- function() {
  if (exists("object", envir = .biprobit_test_cache, inherits = FALSE))
    return(get("object", envir = .biprobit_test_cache, inherits = FALSE))
  set.seed(6041)
  n <- 80L
  x <- rnorm(n)
  z <- rnorm(n)
  errors <- mvtnorm::rmvnorm(
    n,
    sigma = matrix(c(1, 0.4, 0.4, 1), nrow = 2L)
  )
  y1 <- as.integer(0.25 + 0.65*x - 0.15*z + errors[, 1L] > 0)
  y2 <- as.integer(-0.10 + 0.35*x + 0.20*z + errors[, 2L] > 0)
  dat <- data.frame(y1, y2, x, z)
  model <- mvProbit::mvProbit(
    cbind(y1, y2) ~ x + z,
    data = dat,
    algorithm = mvtnorm::TVPACK(),
    method = "BFGS",
    finalHessian = TRUE,
    iterlim = 200,
    reltol = 1e-10,
    printLevel = 0
  )
  object <- list(model = model, data = dat)
  assign("object", object, envir = .biprobit_test_cache)
  object
}

test_that("bivariate probit transforms rho to Stata's athrho parameterization", {
  skip_if_not_installed("mvProbit")
  skip_if_not_installed("mvtnorm")
  object <- biprobit_test_fit()
  model <- object$model
  component <- suest:::.suest_biprobit_components(model)
  rho <- unname(coef(model)[length(coef(model))])

  expect_equal(
    unname(component$parameters["athrho"]),
    atanh(rho),
    tolerance = 1e-12
  )
  expect_equal(
    component$score[, "athrho"],
    model$gradientObs[, ncol(model$gradientObs)]*(1 - rho^2),
    tolerance = 1e-12
  )
  expect_true(all(is.finite(component$bread)))
})

test_that("bivariate probit systems predict joint success with marginaleffects", {
  skip_if_not_installed("mvProbit")
  skip_if_not_installed("mvtnorm")
  object <- biprobit_test_fit()
  model <- object$model
  dat <- object$data
  fit <- suest(model, model, model_names = c("First", "Second"))

  expect_equal(fit$nobs_overlap, nrow(dat))
  expect_true(all(c("First::athrho", "Second::athrho") %in%
                    names(coef(fit))))
  expect_true(all(is.finite(vcov(fit))))
  p <- length(coef(model))
  expect_equal(vcov(fit)[seq_len(p), p + seq_len(p)],
               vcov(fit)[seq_len(p), seq_len(p)], ignore_attr = TRUE,
               tolerance = 1e-12)
  expect_identical(suest:::.suest_data_source(model),
                   suest:::.suest_data_source(model))

  prediction <- marginaleffects::predictions(fit, newdata = dat[1L, ])
  k <- model$nReg
  X <- model.matrix(~x + z, dat[1L, , drop = FALSE])
  eta <- c(
    as.numeric(X %*% coef(model)[seq_len(k)]),
    as.numeric(X %*% coef(model)[k + seq_len(k)])
  )
  rho <- unname(coef(model)[length(coef(model))])
  expected <- as.numeric(mvtnorm::pmvnorm(
    upper = eta,
    sigma = matrix(c(1, rho, rho, 1), nrow = 2L),
    algorithm = mvtnorm::TVPACK()
  ))
  expect_equal(prediction$estimate, rep(expected, 2L), tolerance = 1e-8)

  slopes <- marginaleffects::avg_slopes(
    fit, variables = "x", newdata = dat[seq_len(10L), ]
  )
  expect_true(all(is.finite(slopes$estimate)))
  expect_true(all(is.finite(slopes$std.error)))
})

test_that("bivariate probit requires an observed final Hessian", {
  skip_if_not_installed("mvProbit")
  skip_if_not_installed("mvtnorm")
  model <- biprobit_test_fit()$model
  model$call$finalHessian <- "BHHH"
  expect_error(
    suest(model, model),
    "finalHessian = TRUE",
    fixed = TRUE
  )
})
