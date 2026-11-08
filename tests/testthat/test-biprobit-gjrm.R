gjrm_test_data <- function(n = 400L) {
  set.seed(6043)
  x <- rnorm(n)
  z <- rnorm(n)
  w <- rnorm(n)
  e <- mvtnorm::rmvnorm(n, sigma = matrix(c(1, 0.4, 0.4, 1), 2))
  d <- data.frame(id = seq_len(n),
    y1 = as.integer(0.25 + 0.65*x - 0.3*z + e[, 1] > 0),
    y2 = as.integer(-0.1 + 0.35*x + 0.4*w + e[, 2] > 0),
    x = x, z = z, w = w, f = factor(sample(c("a", "b", "c"), n, TRUE)))
  d$z[7] <- NA
  d
}

# Independent bivariate probit log likelihood: P(Y1 = y1, Y2 = y2).
gjrm_loglik <- function(b, X1, X2, y1, y2) {
  k1 <- ncol(X1)
  eta1 <- drop(X1 %*% b[seq_len(k1)])
  eta2 <- drop(X2 %*% b[k1 + seq_len(ncol(X2))])
  rho <- tanh(b[length(b)])
  vapply(seq_along(eta1), function(i) {
    lower <- ifelse(c(y1[i], y2[i]) == 1, -Inf, c(eta1[i], eta2[i]))
    upper <- ifelse(c(y1[i], y2[i]) == 1, c(eta1[i], eta2[i]), Inf)
    log(as.numeric(mvtnorm::pmvnorm(lower = -upper, upper = -lower,
      corr = matrix(c(1, rho, rho, 1), 2), algorithm = mvtnorm::Miwa(steps = 256))))
  }, numeric(1))
}

test_that("GJRM bivariate probit scores and information are exact", {
  skip_if_not_installed("GJRM")
  skip_if_not_installed("numDeriv")
  d <- gjrm_test_data()
  model <- GJRM::gjrm(list(y1 ~ x + z + f, y2 ~ x + w), data = d,
    margins = c("probit", "probit"), model = "B")
  wrapper <- suest:::.suest_wrap_gjrm(model)
  expect_identical(suest:::.suest_model_adapter(wrapper)$type, "biprobit_gjrm")
  component <- suest:::.suest_gjrm_components(wrapper)
  theta <- component$parameters
  expect_identical(names(theta), c("y1:(Intercept)", "y1:x", "y1:z", "y1:fb", "y1:fc",
    "y2:(Intercept)", "y2:x", "y2:w", "athrho"))
  expect_equal(unname(theta), unname(model$coefficients), tolerance = 1e-12)
  expect_equal(tanh(theta[["athrho"]]), unname(model$theta), tolerance = 1e-10)
  expect_identical(rownames(component$score), rownames(d)[-7])

  # The latent errors enter with opposite sign in the lower-tail form, so
  # P(Y1 = 1, Y2 = 1) = P(-e1 < eta1, -e2 < eta2) with the same correlation.
  loglik <- function(b) gjrm_loglik(b, model$X1, model$X2, model$y1, model$y2)
  # GJRM's bivariate normal probabilities are accurate to about 1e-8 relative.
  expect_equal(sum(loglik(theta)), as.numeric(logLik(model)), tolerance = 1e-7)
  for (i in c(1L, 50L, 300L)) {
    one <- function(b) gjrm_loglik(b, model$X1[i, , drop = FALSE], model$X2[i, , drop = FALSE],
      model$y1[i], model$y2[i])
    expect_equal(unname(component$score[i, ]), numDeriv::grad(one, theta), tolerance = 1e-6)
  }
  # GJRM's covariance is the inverse of its analytic Hessian.
  expect_equal(unname(component$bread/nrow(model$X1)), unname(model$Vb), tolerance = 1e-5)
})

test_that("GJRM and mvProbit give the same bivariate probit system", {
  skip_if_not_installed("GJRM")
  skip_if_not_installed("mvProbit")
  d <- gjrm_test_data()[8:157, ]
  rownames(d) <- NULL
  gjrm <- GJRM::gjrm(list(y1 ~ x + z, y2 ~ x + z), data = d,
    margins = c("probit", "probit"), model = "B")
  mvp <- mvProbit::mvProbit(cbind(y1, y2) ~ x + z, data = d,
    algorithm = mvtnorm::TVPACK(), method = "BFGS", finalHessian = TRUE,
    iterlim = 300, reltol = 1e-12, printLevel = 0)
  probit <- glm(y1 ~ x, binomial("probit"), d)
  a <- suest(gjrm, probit, model_names = c("BP", "P"))
  b <- suest(mvp, probit, model_names = c("BP", "P"))
  expect_equal(unname(coef(a)), unname(coef(b)), tolerance = 1e-4)
  expect_equal(unname(sqrt(diag(vcov(a)))), unname(sqrt(diag(vcov(b)))), tolerance = 1e-3)
  expect_equal(a$nobs_union, nrow(d))
})

test_that("GJRM bivariate probits combine, predict joint success, and find their data", {
  skip_if_not_installed("GJRM")
  skip_if_not_installed("marginaleffects")
  d <- gjrm_test_data()
  model <- GJRM::gjrm(list(y1 ~ x + z, y2 ~ x + w), data = d,
    margins = c("probit", "probit"), model = "B")
  probit <- glm(y2 ~ x + w, binomial("probit"), d)
  fit <- suest(model, probit, model_names = c("BP", "P"))
  expect_identical(unname(fit$model_types), c("biprobit_gjrm", "probit"))
  expect_equal(fit$nobs_union, nrow(d))
  expect_equal(fit$nobs_overlap, nrow(d) - 1L)

  wrapper <- suest:::.suest_wrap_gjrm(model)
  expect_identical(wrapper$call$data, quote(d))
  predicted <- suest:::.suest_predict_values(wrapper, d[-7, ][1:5, ], "response",
    "GJRM::gjrm", "biprobit_gjrm")
  expect_equal(predicted, unname(model$p11[1:5]), tolerance = 1e-6)

  averages <- marginaleffects::avg_predictions(fit, newdata = d[-7, ])
  expect_true(all(is.finite(averages$std.error)))
  changed <- coef(fit)
  changed["BP::athrho"] <- changed["BP::athrho"] + 0.2
  before <- marginaleffects::predictions(fit, newdata = d[1:3, ], vcov = FALSE)
  after <- marginaleffects::predictions(marginaleffects::set_coef(fit, changed),
    newdata = d[1:3, ], vcov = FALSE)
  rows <- before$group == "BP"
  expect_true(all(after$estimate[rows] > before$estimate[rows]))

  # Fit inside a function: the data are found in the function's environment.
  inner <- function(data) GJRM::gjrm(list(y1 ~ x + z, y2 ~ x + w), data = data,
    margins = c("probit", "probit"), model = "B")
  expect_identical(suest:::.suest_wrap_gjrm(inner(d))$call$data, quote(data))
  fit <- suest(model, probit, observation_id = "id", model_names = c("BP", "P"))
  expect_equal(fit$nobs_union, nrow(d))
})

test_that("unsupported GJRM models are refused with a reason", {
  skip_if_not_installed("GJRM")
  d <- gjrm_test_data()[-7, ]
  probit <- glm(y1 ~ x, binomial("probit"), d)
  logit <- GJRM::gjrm(list(y1 ~ x, y2 ~ x), data = d, margins = c("logit", "logit"),
    model = "B")
  expect_error(suest(logit, probit), "bivariate probits")
  clayton <- GJRM::gjrm(list(y1 ~ x, y2 ~ x), data = d, margins = c("probit", "probit"),
    model = "B", copula = "C0")
  expect_error(suest(clayton, probit), "Gaussian copula")
  smooth <- suppressWarnings(GJRM::gjrm(list(y1 ~ s(x), y2 ~ x), data = d,
    margins = c("probit", "probit"), model = "B"))
  expect_error(suest(smooth, probit), "smooth terms")
  correlation <- GJRM::gjrm(list(y1 ~ x, y2 ~ x, ~ z), data = d,
    margins = c("probit", "probit"), model = "B")
  expect_error(suest(correlation, probit), "covariates in the correlation")
})
