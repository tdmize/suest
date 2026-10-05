heckman_test_data <- function() {
  set.seed(1011)
  n <- 600L
  d <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n))
  e1 <- rnorm(n)
  e2 <- 0.5*e1 + sqrt(1 - 0.25)*rnorm(n)
  d$s <- (0.3 + 0.5*d$x + 0.8*d$w + e1) > 0
  d$y <- ifelse(d$s, 1 + 0.6*d$x - 0.4*d$z + 1.2*e2, NA)
  d$yc <- 1 + 0.5*d$x + rnorm(n)
  d
}

# Independent Heckman log likelihood per observation in Stata's parameters.
heckman_hand_loglik <- function(theta, d) {
  g <- theta[1:3]; b <- theta[4:6]
  sigma <- exp(theta[7]); rho <- tanh(theta[8])
  zg <- g[1] + g[2]*d$x + g[3]*d$w
  out <- stats::pnorm(-zg, log.p = TRUE)
  sel <- d$s
  r <- (d$y[sel] - (b[1] + b[2]*d$x[sel] + b[3]*d$z[sel]))/sigma
  out[sel] <- stats::pnorm((zg[sel] + rho*r)/sqrt(1 - rho^2), log.p = TRUE) +
    stats::dnorm(r, log = TRUE) - log(sigma)
  out
}

test_that("Heckman scores and information match an independent likelihood", {
  skip_if_not_installed("sampleSelection")
  d <- heckman_test_data()
  model <- sampleSelection::selection(s ~ x + w, y ~ x + z, data = d, method = "ml")
  expect_identical(suest:::.suest_heckman_adapter(model)$type, "heckman")
  component <- suest:::.suest_heckman_components(suest:::.suest_prepare_heckman(model))
  theta <- component$parameters
  expect_identical(names(theta)[7:8], c("lnsigma", "athrho"))
  expect_equal(sum(heckman_hand_loglik(theta, d)), as.numeric(stats::logLik(model)),
    tolerance = 1e-10)
  numerical <- sapply(seq_along(theta), function(k) {
    up <- down <- theta
    up[k] <- up[k] + 1e-6
    down[k] <- down[k] - 1e-6
    (heckman_hand_loglik(up, d) - heckman_hand_loglik(down, d))/2e-6
  })
  expect_equal(unname(component$score), unname(numerical), tolerance = 1e-6)
  total <- function(t) sum(heckman_hand_loglik(t, d))
  hessian <- sapply(seq_along(theta), function(k) sapply(seq_along(theta), function(m) {
    h <- 1e-4
    e <- function(i) replace(numeric(length(theta)), i, h)
    (total(theta + e(k) + e(m)) - total(theta + e(k) - e(m)) -
      total(theta - e(k) + e(m)) + total(theta - e(k) - e(m)))/(4*h^2)
  }))
  expect_equal(unname(component$bread/nrow(component$score)), unname(solve(-hessian)),
    tolerance = 1e-4)
})

test_that("Heckman models combine with others and predict the outcome equation", {
  skip_if_not_installed("sampleSelection")
  d <- heckman_test_data()
  model <- sampleSelection::selection(s ~ x + w, y ~ x + z, data = d, method = "ml")
  fit <- suest(model, lm(yc ~ x, d), model_names = c("heckman", "regress"))
  expect_identical(fit$comparison_scale, "fitted values")
  rows <- d[1:3, ]
  p <- coef(fit)
  expected <- p[["heckman::outcome:(Intercept)"]] + p[["heckman::outcome:x"]]*rows$x +
    p[["heckman::outcome:z"]]*rows$z
  predicted <- marginaleffects::predictions(fit, newdata = rows)
  expect_equal(predicted$estimate[predicted$group == "heckman"], expected, tolerance = 1e-12)
  slope <- marginaleffects::avg_slopes(fit, variables = "x")
  expect_equal(slope$estimate[1], p[["heckman::outcome:x"]], tolerance = 1e-6)
  expect_equal(slope$std.error[1], sqrt(vcov(fit)["heckman::outcome:x", "heckman::outcome:x"]),
    tolerance = 1e-5)
})

test_that("two-step Heckman models are refused", {
  skip_if_not_installed("sampleSelection")
  d <- heckman_test_data()
  two_step <- sampleSelection::selection(s ~ x + w, y ~ x + z, data = d, method = "2step")
  expect_error(suest(two_step, lm(yc ~ x, d)), "maximum-likelihood")
})
