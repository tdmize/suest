gologit_test_data <- function() {
  set.seed(1010)
  n <- 500L
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  d$yo <- ordered(cut(0.8*d$x + 0.3*d$z + rlogis(n), c(-Inf, -1, 0.5, 2, Inf),
    labels = c("low", "mid", "high", "top")))
  d$yb <- rbinom(n, 1L, plogis(0.5*d$x))
  d
}

# Independent log likelihood for yo ~ x + z with x non-parallel and z parallel:
# P(Y > j) = plogis(a_j + b_j x + c z).
gologit_hand_loglik <- function(theta, d) {
  a <- theta[1:3]; b <- theta[4:6]; c <- theta[7]
  upper <- sapply(1:3, function(j) stats::plogis(a[j] + b[j]*d$x + c*d$z))
  p <- cbind(1, upper) - cbind(upper, 0)
  log(p[cbind(seq_len(nrow(d)), as.integer(d$yo))])
}

test_that("generalized ordered logit scores and information match an independent likelihood", {
  skip_if_not_installed("VGAM")
  d <- gologit_test_data()
  fit <- VGAM::vglm(yo ~ x + z, VGAM::cumulative(parallel = FALSE ~ x, reverse = TRUE),
    data = d)
  wrapper <- suest:::.suest_wrap_vglm(fit)
  expect_identical(suest:::.suest_vglm_adapter(wrapper)$type, "gologit")
  component <- suest:::.suest_vglm_components(wrapper)
  theta <- component$parameters
  expect_equal(sum(gologit_hand_loglik(theta, d)),
    as.numeric(VGAM::logLik(fit)), tolerance = 1e-10)
  step <- 1e-6
  numerical <- sapply(seq_along(theta), function(k) {
    up <- down <- theta
    up[k] <- up[k] + step
    down[k] <- down[k] - step
    (gologit_hand_loglik(up, d) - gologit_hand_loglik(down, d))/(2*step)
  })
  expect_equal(unname(component$score), unname(numerical), tolerance = 1e-6)
  hessian <- sapply(seq_along(theta), function(k) {
    up <- down <- theta
    up[k] <- up[k] + 1e-4
    down[k] <- down[k] - 1e-4
    g <- function(t) colSums(sapply(seq_along(t), function(m) {
      u <- l <- t; u[m] <- u[m] + 1e-6; l[m] <- l[m] - 1e-6
      (gologit_hand_loglik(u, d) - gologit_hand_loglik(l, d))/2e-6
    }))
    (g(up) - g(down))/2e-4
  })
  expect_equal(unname(component$bread/nrow(component$score)),
    unname(solve(-(hessian + t(hessian))/2)), tolerance = 1e-4)
})

test_that("generalized ordered predictions are category probabilities", {
  skip_if_not_installed("VGAM")
  d <- gologit_test_data()
  partial <- VGAM::vglm(yo ~ x + z, VGAM::cumulative(parallel = FALSE ~ x, reverse = TRUE),
    data = d)
  probit <- VGAM::vglm(yo ~ x + z, VGAM::cumulative(link = "probitlink", parallel = TRUE),
    data = d)
  fit <- suest(partial, probit, model_names = c("gologit", "oprobit"))
  expect_identical(fit$comparison_scale, "category probabilities")
  predicted <- marginaleffects::predictions(fit, newdata = d[1:2, ])
  expect_equal(predicted$estimate[startsWith(as.character(predicted$group), "gologit::")],
    as.vector(VGAM::fitted(partial)[1:2, ]), tolerance = 1e-10)
  expect_equal(predicted$estimate[startsWith(as.character(predicted$group), "oprobit::")],
    as.vector(VGAM::fitted(probit)[1:2, ]), tolerance = 1e-10)
  effects <- marginaleffects::avg_slopes(fit, variables = "x")
  expect_identical(nrow(effects), 8L)
  expect_equal(sum(effects$estimate[1:4]), 0, tolerance = 1e-8)
})

test_that("a parallel vglm equals MASS::polr in a joint system", {
  skip_if_not_installed("VGAM")
  d <- gologit_test_data()
  partner <- glm(yb ~ x, binomial, d)
  polr <- MASS::polr(yo ~ x + z, d, Hess = TRUE)
  vglm <- VGAM::vglm(yo ~ x + z, VGAM::cumulative(parallel = TRUE, reverse = TRUE), data = d)
  a <- suest(polr, partner, model_names = c("m", "p"))
  b <- suest(vglm, partner, model_names = c("m", "p"))
  expect_equal(unname(coef(a)[c("m::x", "m::z")]), unname(coef(b)[c("m::x", "m::z")]),
    tolerance = 1e-5)
  keep_a <- c("m::x", "m::z", "p::(Intercept)", "p::x")
  expect_equal(unname(vcov(a)[keep_a, keep_a]), unname(vcov(b)[keep_a, keep_a]),
    tolerance = 1e-4)
})
