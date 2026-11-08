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

# Independent log likelihood for any link: P(Y > j) = upper(a_j + b_j x + c z).
gologit_link_loglik <- function(theta, d, upper) {
  a <- theta[1:3]; b <- theta[4:6]; c <- theta[7]
  above <- matrix(sapply(1:3, function(j) upper(a[j] + b[j]*d$x + c*d$z)), nrow(d))
  p <- cbind(1, above) - cbind(above, 0)
  log(p[cbind(seq_len(nrow(d)), as.integer(d$y))])
}

test_that("cauchit and log-log generalized ordered models have exact scores", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("numDeriv")
  set.seed(1010)
  n <- 500L
  base <- data.frame(x = rnorm(n), z = rnorm(n))
  labels <- c("low", "mid", "high", "top")
  cases <- list(
    # Stata gologit2 link(cauchit): P(Y > j) = pcauchy(index).
    cauchit = list(family = VGAM::cumulative(link = "cauchitlink",
      parallel = FALSE ~ x, reverse = TRUE), upper = stats::pcauchy, sign = 1,
      y = cut(0.8*base$x + 0.3*base$z + rcauchy(n), c(-Inf, -1, 0.5, 2, Inf))),
    # Stata gologit2 link(loglog): P(Y > j) = exp(-exp(-index)), which is
    # VGAM's cloglog with reverse = FALSE and every coefficient negated.
    loglog = list(family = VGAM::cumulative(link = "clogloglink",
      parallel = FALSE ~ x, reverse = FALSE), upper = function(x) exp(-exp(-x)), sign = -1,
      y = cut(0.3*base$x + 0.2*base$z - log(-log(runif(n))), c(-Inf, -0.8, 0, 1, Inf))))
  for (case in names(cases)) {
    spec <- cases[[case]]
    d <- base
    d$y <- ordered(spec$y, labels = labels)
    fit <- VGAM::vglm(y ~ x + z, spec$family, data = d)
    wrapper <- suest:::.suest_wrap_vglm(fit)
    component <- suest:::.suest_vglm_components(wrapper)
    theta <- component$parameters
    stata <- spec$sign*unname(theta)
    expect_equal(sum(gologit_link_loglik(stata, d, spec$upper)),
      as.numeric(VGAM::logLik(fit)), tolerance = 1e-10, label = case)
    for (i in c(3L, 150L, 444L)) {
      gradient <- numDeriv::grad(function(t) gologit_link_loglik(t, d[i, ], spec$upper), stata)
      expect_equal(unname(component$score[i, ]), spec$sign*gradient, tolerance = 1e-7,
        label = case)
    }
    information <- -numDeriv::hessian(function(t) sum(gologit_link_loglik(t, d, spec$upper)),
      stata)
    expect_equal(unname(component$bread/nrow(d)), unname(solve(information)),
      tolerance = 1e-5, label = case)
    predicted <- suest:::.suest_predict_probabilities(wrapper, d[1:4, ], "VGAM::vglm")
    expect_equal(unname(predicted), unname(VGAM::fitted(fit)[1:4, ]), tolerance = 1e-10,
      label = case)
    partner <- glm(I(as.integer(y) > 2) ~ x, binomial, d)
    joint <- suest(fit, partner, model_names = c(case, "logit"))
    expect_identical(unname(joint$model_types), c("gologit", "logit"))
  }
  d <- base
  d$y <- ordered(cases$cauchit$y, labels = labels)
  other <- VGAM::vglm(y ~ x, VGAM::cumulative(parallel = TRUE), data = d)
  other@misc$link[] <- "foldsqrtlink"
  expect_error(suest(other, glm(I(as.integer(y) > 2) ~ x, binomial, d)), "cauchit link")
})
