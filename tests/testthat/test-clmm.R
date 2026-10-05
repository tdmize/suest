clmm_test_data <- function() {
  set.seed(1012)
  groups <- 60L
  size <- 6L
  d <- data.frame(id = factor(rep(seq_len(groups), each = size)),
    x = rnorm(groups*size), z = rnorm(groups*size))
  u <- rep(rnorm(groups, sd = 0.8), each = size)
  d$yo <- ordered(cut(0.7*d$x + 0.3*d$z + u + rlogis(nrow(d)), c(-Inf, -1, 0.5, 2, Inf),
    labels = c("a", "b", "c", "d")))
  d$yb <- rbinom(nrow(d), 1L, plogis(-0.3 + 0.5*d$x + u))
  d
}

test_that("clmm per-group likelihoods reproduce ordinal's log likelihood", {
  skip_if_not_installed("ordinal")
  d <- clmm_test_data()
  for (q in c(1L, 7L)) {
    model <- ordinal::clmm(yo ~ x + z + (1 | id), data = d, nAGQ = q)
    loglik <- suest:::.suest_clmm_group_loglik(model)
    parameters <- suest:::.suest_clmm_parameters(model)
    # clmm's Laplace stops its inner mode search at gradTol = 1e-4.
    expect_equal(sum(loglik(parameters)), as.numeric(stats::logLik(model)),
      tolerance = if (q == 1L) 1e-5 else 1e-9)
    component <- suest:::.suest_clmm_components(model)
    if (q > 1L) expect_lt(max(abs(colSums(component$score))), 1e-3)
  }
})

test_that("a binary clmm equals the matching glmer in a joint system", {
  skip_if_not_installed("ordinal")
  skip_if_not_installed("lme4")
  d <- clmm_test_data()
  d$fb <- ordered(d$yb)
  ordered_fit <- ordinal::clmm(fb ~ x + (1 | id), data = d, nAGQ = 7)
  binary_fit <- lme4::glmer(yb ~ x + (1 | id), family = binomial, data = d, nAGQ = 7)
  partner <- lm(z ~ x, d)
  a <- suest(ordered_fit, partner, model_names = c("m", "p"))
  b <- suest(binary_fit, partner, model_names = c("m", "p"))
  # P(Y = 1) = F(-threshold + x b + u): the threshold is minus the intercept.
  flip <- diag(c(-1, 1, 1, 1, 1, 1))
  expect_equal(drop(flip %*% unname(coef(a))), unname(coef(b)), tolerance = 1e-4)
  expect_equal(unname(flip %*% vcov(a) %*% flip), unname(vcov(b)), tolerance = 1e-3)
})

test_that("clmm predictions are population-averaged category probabilities", {
  skip_if_not_installed("ordinal")
  d <- clmm_test_data()
  rows <- d[1:2, ]
  for (link in c("logit", "probit")) {
    model <- ordinal::clmm(yo ~ x + z + (1 | id), data = d, nAGQ = 7, link = link)
    fit <- suest(model, glm(yb ~ x, binomial, d), model_names = c("ordered", "logit"))
    predicted <- marginaleffects::predictions(fit, newdata = rows)
    estimate <- predicted$estimate[startsWith(as.character(predicted$group), "ordered::")]
    p <- coef(fit)
    eta <- p[["ordered::x"]]*rows$x + p[["ordered::z"]]*rows$z
    sigma <- exp(p[["ordered::log_sigma"]])
    F <- if (link == "logit") stats::plogis else stats::pnorm
    cumulative <- sapply(c("ordered::a|b", "ordered::b|c", "ordered::c|d"), function(t)
      sapply(eta, function(e) stats::integrate(function(u) F(p[[t]] - e - sigma*u)*stats::dnorm(u),
        -30, 30, rel.tol = 1e-12)$value))
    expected <- cbind(cumulative, 1) - cbind(0, cumulative)
    expect_equal(estimate, as.vector(expected), tolerance = 1e-8, label = link)
    expect_identical(fit$n_clusters, 60L)
  }
})
