re_engine_test_data <- function(groups = 80L, size = 6L) {
  set.seed(1105)
  d <- data.frame(id = factor(rep(seq_len(groups), each = size)),
    x = rnorm(groups*size), z = rnorm(groups*size))
  u <- rep(rnorm(groups, sd = 0.7), each = size)
  u1 <- rep(rnorm(groups, sd = 0.4), each = size)
  d$yb <- rbinom(nrow(d), 1L, plogis(-0.2 + 0.6*d$x + u + u1*d$x))
  d$count <- rpois(nrow(d), exp(0.3 + 0.3*d$x + 0.5*u + u1*d$x))
  d$nb <- rnbinom(nrow(d), size = 2, mu = exp(0.5 + 0.3*d$x + u))
  d$yc <- 1 + 0.5*d$x + u + u1*d$x + rnorm(nrow(d))
  d$positive <- rgamma(nrow(d), shape = 3, scale = exp(0.2 + 0.3*d$x + u)/3)
  d
}

test_that("glmmTMB group scores match glmmTMB's own scores for every family", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  d <- re_engine_test_data()
  fit <- function(...) suppressWarnings(glmmTMB::glmmTMB(..., data = d))
  models <- list(
    fit(yb ~ x + z + (1 | id), family = binomial),
    fit(yb ~ x + (1 | id), family = binomial("probit")),
    fit(yb ~ x + (1 | id), family = binomial("cloglog")),
    fit(count ~ x + (1 | id), family = poisson),
    fit(nb ~ x + (1 | id), family = glmmTMB::nbinom2),
    fit(nb ~ x + (1 | id), family = glmmTMB::nbinom1),
    fit(yc ~ x + (1 | id)),
    fit(positive ~ x + (1 | id), family = Gamma("log")),
    fit(yb ~ x + (1 + x | id), family = binomial),
    fit(count ~ x + (1 + x | id), family = poisson),
    fit(count ~ x + (1 + x | id), family = glmmTMB::nbinom1),
    fit(yc ~ x + (1 + x | id)))
  for (model in models) {
    type <- suest:::.suest_model_adapter(model)$type
    ours <- suest:::.suest_glmmtmb_group_scores(model, type)
    panel <- droplevels(factor(suest:::.suest_panel_id(model)))
    theirs <- sandwich::estfun(model, full = TRUE, cluster = panel, rawnames = FALSE)
    expect_equal(unname(ours[rownames(theirs), ]), unname(theirs), tolerance = 1e-7)
  }
})

test_that("clmm group scores are exact derivatives of the group log likelihood", {
  skip_if_not_installed("ordinal")
  d <- re_engine_test_data()
  d$yo <- ordered(cut(d$yc, c(-Inf, 0, 1, 2, Inf), labels = c("a", "b", "c", "d")))
  for (link in c("logit", "probit")) for (q in c(1L, 7L)) {
    model <- ordinal::clmm(yo ~ x + z + (1 | id), data = d, link = link, nAGQ = q)
    evaluate <- suest:::.suest_clmm_group_loglik(model)
    theta <- suest:::.suest_clmm_parameters(model)
    base <- evaluate(theta, scores = TRUE)
    numeric <- suest:::.suest_re_numeric_scores(evaluate, theta, base)
    expect_equal(base$score, numeric, tolerance = 1e-6)
  }
})

test_that("logit-normal means at large standard deviations match direct integration", {
  set.seed(1106)
  eta <- runif(40, -20, 20)
  sigma <- runif(40, 3.01, 40)
  direct <- vapply(seq_along(eta), function(i) stats::integrate(function(t)
    stats::dlogis(t)*stats::pnorm((eta[i] - t)/sigma[i]), -Inf, Inf,
    rel.tol = 1e-12, abs.tol = 1e-13)$value, numeric(1))
  expect_equal(suest:::.suest_logit_normal_mean(eta, sigma), direct, tolerance = 1e-10)
})
