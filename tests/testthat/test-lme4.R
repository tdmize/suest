lme4_test_data <- function(groups = 60L, size = 6L) {
  set.seed(1008)
  d <- data.frame(id = factor(rep(seq_len(groups), each = size)),
    time = rep(seq_len(size), groups), x = rnorm(groups*size), z = rnorm(groups*size))
  u <- rep(rnorm(groups, sd = 0.7), each = size)
  u1 <- rep(rnorm(groups, sd = 0.4), each = size)
  d$yb <- rbinom(nrow(d), 1L, plogis(-0.2 + 0.5*d$x + u))
  d$yc <- 1 + 0.5*d$x - 0.3*d$z + u + u1*d$x + rnorm(nrow(d))
  d$count <- rpois(nrow(d), exp(0.2 + 0.3*d$x + 0.5*u))
  d
}

lme4_block <- function(fit, model) {
  keep <- startsWith(names(coef(fit)), paste0(model, "::"))
  vcov(fit)[keep, keep, drop = FALSE]
}

test_that("lmer systems match the same models fit by glmmTMB", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  d <- lme4_test_data()
  lmer_fit <- lme4::lmer(yc ~ x + z + (1 | id), data = d, REML = FALSE)
  tmb_fit <- glmmTMB::glmmTMB(yc ~ x + z + (1 | id), data = d)
  partner <- glm(yb ~ x, binomial, d)
  a <- suest(lmer_fit, partner, model_names = c("mixed", "logit"))
  b <- suest(tmb_fit, partner, model_names = c("mixed", "logit"))
  expect_identical(a$model_types[["mixed"]], "glmm_gaussian_ri")
  expect_identical(names(coef(a)), names(coef(b)))
  expect_equal(coef(a), coef(b), tolerance = 1e-5)
  expect_equal(vcov(a), vcov(b), tolerance = 1e-4)

  slope <- lme4::lmer(yc ~ x + z + (1 + x | id), data = d, REML = FALSE)
  slope_tmb <- glmmTMB::glmmTMB(yc ~ x + z + (1 + x | id), data = d)
  expect_equal(vcov(suest(slope, partner, model_names = c("m", "p"))),
    vcov(suest(slope_tmb, partner, model_names = c("m", "p"))), tolerance = 1e-3)
})

test_that("lmer group scores match merDeriv's analytic scores", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("merDeriv")
  d <- lme4_test_data()
  fit <- lme4::lmer(yc ~ x + z + (1 | id), data = d, REML = FALSE)
  wrapper <- suest:::.suest_wrap_lme4(fit)
  component <- suest:::.suest_lme4_components(wrapper, "glmm_gaussian_ri")
  ours <- rowsum(component$score, d$id)
  analytic <- merDeriv::estfun.lmerMod(fit, level = 2)
  # merDeriv uses the variances; d variance/d log SD = 2*variance.
  variance <- c(exp(2*component$parameters[["log_sigma"]]),
    exp(2*component$parameters[["log_sigma_e"]]))
  converted <- cbind(analytic[, 1:3], analytic[, 5]*2*variance[2], analytic[, 4]*2*variance[1])
  expect_equal(unname(ours), unname(converted), tolerance = 1e-6)
})

test_that("glmer scores come from the fit's own quadrature", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  d <- lme4_test_data()
  adaptive <- lme4::glmer(yb ~ x + (1 | id), family = binomial, data = d, nAGQ = 7)
  wrapper <- suest:::.suest_wrap_lme4(adaptive)
  loglik <- suest:::.suest_lme4_group_loglik(wrapper, "glmm_logit_ri")
  parameters <- suest:::.suest_lme4_parameters(wrapper, "glmm_logit_ri")
  expect_equal(sum(loglik(parameters)), as.numeric(stats::logLik(adaptive)),
    tolerance = 1e-9)
  component <- suest:::.suest_lme4_components(wrapper, "glmm_logit_ri")
  expect_lt(max(abs(colSums(component$score))), 1e-3)
  model_based <- component$bread/nrow(component$score)
  expect_equal(unname(model_based[1:2, 1:2]), unname(as.matrix(vcov(adaptive))),
    tolerance = 1e-3)

  laplace <- lme4::glmer(yb ~ x + (1 | id), family = binomial, data = d,
    control = lme4::glmerControl(tolPwrss = 1e-12))
  tmb <- glmmTMB::glmmTMB(yb ~ x + (1 | id), family = binomial, data = d)
  partner <- lm(yc ~ x, d)
  expect_equal(vcov(suest(laplace, partner, model_names = c("m", "p"))),
    vcov(suest(tmb, partner, model_names = c("m", "p"))), tolerance = 1e-4)
})

test_that("lme4 models combine with other models and predict population averages", {
  skip_if_not_installed("lme4")
  d <- lme4_test_data()
  logit <- lme4::glmer(yb ~ x + (1 | id), family = binomial, data = d, nAGQ = 7)
  count <- lme4::glmer(count ~ x + (1 | id), family = poisson, data = d, nAGQ = 7)
  fit <- suest(logit, count, model_names = c("melogit", "mepoisson"))
  expect_identical(fit$n_clusters, 60L)
  rows <- d[1:3, ]
  p <- coef(fit)
  sigma <- exp(p[["mepoisson::log_sigma"]])
  expected <- exp(p[["mepoisson::(Intercept)"]] + p[["mepoisson::x"]]*rows$x + sigma^2/2)
  predicted <- marginaleffects::predictions(fit, newdata = rows)
  expect_equal(predicted$estimate[predicted$group == "mepoisson"], expected, tolerance = 1e-10)
  effects <- marginaleffects::avg_comparisons(fit, variables = "x")
  expect_true(all(is.finite(effects$std.error)))
})

test_that("unsupported lme4 models are refused with a reason", {
  skip_if_not_installed("lme4")
  d <- lme4_test_data()
  reml <- lme4::lmer(yc ~ x + (1 | id), data = d)
  expect_error(suest(reml, reml), "REML = FALSE")
  d$school <- factor((as.integer(d$id) - 1L) %/% 5L)
  nested <- suppressMessages(lme4::lmer(yc ~ x + (1 | school/id), data = d, REML = FALSE))
  expect_error(suest(nested, nested), "one grouping variable")
  d$positive <- exp(d$z)
  gamma <- suppressWarnings(lme4::glmer(positive ~ x + (1 | id),
    family = Gamma("log"), data = d))
  expect_error(suest(gamma, gamma), "glmer support is currently limited")
})

test_that("lmer and Gaussian glmmTMB use Stata mixed's covariance layout", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  set.seed(2207)
  d <- data.frame(id = rep(seq_len(40), each = 6), x = rnorm(240))
  u <- rep(rnorm(40), each = 6)
  d$y <- 1 + 0.5*d$x + u + rnorm(240)
  d$y2 <- -0.5 + 0.3*d$x + 0.7*u + rnorm(240)
  first <- lme4::lmer(y ~ x + (1 | id), data = d, REML = FALSE)
  second <- lme4::lmer(y2 ~ x + (1 | id), data = d, REML = FALSE)
  wrapper <- suest:::.suest_wrap_lme4(first)
  component <- suest:::.suest_model_components(wrapper, "glmm_gaussian_ri", "lme4::lmer")
  bread <- component$bread/nrow(component$score)
  # Fixed effects: (X'V^-1 X)^-1, which is lmer's own ML covariance.
  expect_equal(unname(bread[1:2, 1:2]), unname(as.matrix(stats::vcov(first))),
    tolerance = 1e-6)
  expect_true(all(bread[1:2, 3:4] == 0))

  fit <- suest(first, second, model_names = c("A", "B"))
  tmb <- suest(glmmTMB::glmmTMB(y ~ x + (1 | id), data = d),
    glmmTMB::glmmTMB(y2 ~ x + (1 | id), data = d), model_names = c("A", "B"))
  scale <- sqrt(outer(diag(vcov(fit)), diag(vcov(fit))))
  expect_lt(max(abs(vcov(fit) - vcov(tmb))/scale), 1e-4)
  # Fixed-effect block: the usual cluster-robust sandwich around lmer's covariance.
  groups <- nlevels(factor(d$id))
  score <- rowsum(component$score, d$id)[, 1:2]
  robust <- groups/(groups - 1)*bread[1:2, 1:2] %*% crossprod(score) %*% bread[1:2, 1:2]
  expect_equal(unname(vcov(fit)[1:2, 1:2]), unname(robust), tolerance = 1e-8)
})
