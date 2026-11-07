twopart_data <- function() {
  set.seed(2026)
  g <- rep(1:60, each = 10)
  u <- rnorm(60, 0, 0.5)[g]
  dat <- data.frame(g = g, x = rnorm(600), z = rbinom(600, 1, 0.5), w = runif(600))
  dat$pos <- rbinom(600, 1, plogis(-0.2 + 0.7 * dat$x + 0.5 * dat$z + u))
  mu <- exp(1 + 0.4 * dat$x + 0.3 * dat$z + 0.6 * dat$w + 0.5 * u)
  dat$y <- ifelse(dat$pos == 1, rgamma(600, shape = 2, rate = 2 / mu), 0)
  dat
}

test_that("two-part predictions and effects match a by-hand delta method", {
  dat <- twopart_data()
  any <- glm(I(y > 0) ~ x + z, family = binomial(), data = dat)
  amount <- glm(y ~ x + z, family = Gamma(link = "log"), data = dat, subset = y > 0)
  fit <- suest_twopart(any, amount, cluster = "g")
  expect_s3_class(fit, c("suest_twopart", "suest_model"))
  expect_equal(fit$model_names, c("any", "amount"))

  b <- coef(fit)
  X <- model.matrix(~ x + z, dat)
  expected <- function(b, M = X) mean(plogis(M %*% b[1:3]) * exp(M %*% b[4:6]))
  plus <- X
  plus[, "x"] <- plus[, "x"] + 1
  effect <- function(b) expected(b, plus) - expected(b)
  se <- function(f) {
    gr <- vapply(seq_along(b), function(j) {
      h <- 1e-6 * max(1, abs(b[j]))
      up <- b; up[j] <- up[j] + h
      dn <- b; dn[j] <- dn[j] - h
      (f(up) - f(dn)) / (2 * h)
    }, numeric(1))
    sqrt(drop(gr %*% vcov(fit) %*% gr))
  }

  p <- marginaleffects::avg_predictions(fit)
  expect_equal(p$estimate, expected(b), tolerance = 1e-10)
  expect_equal(p$std.error, se(expected), tolerance = 1e-5)

  e <- marginaleffects::avg_comparisons(fit, variables = list(x = 1))
  expect_equal(e$estimate, effect(b), tolerance = 1e-10)
  expect_equal(e$std.error, se(effect), tolerance = 1e-5)

  expect_equal(marginaleffects::avg_predictions(fit, type = "prob")$estimate,
               mean(plogis(X %*% b[1:3])), tolerance = 1e-10)
  expect_equal(marginaleffects::avg_predictions(fit, type = "conditional")$estimate,
               mean(exp(X %*% b[4:6])), tolerance = 1e-10)
  expect_setequal(unique(marginaleffects::avg_slopes(fit)$term), c("x", "z"))
  expect_equal(nrow(marginaleffects::avg_comparisons(fit, variables = "x", by = "z")), 2)
  expect_equal(nrow(marginaleffects::predictions(fit,
    newdata = marginaleffects::datagrid(x = c(0, 1)))), 2)
})

test_that("the second part can use its own predictors", {
  dat <- twopart_data()
  any <- glm(I(y > 0) ~ x + factor(z), family = binomial("probit"), data = dat)
  amount <- glm(y ~ x + factor(z) + w, family = quasipoisson(), data = dat, subset = y > 0)
  fit <- suest_twopart(any, amount, model_names = c("Any", "Amount"))
  b <- coef(fit)
  X1 <- model.matrix(~ x + factor(z), dat)
  X2 <- model.matrix(~ x + factor(z) + w, dat)
  expect_equal(marginaleffects::avg_predictions(fit)$estimate,
               mean(pnorm(X1 %*% b[1:3]) * exp(X2 %*% b[4:7])), tolerance = 1e-10)
  expect_setequal(unique(marginaleffects::avg_slopes(fit)$term), c("x", "z", "w"))
})

test_that("weighted and survey two-part models agree", {
  skip_if_not_installed("survey")
  dat <- twopart_data()
  dat$wt <- 0.5 + dat$w
  dat$y <- ceiling(100 * dat$y)
  dat$id <- seq_len(nrow(dat))
  any <- suppressWarnings(glm(I(y > 0) ~ x + z, family = binomial(), data = dat, weights = wt))
  amount <- suppressWarnings(glm(y ~ x + z, family = poisson(), data = dat,
    subset = y > 0, weights = wt))
  weighted <- suest_twopart(any, amount, weight_type = "pweight", cluster = "g")
  design <- survey::svydesign(ids = ~g, weights = ~wt, data = dat)
  s_any <- survey::svyglm(I(y > 0) ~ x + z, design = design, family = quasibinomial())
  s_amount <- survey::svyglm(y ~ x + z, design = subset(design, y > 0), family = quasipoisson())
  surveyed <- suest_twopart(s_any, s_amount, survey_design = design, observation_id = "id")
  a <- marginaleffects::avg_comparisons(weighted, variables = "z", wts = ".suest_weight")
  s <- marginaleffects::avg_comparisons(surveyed, variables = "z", wts = ".suest_weight")
  expect_equal(unname(coef(weighted)), unname(coef(surveyed)), tolerance = 1e-8)
  expect_equal(a$estimate, s$estimate, tolerance = 1e-8)
  expect_equal(a$std.error, s$std.error, tolerance = 1e-6)
})

test_that("suest_twopart() checks the parts", {
  dat <- twopart_data()
  any <- glm(I(y > 0) ~ x, family = binomial(), data = dat)
  amount <- glm(y ~ x, family = Gamma(link = "log"), data = dat, subset = y > 0)
  expect_error(suest_twopart(amount, any), "first part must be a binary")
  expect_error(suest_twopart(any, glm(y ~ x, family = Gamma(link = "log"), data = dat,
    subset = y > 0 & z == 1)), "positive outcome")
  expect_error(suest_twopart(any, lm(x ~ z, data = dat)), "positive outcome")
  fit <- suest_twopart(any, amount)
  expect_error(marginaleffects::avg_predictions(fit, newdata = suest_newdata(fit)),
               "one row per observation")
  expect_error(marginaleffects::avg_predictions(fit, type = "link"), "'type' must be")
  expect_output(print(fit), "Two-part model")
})
