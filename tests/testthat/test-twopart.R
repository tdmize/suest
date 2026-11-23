twopart_data <- function() {
  set.seed(2026)
  g <- rep(1:60, each = 10)
  u <- rnorm(60, 0, 0.5)[g]
  dat <- data.frame(g = g, x = rnorm(600), z = rbinom(600, 1, 0.5), w = runif(600))
  dat$pos <- rbinom(600, 1, plogis(-0.2 + 0.7 * dat$x + 0.5 * dat$z + u))
  mu <- exp(1 + 0.4 * dat$x + 0.3 * dat$z + 0.6 * dat$w + 0.5 * u)
  dat$y <- ifelse(dat$pos == 1, rgamma(600, shape = 2, rate = 2 / mu), 0)
  dat$id <- seq_len(nrow(dat))
  dat
}

delta_se <- function(f, b, V) {
  gr <- vapply(seq_along(b), function(j) {
    h <- 1e-6 * max(1, abs(b[j]))
    up <- b; up[j] <- up[j] + h
    dn <- b; dn[j] <- dn[j] - h
    (f(up) - f(dn)) / (2 * h)
  }, numeric(1))
  sqrt(drop(gr %*% V %*% gr))
}

test_that("two-part predictions and effects match a by-hand delta method", {
  dat <- twopart_data()
  any <- glm(I(y > 0) ~ x + z, family = binomial(), data = dat)
  amount <- glm(y ~ x + z, family = Gamma(link = "log"), data = dat, subset = y > 0)
  fit <- suest(any, amount, model_names = c("Any", "Amount"), cluster = "g", twopart = TRUE)
  b <- coef(fit)
  X <- model.matrix(~ x + z, dat)
  expected <- function(b, M = X) mean(plogis(M %*% b[1:3]) * exp(M %*% b[4:6]))
  plus <- X
  plus[, "x"] <- plus[, "x"] + 1
  effect <- function(b) expected(b, plus) - expected(b)

  p <- marginaleffects::avg_predictions(fit)
  expect_false("group" %in% names(p))
  expect_equal(p$estimate, expected(b), tolerance = 1e-10)
  expect_equal(p$std.error, delta_se(expected, b, vcov(fit)), tolerance = 1e-5)
  e <- marginaleffects::avg_comparisons(fit, variables = list(x = 1))
  expect_equal(e$estimate, effect(b), tolerance = 1e-10)
  expect_equal(e$std.error, delta_se(effect, b, vcov(fit)), tolerance = 1e-5)

  expect_equal(marginaleffects::avg_predictions(fit, type = "positive")$estimate,
               mean(plogis(X %*% b[1:3])), tolerance = 1e-10)
  expect_equal(marginaleffects::avg_predictions(fit, type = "amount")$estimate,
               mean(exp(X %*% b[4:6])), tolerance = 1e-10)
  expect_equal(marginaleffects::avg_predictions(fit, type = "amount",
    newdata = subset(y > 0))$estimate, mean(fitted(amount)), tolerance = 1e-10)
  expect_setequal(unique(marginaleffects::avg_slopes(fit)$term), c("x", "z"))
  expect_equal(nrow(marginaleffects::avg_comparisons(fit, variables = "x", by = "z")), 2)
  expect_equal(nrow(marginaleffects::predictions(fit,
    newdata = marginaleffects::datagrid(x = c(0, 1)))), 2)
  expect_equal(nrow(marginaleffects::predictions(fit)), nrow(dat))
  expect_error(marginaleffects::avg_predictions(fit, type = "link"), "'type' must be")
  expect_output(print(fit), "Two-part models: Any\\+Amount")
})

test_that("the second part can use its own predictors and other model types", {
  dat <- twopart_data()
  any <- glm(I(y > 0) ~ x + factor(z), family = binomial("probit"), data = dat)
  amount <- glm(y ~ x + factor(z) + w, family = quasipoisson(), data = dat, subset = y > 0)
  fit <- suest(any, amount, model_names = c("Any", "Amount"), twopart = TRUE)
  b <- coef(fit)
  X1 <- model.matrix(~ x + factor(z), dat)
  X2 <- model.matrix(~ x + factor(z) + w, dat)
  expect_equal(marginaleffects::avg_predictions(fit)$estimate,
               mean(pnorm(X1 %*% b[1:3]) * exp(X2 %*% b[4:7])), tolerance = 1e-10)
  expect_setequal(unique(marginaleffects::avg_slopes(fit)$term), c("x", "z", "w"))

  dat$count <- round(10 * dat$y)
  any <- glm(I(count > 0) ~ x + z, family = binomial("cloglog"), data = dat)
  nb <- MASS::glm.nb(count ~ x + z, data = dat, subset = count > 0)
  fit <- suest(any, nb, model_names = c("Any", "NegBin"), twopart = TRUE)
  expect_equal(marginaleffects::avg_predictions(fit)$estimate,
               mean(fitted(any) * predict(nb, newdata = dat, type = "response")), tolerance = 1e-10)
  lin <- lm(y ~ x + z, data = dat, subset = y > 0)
  fit <- suest(glm(I(y > 0) ~ x + z, family = binomial(), data = dat), lin,
               model_names = c("Any", "Linear"), twopart = TRUE)
  expect_equal(marginaleffects::avg_predictions(fit, type = "amount")$estimate,
               mean(predict(lin, newdata = dat)), tolerance = 1e-10)
})

test_that("a log-outcome second part is retransformed", {
  dat <- twopart_data()
  dat$lny <- ifelse(dat$y > 0, log(dat$y), NA)
  any <- glm(I(y > 0) ~ x + z, family = binomial(), data = dat)
  ln <- lm(lny ~ x + z, data = dat, subset = y > 0)
  p <- fitted(any)
  xb <- predict(ln, newdata = dat)
  duan <- suest(any, ln, model_names = c("Any", "LnAmount"), twopart = TRUE, retransform = "duan")
  expect_equal(marginaleffects::avg_predictions(duan)$estimate,
               mean(p * exp(xb)) * mean(exp(residuals(ln))), tolerance = 1e-10)
  normal <- suest(any, ln, model_names = c("Any", "LnAmount"), twopart = TRUE,
                  retransform = "normal")
  s2 <- mean(residuals(ln)^2)
  expect_equal(marginaleffects::avg_predictions(normal)$estimate,
               mean(p * exp(xb + s2 / 2)), tolerance = 1e-10)
  X <- model.matrix(~ x + z, dat)
  b <- coef(normal)
  f <- function(b) mean(plogis(X %*% b[1:3]) * exp(X %*% b[4:6] + s2 / 2))
  expect_equal(marginaleffects::avg_predictions(normal)$std.error,
               delta_se(f, b, vcov(normal)), tolerance = 1e-5)
  smear <- mean(exp(residuals(ln)))
  f <- function(b) mean(plogis(X %*% b[1:3]) * exp(X %*% b[4:6])) * smear
  expect_equal(marginaleffects::avg_predictions(duan)$std.error,
               delta_se(f, b, vcov(duan)), tolerance = 1e-5)
  expect_output(print(duan), "Retransformation: duan")
  expect_error(suest(any, glm(y ~ x, family = Gamma(link = "log"), data = dat, subset = y > 0),
    twopart = TRUE, retransform = "duan"), "lm\\(\\) second parts")
  expect_error(suest(any, ln, retransform = "duan"), "specify 'twopart'")
})

test_that("pairs can share a first part and sit beside ordinary models", {
  dat <- twopart_data()
  any <- glm(I(y > 0) ~ x + z, family = binomial(), data = dat)
  gam <- glm(y ~ x + z, family = Gamma(link = "log"), data = dat, subset = y > 0)
  pois <- glm(y ~ x + z, family = quasipoisson(), data = dat, subset = y > 0)
  one <- glm(y ~ x + z, family = quasipoisson(), data = dat)
  fit <- suest(any, gam, pois, model_names = c("Any", "Gamma", "Poisson"),
               twopart = c("Any", "Gamma", "Any", "Poisson"))
  e <- marginaleffects::avg_comparisons(fit, variables = "z")
  expect_equal(as.character(e$group), c("Any+Gamma", "Any+Poisson"))
  same <- suest(any, gam, any, pois, model_names = c("A1", "Gamma", "A2", "Poisson"),
                twopart = TRUE)
  expect_equal(marginaleffects::avg_comparisons(same, variables = "z")$estimate, e$estimate,
               tolerance = 1e-10)

  fit <- suest(any, pois, one, model_names = c("Any", "Poisson", "OnePart"),
               twopart = list("Two-part" = c("Any", "Poisson")))
  e <- marginaleffects::avg_comparisons(fit, variables = "z")
  expect_equal(as.character(e$group), c("Two-part", "OnePart"))
  alone <- marginaleffects::avg_comparisons(one, variables = "z", vcov = "HC0")
  expect_equal(e$estimate[2], alone$estimate, tolerance = 1e-10)
  h <- marginaleffects::hypotheses(e, hypothesis = difference ~ pairwise)
  expect_equal(nrow(h), 1)
  expect_message(a <- marginaleffects::avg_comparisons(fit, variables = "z", type = "amount"),
                 "left out: OnePart")
  expect_equal(as.character(a$group), "Two-part")

  # A model on a different sample: one block of rows per model
  part <- glm(y ~ x, family = quasipoisson(), data = dat, subset = z == 1)
  fit <- suest(any, pois, part, model_names = c("Any", "Poisson", "Subset"),
               twopart = list(TP = c("Any", "Poisson")))
  nd <- suest_newdata(fit)
  expect_equal(as.vector(table(nd$.suest_model)[c("TP", "Subset")]), c(600, sum(dat$z == 1)))
  expect_message(p <- marginaleffects::avg_predictions(fit, type = "positive"), "left out")
  expect_equal(p$estimate, mean(fitted(any)), tolerance = 1e-10)
  g <- marginaleffects::predictions(fit, newdata = marginaleffects::datagrid(x = 0,
    .suest_model = c("TP", "Subset")))
  expect_equal(nrow(g), 2)
})

test_that("twopart checks the pairs", {
  dat <- twopart_data()
  any <- glm(I(y > 0) ~ x, family = binomial(), data = dat)
  amount <- glm(y ~ x, family = Gamma(link = "log"), data = dat, subset = y > 0)
  expect_error(suest(amount, any, twopart = TRUE), "first part of 'amount\\+any' must be a binary")
  expect_error(suest(any, glm(y ~ x, family = Gamma(link = "log"), data = dat,
    subset = y > 0 & z == 1), twopart = TRUE), "positive outcome")
  expect_error(suest(any, amount, twopart = "any"), "in pairs")
  expect_error(suest(any, amount, twopart = c("any", "other")), "not in the system")
  expect_error(suest(any, amount, twopart = 3), "must be TRUE")
  dat$k <- ordered(dat$z + rbinom(600, 1, 0.5))
  any <- glm(I(y > 0) ~ x, family = binomial(), data = dat)
  amount <- glm(y ~ x, family = Gamma(link = "log"), data = dat, subset = y > 0)
  ordered <- MASS::polr(k ~ x, data = dat, Hess = TRUE)
  expect_error(suest(any, amount, ordered, twopart = list(c("any", "amount"))),
               "single outcome")
})

test_that("weighted and survey two-part models agree, and wts = TRUE uses the weights", {
  skip_if_not_installed("survey")
  dat <- twopart_data()
  dat$wt <- 0.5 + dat$w
  any <- glm(I(y > 0) ~ x + z, family = quasibinomial(), data = dat, weights = wt)
  amount <- glm(y ~ x + z, family = quasipoisson(), data = dat, subset = y > 0, weights = wt)
  weighted <- suest(any, amount, weight_type = "pweight", cluster = "g", twopart = TRUE)
  design <- survey::svydesign(ids = ~g, weights = ~wt, data = dat)
  s_any <- survey::svyglm(I(y > 0) ~ x + z, design = design, family = quasibinomial())
  s_amount <- survey::svyglm(y ~ x + z, design = subset(design, y > 0), family = quasipoisson())
  surveyed <- suest(s_any, s_amount, survey_design = design, observation_id = "id",
                    twopart = TRUE)
  a <- marginaleffects::avg_comparisons(weighted, variables = "z", wts = TRUE)
  s <- marginaleffects::avg_comparisons(surveyed, variables = "z", wts = TRUE)
  expect_equal(unname(coef(weighted)), unname(coef(surveyed)), tolerance = 1e-8)
  expect_equal(a$estimate, s$estimate, tolerance = 1e-8)
  expect_equal(a$std.error, s$std.error, tolerance = 1e-6)
  b <- coef(weighted)
  X <- model.matrix(~ x + z, dat)
  ey <- plogis(X %*% b[1:3]) * exp(X %*% b[4:6])
  expect_equal(marginaleffects::avg_predictions(weighted, wts = TRUE)$estimate,
               weighted.mean(ey, dat$wt), tolerance = 1e-10)

  gam <- glm(y ~ x + z, family = Gamma(link = "log"), data = dat, subset = y > 0, weights = wt)
  fit <- suest(any, gam, weight_type = "pweight", twopart = TRUE)
  expect_equal(marginaleffects::avg_predictions(fit, wts = TRUE)$estimate,
               weighted.mean(plogis(X %*% coef(any)) * exp(X %*% coef(gam)), dat$wt),
               tolerance = 1e-10)
  n <- nrow(dat)
  qp <- vcov(weighted)[4:6, 4:6]
  expect_equal(unname(qp), unname(sandwich::vcovCL(amount, cluster = ~g, type = "HC0",
    cadjust = FALSE) * 60 / 59), tolerance = 1e-8)
})

test_that("get_vcov() honors vcov = FALSE and says when it ignores vcov", {
  dat <- twopart_data()
  any <- glm(I(y > 0) ~ x, family = binomial(), data = dat)
  amount <- glm(y ~ x, family = Gamma(link = "log"), data = dat, subset = y > 0)
  fit <- suest(any, amount, twopart = TRUE)
  expect_message(marginaleffects::avg_predictions(fit, vcov = "HC3"), "is ignored")
  expect_false("std.error" %in% names(marginaleffects::avg_predictions(fit, vcov = FALSE)))
})
