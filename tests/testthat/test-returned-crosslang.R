# High-precision full matrices from the user-returned Stata 19.5 logs.
test_that("panel FE, BE, and balanced RE match returned Stata matrices", {
  skip_if_not_installed("plm")
  fixture <- readRDS(test_path("fixtures", "returned-crosslang.rds"))
  dat <- fixture$panel
  dat$time <- factor(dat$time)
  for (route in c("within", "between", "random")) {
    for (partial in c(FALSE, TRUE)) {
      if (route == "random" && partial) next # upstream singular between design
      ys <- if (partial) c("y1_partial", "y2_partial") else c("y1", "y2")
      rhs <- if (route == "between") c("x", "z") else c("x", "z", "time")
      models <- lapply(ys, function(y) plm::plm(reformulate(rhs, y), data = dat,
        index = c("id", "time"), model = route))
      fit <- suest(models[[1]], models[[2]], model_names = c("First", "Second"),
        cluster = if (partial && route == "between") "higher" else NULL)
      local <- c("x", "z", if (route != "between") paste0("time", 2:5), "(Intercept)")
      order <- unlist(lapply(c("First", "Second"), function(x) paste0(x, "::", local)))
      ref <- fixture$reference[[paste(route, partial)]]
      expect_lt(max(abs(coef(fit)[order] - ref$b)), 2e-8)
      expect_lt(max(abs(vcov(fit)[order, order] - ref$V)), 3e-9)
      if (route == "within" && partial) {
        fit <- suest(models[[1]], models[[2]], model_names = c("First", "Second"), cluster = "higher")
        expect_lt(max(abs(vcov(fit)[order, order] - fixture$reference[["within higher"]]$V)), 3e-9)
      }
    }
  }
})

test_that("unbalanced RE remains within the documented cross-engine tolerance", {
  skip_if_not_installed("plm")
  fixture <- readRDS(test_path("fixtures", "returned-crosslang.rds"))
  dat <- fixture$panel
  dat$time <- factor(dat$time)
  models <- lapply(c("y1_partial", "y2_partial"), function(y)
    plm::plm(reformulate(c("x", "z"), y), data = dat,
      index = c("id", "time"), model = "random"))
  fit <- suest(models[[1]], models[[2]], model_names = c("First", "Second"),
    cluster = "higher")
  order <- unlist(lapply(c("First", "Second"), function(x)
    paste0(x, "::", c("x", "z", "(Intercept)"))))
  ref <- fixture$reference[["random unbalanced simple"]]
  expect_lt(max(abs(coef(fit)[order] - ref$b)), 2e-4)
  expect_lt(max(abs(vcov(fit)[order, order] - ref$V)), 1e-5)
  expect_lt(max(abs(diag(vcov(fit)[order, order])/diag(ref$V) - 1)), 0.002)
})

test_that("panel ML matches returned Stata within native engine precision", {
  skip_if_not_installed("nlme")
  fixture <- readRDS(test_path("fixtures", "returned-crosslang.rds"))
  dat <- fixture$panel
  dat$time <- factor(dat$time)
  for (partial in c(FALSE, TRUE)) {
    ys <- if (partial) c("y1_partial", "y2_partial") else c("y1", "y2")
    models <- lapply(ys, function(y) nlme::lme(
      reformulate(c("x", "z", "time"), y), random = ~1 | id, data = dat,
      method = "ML", na.action = na.omit))
    fit <- suest(models[[1]], models[[2]], model_names = c("First", "Second"),
      cluster = if (partial) "higher" else NULL)
    local <- c("x", "z", paste0("time", 2:5), "(Intercept)", "sigma_u", "sigma_e")
    order <- unlist(lapply(c("First", "Second"), function(x) paste0(x, "::", local)))
    ref <- fixture$reference[[paste("ML", partial)]]
    expect_lt(max(abs(coef(fit)[order] - ref$b)), 1e-8)
    expect_lt(max(abs(vcov(fit)[order, order] - ref$V)), 1e-6)
    expect_lt(max(abs(diag(vcov(fit)[order, order])/diag(ref$V) - 1)), 1e-4)
  }
})

test_that("ten GEE cases match returned Stata full covariance matrices", {
  skip_if_not_installed("geepack")
  fixture <- readRDS(test_path("fixtures", "returned-crosslang.rds"))
  dat <- fixture$gee
  families <- list(gaussian(), binomial(), binomial("probit"), binomial("cloglog"), poisson())
  outcomes <- c("y", "binary", "binary", "binary", "count")
  order <- c("Base::x", "Base::(Intercept)", "Adjusted::x", "Adjusted::z", "Adjusted::(Intercept)")
  for (correlation in c("independence", "exchangeable")) {
    for (j in seq_along(families)) {
      models <- lapply(list("x", c("x", "z")), function(rhs)
        geepack::geeglm(reformulate(rhs, outcomes[j]), id = id, data = dat,
          family = families[[j]], corstr = correlation,
          control = geepack::geese.control(epsilon = 1e-10, maxit = 100)))
      fit <- suest(models[[1]], models[[2]], model_names = c("Base", "Adjusted"))
      ref <- fixture$reference[[paste("GEE", correlation, outcomes[j], families[[j]]$link)]]
      expect_lt(max(abs(coef(fit)[order] - ref$b)), 1e-7)
      expect_lt(max(abs(vcov(fit)[order, order] - ref$V)), 1e-9)
    }
  }
})

test_that("biprobit curvature correction matches the returned Stata covariance", {
  skip_if_not_installed("mvProbit")
  skip_if_not_installed("mvtnorm")
  fixture <- readRDS(test_path("fixtures", "returned-crosslang.rds"))
  # Retain fitted objects to avoid re-running the slow optimizer on every test.
  fit <- do.call(suest, fixture$biprobit_models)
  order <- c(2, 1, 4, 3, 5, 7, 8, 6, 10, 11, 9, 12)
  ref <- fixture$reference$biprobit
  expect_lt(max(abs(coef(fit)[order] - ref$b)), 3e-8)
  expect_lt(max(abs(vcov(fit)[order, order] - ref$V)), 1e-9)
})

test_that("IV-probit curvature correction matches returned Stata full covariance", {
  skip_if_not_installed("Rchoice")
  fixture <- readRDS(test_path("fixtures", "returned-crosslang.rds"))
  dat <- fixture$endogenous
  first <- Rchoice::ivpml(binary_iv ~ x + endogenous | x + z1, data = dat,
    messages = FALSE, printLevel = 0, reltol = 1e-12)
  second <- Rchoice::ivpml(binary_iv ~ x + endogenous | x + z1 + z2, data = dat,
    messages = FALSE, printLevel = 0, reltol = 1e-12)
  fit <- suest(first, second)
  order <- c(3, 2, 1, 5, 6, 4, 8, 7, 11, 10, 9, 13, 14, 15, 12, 17, 16)
  ref <- fixture$reference$ivprobit
  expect_lt(max(abs(coef(fit)[order] - ref$b)), 1e-6)
  expect_lt(max(abs(vcov(fit)[order, order] - ref$V)), 1e-8)
})

test_that("pweight ancillary scaling reproduces the returned Stata diagnostic", {
  fixture <- readRDS(test_path("fixtures", "returned-crosslang.rds"))
  dat <- fixture$cluster
  fit_scale <- function(scale) {
    w <- dat$pw * scale
    first <- lm(y ~ x + z, data = dat, weights = w)
    second <- suppressWarnings(glm(b ~ x + z, data = dat, weights = w, family = binomial()))
    suest(first, second, weight_type = "pweight", cluster = "cluster4")
  }
  first <- fit_scale(1)
  second <- fit_scale(10)
  expect_lt(abs(unname(coef(second)[4] - coef(first)[4]) -
    -0.00222088538830051), 2e-12)
  expect_lt(abs(unname(vcov(second)[4, 4] - vcov(first)[4, 4]) -
    7.04424779679382e-6), 3e-13)
})
