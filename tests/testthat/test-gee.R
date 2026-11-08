gee_test_data <- function() {
  set.seed(7629)
  id <- rep(1:80, each = 4)
  x <- rnorm(length(id))
  z <- rnorm(length(id))
  u <- rnorm(80, sd = 0.5)[id]
  data.frame(id, time = rep(1:4, 80), higher = ceiling(id / 4), x, z,
    y = 0.4 + 0.6*x - 0.3*z + u + rnorm(length(id)),
    binary = rbinom(length(id), 1, plogis(-0.2 + 0.4*x + u)),
    count = rpois(length(id), exp(0.2 + 0.3*x + 0.3*u)))
}

test_that("GEE reconstructed bread and grouped scores match geepack", {
  skip_if_not_installed("geepack")
  dat <- gee_test_data()
  families <- list(gaussian(), binomial(), binomial("probit"),
    binomial("cloglog"), poisson())
  responses <- c("y", "binary", "binary", "binary", "count")
  for (correlation in c("independence", "exchangeable", "ar1", "unstructured")) {
    for (j in seq_along(families)) {
      formula <- reformulate(c("x", "z"), responses[j])
      model <- geepack::geeglm(formula, id = id, data = dat,
        family = families[[j]], corstr = correlation,
        control = geepack::geese.control(epsilon = 1e-10, maxit = 100))
      component <- suest:::.suest_gee_components(model)
      B <- component$bread / nrow(dat)
      S <- rowsum(component$score, dat$id)
      V <- B %*% crossprod(S) %*% B
      expect_equal(unname(B), unname(model$geese$vbeta.naiv), tolerance = 1e-8)
      expect_equal(unname(V), unname(vcov(model)), tolerance = 1e-8)
      fit <- suest(model, model)
      expect_equal(unname(vcov(fit)[1:3, 1:3]), unname(V * 80/79), tolerance = 1e-8)
      expect_equal(vcov(fit)[1:3, 1:3], vcov(fit)[1:3, 4:6],
        ignore_attr = TRUE, tolerance = 1e-12)
      prediction <- marginaleffects::predictions(fit, newdata = dat[1:3, ])
      expect_equal(prediction$estimate,
        rep(unname(predict(model, dat[1:3, ], type = "response")), 2), tolerance = 1e-8)
      expect_true(all(is.finite(prediction$std.error)))
      X <- model.matrix(~x + z, dat[1:3, ])
      eta <- drop(X %*% coef(model))
      expected_se <- abs(families[[j]]$mu.eta(eta)) *
        sqrt(rowSums((X %*% (V * 80/79)) * X))
      expect_equal(prediction$std.error, rep(unname(expected_se), 2), tolerance = 1e-6)
    }
  }
})

test_that("GEE partial samples, higher clusters, interactions, and offsets work", {
  skip_if_not_installed("geepack")
  dat <- gee_test_data()
  first <- geepack::geeglm(y ~ x * factor(time) + offset(z), id = id,
    data = dat, subset = id <= 65, corstr = "exchangeable")
  second <- geepack::geeglm(y ~ x + z, id = id, data = dat,
    subset = id >= 16 & time != 4, corstr = "exchangeable")
  fit <- suest(first, second, cluster = "higher")
  expect_equal(fit$nobs_overlap, 150L)
  expect_equal(fit$n_clusters, 20L)
  components <- lapply(list(first, second), suest:::.suest_gee_components)
  influences <- lapply(seq_along(components), function(j) {
    model <- list(first, second)[[j]]
    rows <- as.integer(rownames(model$model))
    group <- dat$higher[rows]
    g <- length(unique(group))
    out <- matrix(0, 20, ncol(components[[j]]$score))
    summed <- rowsum(components[[j]]$score, group)
    out[as.integer(rownames(summed)), ] <- summed %*%
      components[[j]]$bread / length(rows) * sqrt(g/(g - 1))
    out
  })
  expected <- crossprod(do.call(cbind, influences))
  expect_equal(unname(vcov(fit)), expected, tolerance = 1e-10)
  pred <- marginaleffects::predictions(fit, newdata = dat[1:3, ])
  expect_equal(pred$estimate, c(unname(predict(first, dat[1:3, ])),
    unname(predict(second, dat[1:3, ]))), tolerance = 1e-8)
  slopes <- marginaleffects::avg_slopes(fit, variables = "x", newdata = dat)
  expect_true(all(is.finite(slopes$std.error)))
  expect_error(suest(first, second, cluster = "time"), "nested within clusters")
})

test_that("GEE ar1 and unstructured follow geepack's waves", {
  skip_if_not_installed("geepack")
  dat <- gee_test_data()
  control <- geepack::geese.control(epsilon = 1e-10, maxit = 100)
  check <- function(model, panels) {
    component <- suest:::.suest_gee_components(model)
    B <- component$bread/nrow(model$model)
    S <- rowsum(component$score, model$id)
    expect_equal(unname(B), unname(model$geese$vbeta.naiv), tolerance = 1e-8)
    expect_equal(unname(B %*% crossprod(S) %*% B), unname(vcov(model)), tolerance = 1e-8)
    fit <- suest(model, model)
    expect_identical(fit$n_clusters, panels)
  }
  # Interior gaps: ar1 uses the distance between wave codes, not positions.
  gaps <- dat[!(dat$id %% 3 == 0 & dat$time == 2) & !(dat$id %% 5 == 0 & dat$time == 4), ]
  ar1 <- geepack::geeglm(binary ~ x + z, id = id, waves = time, data = gaps,
    family = binomial, corstr = "ar1", control = control)
  check(ar1, 80L)
  waves <- suest:::.suest_gee_waves(ar1)
  expect_identical(waves, as.integer(gaps$time))
  R <- suest:::.suest_gee_correlations(ar1)[["3"]]
  expect_equal(R[1, 2], unname(ar1$geese$alpha)^2, tolerance = 1e-12)
  positions <- geepack::geeglm(binary ~ x + z, id = id, data = gaps,
    family = binomial, corstr = "ar1", control = control)
  check(positions, 80L)
  expect_false(isTRUE(all.equal(unname(ar1$geese$alpha), unname(positions$geese$alpha))))

  # Waves recovered from a transformed expression and a subset.
  spaced <- geepack::geeglm(count ~ x + z, id = id, waves = 10*time, data = gaps,
    subset = id > 10, family = poisson, corstr = "ar1", control = control)
  check(spaced, 70L)

  # Unstructured with panels that end early.
  short <- dat[!(dat$id %% 4 == 0 & dat$time == 4), ]
  unstructured <- geepack::geeglm(y ~ x + z, id = id, waves = time, data = short,
    corstr = "unstructured", control = control)
  check(unstructured, 80L)
  bad <- unstructured
  bad$data$time[bad$data$id == 1 & bad$data$time == 2] <- 4
  expect_error(suest(bad, bad), "without gaps")

  combined <- suest(ar1, unstructured, positions, observation_id = c("id", "time"))
  expect_identical(unname(combined$model_types), rep("panel_gee", 3))
  expect_equal(combined$n_clusters, 80L)
  expect_equal(combined$nobs_overlap, sum(rownames(gaps) %in% rownames(short)))
})

test_that("GEE refuses unsupported inputs before falling through to GLM", {
  skip_if_not_installed("geepack")
  dat <- gee_test_data()
  pairs <- sum(table(dat$id)*(table(dat$id) - 1)/2)
  model <- geepack::geeglm(y ~ x, id = id, data = dat, corstr = "fixed",
    zcor = rep(0.2, pairs))
  expect_error(suest(model, model), "unstructured correlation")
  model <- geepack::geeglm(y ~ x, id = id, data = dat, corstr = "ar1")
  bad <- model
  bad$geese$alpha[] <- 1.2
  expect_error(suest(bad, bad), "inside")
  model <- geepack::geeglm(y ~ x, id = id, data = dat)
  expect_error(suest(model, lm(y ~ x, dat)), "same type")
  bad <- model
  bad$geese$error <- 1L
  expect_error(suest(bad, bad), "did not converge")
  bad <- model
  bad$prior.weights[1] <- 2
  expect_error(suest(bad, bad), "Weighted GEE")
  bad <- model
  bad$id[1] <- 80L
  expect_error(suest(bad, bad), "contiguous panels")
})
