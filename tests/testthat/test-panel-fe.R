panel_fe_test_data <- function(groups = 120L, periods = 5L) {
  set.seed(7621)
  id <- rep(seq_len(groups), each = periods)
  time <- rep(seq_len(periods), groups)
  higher <- ceiling(id / 4)
  person1 <- rnorm(groups)[id]
  person2 <- rnorm(groups)[id]
  common <- rnorm(max(higher))[higher]
  x <- rnorm(length(id)) + 0.2*person1
  z <- rnorm(length(id))
  y1 <- 1.1 + 0.55*x - 0.30*z + 0.15*time + person1 + 0.4*common + rnorm(length(id))
  y2 <- -0.4 + 0.35*x + 0.20*z - 0.08*time + person2 + 0.3*person1 + rnorm(length(id))
  data.frame(id, time, higher, x, z, y1, y2)
}

panel_fe_reference <- function(models, cluster) {
  rows <- lapply(models, function(model) rownames(stats::model.frame(model)))
  union_rows <- unique(unlist(rows, use.names = FALSE))
  union_cluster <- as.character(cluster[match(union_rows, names(cluster))])

  transformed <- lapply(seq_along(models), function(i) {
    model <- models[[i]]
    X <- stats::model.matrix(model, model = "within")
    beta <- stats::coef(model)
    X <- X[, names(beta), drop = FALSE]
    X_pooling <- stats::model.matrix(model, model = "pooling")
    means <- colMeans(X_pooling[, names(beta), drop = FALSE])
    mapping <- rbind(-means, diag(length(beta)))
    inverse <- solve(crossprod(X))
    influence <- (X * as.numeric(stats::residuals(model))) %*%
      inverse %*% t(mapping)

    model_cluster <- cluster[match(rows[[i]], names(cluster))]
    g <- length(unique(model_cluster))
    n <- nrow(X)
    correction <- (g / (g - 1)) * ((n - 1) / (n - length(beta) - 1))
    aligned <- matrix(0, nrow = length(union_rows), ncol = ncol(influence))
    aligned[match(rows[[i]], union_rows), ] <- influence * sqrt(correction)
    aligned
  })

  influence <- do.call(cbind, transformed)
  clustered <- rowsum(influence, union_cluster, reorder = FALSE)
  V <- crossprod(clustered)
  (V + t(V)) / 2
}

test_that("plm within systems match the xtreg FE sandwich construction", {
  skip_if_not_installed("plm")
  dat <- panel_fe_test_data()
  first <- plm::plm(
    y1 ~ x + z + factor(time), data = dat,
    index = c("id", "time"), model = "within"
  )
  second <- plm::plm(
    y2 ~ x + z + factor(time), data = dat,
    index = c("id", "time"), model = "within"
  )
  fit <- suest(first, second, model_names = c("First", "Second"))
  cluster <- setNames(dat$id, rownames(dat))
  expected <- panel_fe_reference(list(first, second), cluster)

  expect_equal(fit$n_clusters, 120L)
  expect_equal(unname(vcov(fit)), unname(expected), tolerance = 1e-12)
  expect_equal(
    unname(coef(fit)["First::(Intercept)"]),
    mean(dat$y1) - sum(
      colMeans(stats::model.matrix(first, model = "pooling")[, names(coef(first))]) *
        coef(first)
    ),
    tolerance = 1e-12
  )
  expect_gt(max(abs(offdiag_vcov(fit))), 0)
})

test_that("plm within systems support unbalanced overlap and higher clusters", {
  skip_if_not_installed("plm")
  dat <- panel_fe_test_data()
  dat$y1[dat$time == 5 & dat$id %% 4 == 0] <- NA
  dat$y2[dat$time == 1 & dat$id %% 5 == 0] <- NA
  first <- plm::plm(
    y1 ~ x + z, data = dat, index = c("id", "time"), model = "within"
  )
  second <- plm::plm(
    y2 ~ x + z, data = dat, index = c("id", "time"), model = "within"
  )
  fit <- suest(
    first, second, model_names = c("First", "Second"), cluster = "higher"
  )
  higher <- setNames(dat$higher, rownames(dat))
  expected <- panel_fe_reference(list(first, second), higher)

  expect_equal(fit$n_clusters, 30L)
  expect_lt(fit$nobs_overlap, fit$nobs_union)
  expect_equal(unname(vcov(fit)), unname(expected), tolerance = 1e-12)
})

test_that("plm within systems work with marginaleffects", {
  skip_if_not_installed("plm")
  dat <- panel_fe_test_data()
  first <- plm::plm(
    y1 ~ x + z, data = dat, index = c("id", "time"), model = "within"
  )
  second <- plm::plm(
    y2 ~ x + z, data = dat, index = c("id", "time"), model = "within"
  )
  fit <- suest(first, second, model_names = c("First", "Second"))
  effects <- marginaleffects::avg_slopes(fit, variables = "x", newdata = dat)
  effects <- effects[match(c("First", "Second"), effects$group), ]

  expect_equal(
    effects$estimate,
    unname(c(coef(first)["x"], coef(second)["x"]))
  )
  expect_true(all(is.finite(effects$std.error)))

  predictions <- marginaleffects::avg_predictions(fit, newdata = dat)
  predictions <- predictions[match(c("First", "Second"), predictions$group), ]
  expect_equal(
    predictions$estimate,
    c(mean(dat$y1), mean(dat$y2)),
    tolerance = 1e-10
  )
})

test_that("plm within validation rejects unsupported panel specifications", {
  skip_if_not_installed("plm")
  dat <- panel_fe_test_data()
  within <- plm::plm(
    y1 ~ x + z, data = dat, index = c("id", "time"), model = "within"
  )
  random <- plm::plm(
    y1 ~ x + z, data = dat, index = c("id", "time"), model = "random"
  )
  time_fe <- plm::plm(
    y1 ~ x + z, data = dat, index = c("id", "time"),
    model = "within", effect = "time"
  )
  weighted <- plm::plm(
    y1 ~ x + z, data = dat, index = c("id", "time"),
    model = "within", weights = rep(1, nrow(dat))
  )
  linear <- lm(y1 ~ x + z, data = dat)

  expect_error(suest(within, random), "only with other panel models of the same type")
  expect_error(suest(within, time_fe), "Each model must be a supported")
  expect_error(suest(within, weighted), "Weighted plm")
  expect_error(suest(within, linear), "only with other panel models of the same type")

  dat$bad_cluster <- dat$higher
  dat$bad_cluster[dat$time == 5] <- dat$bad_cluster[dat$time == 5] + 100L
  expect_error(
    suest(within, within, cluster = "bad_cluster"),
    "Panel IDs must be nested within clusters"
  )
})

panel_be_reference <- function(models, cluster) {
  rows <- lapply(models, function(model) rownames(stats::model.frame(model)))
  union_rows <- unique(unlist(rows, use.names = FALSE))
  union_cluster <- as.character(cluster[match(union_rows, names(cluster))])

  transformed <- lapply(seq_along(models), function(i) {
    model <- models[[i]]
    X <- stats::model.matrix(model)
    residual <- as.numeric(stats::residuals(model))
    influence_panel <- (X * residual) %*% solve(crossprod(X))
    panel <- as.character(plm::index(model)[[1L]])
    first <- match(rownames(X), panel)
    influence <- matrix(0, nrow = length(panel), ncol = ncol(X))
    influence[first, ] <- influence_panel

    model_cluster <- cluster[match(rows[[i]], names(cluster))]
    panel_cluster <- model_cluster[first]
    g <- length(unique(panel_cluster))
    n_panel <- nrow(X)
    k <- ncol(X)
    correction <- (g / (g - 1)) * ((n_panel - 1) / (n_panel - k))
    aligned <- matrix(0, nrow = length(union_rows), ncol = ncol(influence))
    aligned[match(rows[[i]], union_rows), ] <- influence * sqrt(correction)
    aligned
  })

  influence <- do.call(cbind, transformed)
  clustered <- rowsum(influence, union_cluster, reorder = FALSE)
  V <- crossprod(clustered)
  (V + t(V)) / 2
}

test_that("plm between systems match the xtreg BE sandwich construction", {
  skip_if_not_installed("plm")
  dat <- panel_fe_test_data()
  first <- plm::plm(
    y1 ~ x + z, data = dat, index = c("id", "time"), model = "between"
  )
  second <- plm::plm(
    y2 ~ x + z, data = dat, index = c("id", "time"), model = "between"
  )
  fit <- suest(first, second, model_names = c("First", "Second"))
  cluster <- setNames(dat$id, rownames(dat))
  expected <- panel_be_reference(list(first, second), cluster)

  expect_equal(fit$n_clusters, 120L)
  expect_equal(unname(vcov(fit)), unname(expected), tolerance = 1e-12)
  expect_equal(
    coef(fit),
    c(
      setNames(coef(first), paste0("First::", names(coef(first)))),
      setNames(coef(second), paste0("Second::", names(coef(second))))
    )
  )
})

test_that("plm between supports unbalanced samples and higher clusters", {
  skip_if_not_installed("plm")
  dat <- panel_fe_test_data()
  dat$y1[dat$time == 5 & dat$id %% 4 == 0] <- NA
  dat$y2[dat$time == 1 & dat$id %% 5 == 0] <- NA
  first <- plm::plm(
    y1 ~ x + z, data = dat, index = c("id", "time"), model = "between"
  )
  second <- plm::plm(
    y2 ~ x + z, data = dat, index = c("id", "time"), model = "between"
  )
  fit <- suest(
    first, second, model_names = c("First", "Second"), cluster = "higher"
  )
  higher <- setNames(dat$higher, rownames(dat))
  expected <- panel_be_reference(list(first, second), higher)

  expect_equal(fit$n_clusters, 30L)
  expect_lt(fit$nobs_overlap, fit$nobs_union)
  expect_equal(unname(vcov(fit)), unname(expected), tolerance = 1e-12)

  effects <- marginaleffects::avg_slopes(fit, variables = "x", newdata = dat)
  effects <- effects[match(c("First", "Second"), effects$group), ]
  expect_equal(
    effects$estimate,
    unname(c(coef(first)["x"], coef(second)["x"]))
  )
  expect_true(all(is.finite(effects$std.error)))
})

panel_re_reference <- function(models, cluster) {
  rows <- lapply(models, function(model) rownames(stats::model.frame(model)))
  union_rows <- unique(unlist(rows, use.names = FALSE))
  union_cluster <- as.character(cluster[match(union_rows, names(cluster))])

  transformed <- lapply(seq_along(models), function(i) {
    model <- models[[i]]
    X <- stats::model.matrix(model)
    influence <- (X * as.numeric(stats::residuals(model))) %*%
      solve(crossprod(X))
    model_cluster <- cluster[match(rows[[i]], names(cluster))]
    g <- length(unique(model_cluster))
    n <- nrow(X)
    k <- ncol(X)
    correction <- (g / (g - 1)) * ((n - 1) / (n - k))
    aligned <- matrix(0, nrow = length(union_rows), ncol = ncol(influence))
    aligned[match(rows[[i]], union_rows), ] <- influence * sqrt(correction)
    aligned
  })

  influence <- do.call(cbind, transformed)
  clustered <- rowsum(influence, union_cluster, reorder = FALSE)
  V <- crossprod(clustered)
  (V + t(V)) / 2
}

test_that("plm random systems match the xtreg RE sandwich construction", {
  skip_if_not_installed("plm")
  dat <- panel_fe_test_data()
  first <- plm::plm(
    y1 ~ x + z, data = dat, index = c("id", "time"),
    model = "random", random.method = "swar"
  )
  second <- plm::plm(
    y2 ~ x + z, data = dat, index = c("id", "time"),
    model = "random", random.method = "swar"
  )
  fit <- suest(first, second, model_names = c("First", "Second"))
  cluster <- setNames(dat$id, rownames(dat))
  expected <- panel_re_reference(list(first, second), cluster)

  expect_equal(fit$n_clusters, 120L)
  expect_equal(unname(vcov(fit)), unname(expected), tolerance = 1e-12)
  expect_equal(
    coef(fit),
    c(
      setNames(coef(first), paste0("First::", names(coef(first)))),
      setNames(coef(second), paste0("Second::", names(coef(second))))
    )
  )
})

test_that("plm random supports unbalanced samples and higher clusters", {
  skip_if_not_installed("plm")
  dat <- panel_fe_test_data()
  dat$y1[dat$time == 5 & dat$id %% 4 == 0] <- NA
  dat$y2[dat$time == 1 & dat$id %% 5 == 0] <- NA
  first <- plm::plm(
    y1 ~ x + z, data = dat, index = c("id", "time"), model = "random"
  )
  second <- plm::plm(
    y2 ~ x + z, data = dat, index = c("id", "time"), model = "random"
  )
  fit <- suest(
    first, second, model_names = c("First", "Second"), cluster = "higher"
  )
  higher <- setNames(dat$higher, rownames(dat))
  expected <- panel_re_reference(list(first, second), higher)

  expect_equal(fit$n_clusters, 30L)
  expect_lt(fit$nobs_overlap, fit$nobs_union)
  expect_equal(unname(vcov(fit)), unname(expected), tolerance = 1e-12)

  effects <- marginaleffects::avg_slopes(fit, variables = "x", newdata = dat)
  effects <- effects[match(c("First", "Second"), effects$group), ]
  expect_equal(
    effects$estimate,
    unname(c(coef(first)["x"], coef(second)["x"]))
  )
  expect_true(all(is.finite(effects$std.error)))
})

test_that("plm random requires the Stata-compatible Swamy-Arora method", {
  skip_if_not_installed("plm")
  dat <- panel_fe_test_data()
  swar <- plm::plm(
    y1 ~ x + z, data = dat, index = c("id", "time"), model = "random"
  )
  walhus <- plm::plm(
    y2 ~ x + z, data = dat, index = c("id", "time"),
    model = "random", random.method = "walhus"
  )

  expect_error(suest(swar, walhus), "Each model must be a supported")
})
