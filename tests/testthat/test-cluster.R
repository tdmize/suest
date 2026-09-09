cluster_system_reference <- function(models, cluster, types = NULL) {
  if (is.null(types))
    types <- vapply(models, function(model) {
      suest:::.suest_model_adapter(model)$type
    }, character(1))
  engines <- lapply(models, function(model) {
    suest:::.suest_model_adapter(model)$engine
  })
  components <- Map(
    function(model, type, engine) {
      suest:::.suest_model_components(model, type, engine)
    },
    models,
    types,
    engines
  )
  frames <- Map(
    function(model, engine) suest:::.suest_model_frame(model, engine),
    models,
    engines
  )
  rows <- lapply(frames, rownames)
  union_rows <- unique(unlist(rows, use.names = FALSE))
  n_union <- length(union_rows)
  union_cluster <- as.character(cluster[match(union_rows, names(cluster))])
  g <- length(unique(union_cluster))

  aligned <- lapply(seq_along(models), function(i) {
    U <- components[[i]]$score
    correction <- if (types[i] == "ivreg") 1 else g / (g - 1)
    U <- U * sqrt(correction)
    out <- matrix(0, nrow = n_union, ncol = ncol(U))
    out[match(rows[[i]], union_rows), ] <- U
    out
  })
  breads <- lapply(seq_along(models), function(i) {
    components[[i]]$bread * n_union / nrow(components[[i]]$score)
  })
  U <- do.call(cbind, aligned)
  B <- suest:::.suest_block_diag(breads)
  C <- rowsum(U, union_cluster, reorder = FALSE)
  V <- B %*% crossprod(C) %*% B / n_union^2
  (V + t(V)) / 2
}

cluster_test_data <- function(n = 600L) {
  set.seed(6109)
  cluster <- rep(seq_len(n / 4L), each = 4L)
  cluster_effect <- rnorm(max(cluster))[cluster]
  x <- rnorm(n) + 0.25*cluster_effect
  z <- rnorm(n)
  y <- 1 + 0.45*x - 0.30*z + 0.60*cluster_effect + rnorm(n)
  b <- rbinom(
    n,
    1,
    plogis(-0.25 + 0.50*x - 0.20*z + 0.45*cluster_effect)
  )
  data.frame(id = seq_len(n), cluster = cluster, x = x, z = z, y = y, b = b)
}

test_that("clustered systems match direct grouped-score calculations", {
  dat <- cluster_test_data()
  linear <- lm(y ~ x + z, data = dat)
  logit <- glm(b ~ x + z, family = binomial(), data = dat)
  fit <- suest(
    linear,
    logit,
    model_names = c("Linear", "Logit"),
    cluster = "cluster"
  )
  cluster <- setNames(dat$cluster, rownames(dat))
  expected <- cluster_system_reference(list(linear, logit), cluster)

  expect_equal(fit$n_clusters, 150L)
  expect_equal(unname(vcov(fit)), unname(expected), tolerance = 1e-12)
  expect_true("Linear::lnvar" %in% names(coef(fit)))

  effects <- marginaleffects::avg_slopes(
    fit,
    variables = "x",
    newdata = dat
  )
  expect_equal(nrow(effects), 2L)
  expect_true(all(is.finite(effects$std.error)))
})

test_that("clusters align across partially overlapping model samples", {
  dat <- cluster_test_data()
  linear <- lm(y ~ x + z, data = dat, subset = id <= 450)
  logit <- glm(
    b ~ x + z,
    family = binomial(),
    data = dat,
    subset = id >= 151
  )
  fit <- suest(linear, logit, cluster = "cluster")
  cluster <- setNames(dat$cluster, rownames(dat))
  expected <- cluster_system_reference(list(linear, logit), cluster)

  expect_equal(fit$nobs_overlap, 300L)
  expect_equal(fit$n_clusters, 150L)
  expect_equal(unname(vcov(fit)), unname(expected), tolerance = 1e-12)
})

test_that("disjoint observations can contribute to shared clusters", {
  dat <- cluster_test_data()
  dat$paired_cluster <- rep(seq_len(nrow(dat) / 2L), each = 2L)
  dat$side <- rep(c("left", "right"), nrow(dat) / 2L)
  dat$unique_cluster <- interaction(
    dat$side,
    dat$paired_cluster,
    drop = TRUE
  )
  linear <- lm(y ~ x + z, data = dat, subset = side == "left")
  logit <- glm(
    b ~ x + z,
    family = binomial(),
    data = dat,
    subset = side == "right"
  )
  shared <- suest(linear, logit, cluster = "paired_cluster")
  separate <- suest(linear, logit, cluster = "unique_cluster")

  expect_equal(shared$nobs_overlap, 0L)
  expect_gt(max(abs(offdiag_vcov(shared))), 0)
  expect_true(all(offdiag_vcov(separate) == 0))
})

test_that("cluster IDs can be supplied directly and are validated", {
  dat <- cluster_test_data()
  linear <- lm(y ~ x, data = dat)
  logit <- glm(b ~ x, family = binomial(), data = dat)
  from_column <- suest(
    linear,
    logit,
    model_names = c("Linear", "Logit"),
    cluster = "cluster"
  )
  from_list <- suest(
    linear,
    logit,
    model_names = c("Linear", "Logit"),
    cluster = list(Linear = dat$cluster, Logit = dat$cluster)
  )
  expect_equal(vcov(from_list), vcov(from_column), tolerance = 1e-12)

  bad <- dat$cluster
  bad[1L] <- max(bad) + 1L
  expect_error(
    suest(
      linear,
      logit,
      model_names = c("Linear", "Logit"),
      cluster = list(Linear = dat$cluster, Logit = bad)
    ),
    "Cluster IDs disagree"
  )
  expect_error(
    suest(linear, logit, cluster = list(dat$cluster)),
    "one cluster vector"
  )
  expect_error(
    suest(linear, logit, cluster = list(rep(NA, nrow(dat)), dat$cluster)),
    "cannot be missing"
  )
  expect_error(
    suest(linear, logit, cluster = list(rep(1, nrow(dat)), rep(1, nrow(dat)))),
    "at least two clusters"
  )
})

test_that("clustered fixest 2SLS uses raw grouped influences", {
  skip_if_not_installed("fixest")
  set.seed(927)
  n <- 600L
  dat <- data.frame(
    id = seq_len(n),
    cluster = rep(seq_len(n / 4L), each = 4L),
    x = rnorm(n),
    z1 = rnorm(n),
    z2 = rnorm(n)
  )
  disturbance <- rnorm(n)
  dat$endogenous <- 0.8*dat$z1 + 0.2*dat$z2 +
    0.4*disturbance + rnorm(n)
  dat$y <- 0.5 + 0.6*dat$x + 1.1*dat$endogenous + disturbance
  base <- fixest::feols(y ~ x | endogenous ~ z1, data = dat)
  adjusted <- fixest::feols(y ~ x | endogenous ~ z1 + z2, data = dat)
  fit <- suest(base, adjusted, cluster = "cluster")
  cluster <- setNames(dat$cluster, rownames(dat))
  expected <- cluster_system_reference(
    list(base, adjusted),
    cluster,
    types = c("ivreg", "ivreg")
  )

  expect_equal(unname(vcov(fit)), unname(expected), tolerance = 1e-12)
})
