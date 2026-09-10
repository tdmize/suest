panel_logit_re_test_data <- function(groups = 70L, periods = 6L) {
  set.seed(10427)
  id <- rep(seq_len(groups), each = periods)
  time <- rep(seq_len(periods), groups)
  higher <- ceiling(id / 5)
  x <- rnorm(length(id))
  z <- rnorm(length(id))
  woman <- factor(rbinom(length(id), 1, 0.48), labels = c("Men", "Women"))
  common <- rnorm(groups, sd = 1.15)[id]
  second <- rnorm(groups, sd = 0.55)[id]
  eta1 <- -0.35 + 0.65*x - 0.30*z + 0.20*(woman == "Women") + common
  eta2 <- 0.15 + 0.40*x + 0.25*z - 0.15*(woman == "Women") +
    0.75*common + second
  data.frame(
    id, time, higher, x, z, woman,
    y1 = rbinom(length(id), 1, plogis(eta1)),
    y2 = rbinom(length(id), 1, plogis(eta2))
  )
}

panel_logit_re_loglik <- function(parameters, model, panels = NULL) {
  fit_environment <- environment(model$objectiveFn)
  X <- get("X", fit_environment, inherits = FALSE)
  y <- as.numeric(get("y", fit_environment, inherits = FALSE))
  id <- as.character(get("id", fit_environment, inherits = FALSE))
  quadrature <- get("rn", fit_environment, inherits = FALSE)
  beta_names <- setdiff(names(parameters), "sigma")
  eta <- as.numeric(X[, beta_names, drop = FALSE] %*% parameters[beta_names])
  sigma <- unname(parameters["sigma"])
  q <- 2*y - 1
  if (is.null(panels)) panels <- unique(id)

  sum(vapply(panels, function(panel) {
    rows <- id == panel
    likelihood <- sum(vapply(seq_along(quadrature$nodes), function(j) {
      quadrature$weights[j] * prod(plogis(
        q[rows] * (eta[rows] + sqrt(2)*sigma*quadrature$nodes[j])
      ))
    }, numeric(1))) / sqrt(pi)
    log(likelihood)
  }, numeric(1)))
}

panel_logit_re_start <- function(formula, data, link = "logit", random = TRUE) {
  start <- coef(glm(formula, data = data, family = binomial(link)))
  if (random) c(start, sigma = 1) else start
}

panel_logit_re_stata_model <- function(model) {
  parameters <- coef(model)
  sigma_sign <- sign(parameters["sigma"])
  parameters["sigma"] <- abs(parameters["sigma"])
  covariance <- vcov(model)
  canonical <- diag(length(parameters))
  canonical[nrow(canonical), ncol(canonical)] <- sigma_sign
  covariance <- t(canonical) %*% covariance %*% canonical
  jacobian <- diag(length(parameters))
  jacobian[nrow(jacobian), ncol(jacobian)] <- 2/parameters["sigma"]
  parameters["sigma"] <- 2*log(parameters["sigma"])
  names(parameters)[length(parameters)] <- "lnsig2u"
  covariance <- jacobian %*% covariance %*% t(jacobian)
  dimnames(covariance) <- list(names(parameters), names(parameters))
  list(coefficients = parameters, vcov = covariance,
    logLik = as.numeric(logLik(model)))
}

panel_logit_re_stata_system <- function(fit) {
  parameters <- coef(fit)
  sigma <- grep("::sigma$", names(parameters))
  jacobian <- diag(length(parameters))
  jacobian[cbind(sigma, sigma)] <- 2/parameters[sigma]
  parameters[sigma] <- 2*log(parameters[sigma])
  names(parameters)[sigma] <- sub("sigma$", "lnsig2u", names(parameters)[sigma])
  covariance <- jacobian %*% vcov(fit) %*% t(jacobian)
  dimnames(covariance) <- list(names(parameters), names(parameters))
  list(coefficients = parameters, vcov = covariance,
    nobs_union = fit$nobs_union, nobs_overlap = fit$nobs_overlap,
    n_clusters = fit$n_clusters)
}

test_that("pglm random-effects logit scores and bread match its likelihood", {
  skip_if_not_installed("pglm")
  maxLik <- maxLik::maxLik
  dat <- panel_logit_re_test_data()
  model <- pglm::pglm(
    y1 ~ x + z, data = dat, family = binomial("logit"),
    model = "random", effect = "individual", index = c("id", "time"),
    R = 12, method = "bfgs", print.level = 0, reltol = 1e-12,
    start = panel_logit_re_start(y1 ~ x + z, dat)
  )
  component <- suest:::.suest_pglm_logit_re_components(model)
  parameters <- component$parameters
  panel <- as.character(suest:::.suest_panel_id(model))

  expect_equal(
    panel_logit_re_loglik(parameters, model),
    as.numeric(logLik(model)),
    tolerance = 1e-9
  )
  expect_lt(max(abs(colSums(component$score))), 2e-5)
  expect_equal(
    unname(component$bread / nrow(component$score)),
    unname({
      sign_matrix <- diag(length(coef(model)))
      sign_matrix[nrow(sign_matrix), ncol(sign_matrix)] <- sign(coef(model)["sigma"])
      t(sign_matrix) %*% vcov(model) %*% sign_matrix
    }),
    tolerance = 1e-10
  )

  first_panel <- unique(panel)[1L]
  rows <- which(panel == first_panel)
  step <- 1e-5*pmax(1, abs(parameters))
  numerical_score <- vapply(seq_along(parameters), function(j) {
    plus <- minus <- parameters
    plus[j] <- plus[j] + step[j]
    minus[j] <- minus[j] - step[j]
    (panel_logit_re_loglik(plus, model, first_panel) -
       panel_logit_re_loglik(minus, model, first_panel))/(2*step[j])
  }, numeric(1))
  names(numerical_score) <- names(parameters)
  expect_equal(
    colSums(component$score[rows, , drop = FALSE]),
    numerical_score,
    tolerance = 2e-6
  )
  expect_equal(sum(rowSums(abs(component$score)) > 0), length(unique(panel)))
})

test_that("pglm random-effects logit systems include sigma and panel clustering", {
  skip_if_not_installed("pglm")
  maxLik <- maxLik::maxLik
  dat <- panel_logit_re_test_data()
  first <- pglm::pglm(
    y1 ~ x + z, data = dat, family = binomial("logit"),
    model = "random", index = c("id", "time"), R = 12,
    method = "bfgs", print.level = 0, reltol = 1e-12,
    start = panel_logit_re_start(y1 ~ x + z, dat)
  )
  second <- pglm::pglm(
    y2 ~ x + z, data = dat, family = binomial("logit"),
    model = "random", index = c("id", "time"), R = 12,
    method = "bfgs", print.level = 0, reltol = 1e-12,
    start = panel_logit_re_start(y2 ~ x + z, dat)
  )
  fit <- suest(
    first, second, model_names = c("First", "Second"),
    observation_id = c("id", "time")
  )

  expect_equal(fit$n_clusters, length(unique(dat$id)))
  expect_true(all(c("First::sigma", "Second::sigma") %in% names(coef(fit))))
  expect_gt(unname(coef(fit)["First::sigma"]), 0)
  expect_gt(unname(coef(fit)["Second::sigma"]), 0)
  expect_true(all(is.finite(vcov(fit))))
  expect_equal(vcov(fit), t(vcov(fit)), tolerance = 1e-12)
  expect_gt(max(abs(offdiag_vcov(fit))), 0)
})

test_that("pglm random-effects logit supports partial and disjoint panel systems", {
  skip_if_not_installed("pglm")
  maxLik <- maxLik::maxLik
  dat <- panel_logit_re_test_data(groups = 80L)
  dat$y1[dat$time == 6 & dat$id %% 4 == 0] <- NA
  dat$y2[dat$time == 1 & dat$id %% 5 == 0] <- NA
  first <- pglm::pglm(
    y1 ~ x + z, data = dat, family = binomial("logit"),
    model = "random", index = c("id", "time"), R = 12,
    method = "bfgs", print.level = 0, reltol = 1e-11,
    start = panel_logit_re_start(y1 ~ x + z, dat)
  )
  second <- pglm::pglm(
    y2 ~ x + z, data = dat, family = binomial("logit"),
    model = "random", index = c("id", "time"), R = 12,
    method = "bfgs", print.level = 0, reltol = 1e-11,
    start = panel_logit_re_start(y2 ~ x + z, dat)
  )
  partial <- suest(
    first, second, model_names = c("First", "Second"),
    observation_id = c("id", "time"), cluster = "higher"
  )
  expect_equal(partial$n_clusters, length(unique(dat$higher)))
  expect_lt(partial$nobs_overlap, partial$nobs_union)
  expect_true(all(is.finite(vcov(partial))))

  low <- dat[dat$id <= 40, ]
  high <- dat[dat$id > 40, ]
  low_fit <- pglm::pglm(
    y1 ~ x + z, data = low, family = binomial("logit"),
    model = "random", index = c("id", "time"), R = 12,
    method = "bfgs", print.level = 0, reltol = 1e-11,
    start = panel_logit_re_start(y1 ~ x + z, low)
  )
  high_fit <- pglm::pglm(
    y2 ~ x + z, data = high, family = binomial("logit"),
    model = "random", index = c("id", "time"), R = 12,
    method = "bfgs", print.level = 0, reltol = 1e-11,
    start = panel_logit_re_start(y2 ~ x + z, high)
  )
  disjoint <- suest(
    low_fit, high_fit, model_names = c("Low", "High"),
    observation_id = c("id", "time")
  )
  expect_equal(disjoint$nobs_overlap, 0L)
  expect_equal(
    unname(offdiag_vcov(disjoint)), matrix(0, 4L, 4L), tolerance = 1e-12
  )
})

test_that("pglm integrated predictions and marginal effects use sigma", {
  skip_if_not_installed("pglm")
  maxLik <- maxLik::maxLik
  dat <- panel_logit_re_test_data()
  first <- pglm::pglm(
    y1 ~ x*woman + z, data = dat, family = binomial("logit"),
    model = "random", index = c("id", "time"), R = 12,
    method = "bfgs", print.level = 0, reltol = 1e-12,
    start = panel_logit_re_start(y1 ~ x*woman + z, dat)
  )
  second <- pglm::pglm(
    y2 ~ x*woman + z, data = dat, family = binomial("logit"),
    model = "random", index = c("id", "time"), R = 12,
    method = "bfgs", print.level = 0, reltol = 1e-12,
    start = panel_logit_re_start(y2 ~ x*woman + z, dat)
  )
  fit <- suest(
    first, second, model_names = c("First", "Second"),
    observation_id = c("id", "time")
  )
  predictions <- marginaleffects::avg_predictions(fit, newdata = dat)
  slopes <- marginaleffects::avg_slopes(fit, variables = "x", newdata = dat)

  expected <- vapply(list(first, second), function(model) {
    mean(suest:::.suest_predict_values(
      model, dat, "response", "pglm::pglm", "panel_logit_re"
    ))
  }, numeric(1))
  expect_equal(predictions$estimate, expected, tolerance = 1e-8)
  expect_true(all(is.finite(predictions$std.error)))
  expect_true(all(is.finite(slopes$estimate)))
  expect_true(all(is.finite(slopes$std.error)))

  changed <- coef(fit)
  changed["First::sigma"] <- changed["First::sigma"]*1.5
  altered <- marginaleffects::set_coef(fit, changed)
  original_prediction <- marginaleffects::predictions(
    fit, newdata = dat[1:5, ], vcov = FALSE
  )
  altered_prediction <- marginaleffects::predictions(
    altered, newdata = dat[1:5, ], vcov = FALSE
  )
  first_rows <- original_prediction$group == "First"
  expect_gt(max(abs(
    original_prediction$estimate[first_rows] -
      altered_prediction$estimate[first_rows]
  )), 1e-5)
})

test_that("pglm random-effects logit validation is deliberately narrow", {
  skip_if_not_installed("pglm")
  maxLik <- maxLik::maxLik
  dat <- panel_logit_re_test_data(groups = 45L)
  wrong_points <- pglm::pglm(
    y1 ~ x + z, data = dat, family = binomial("logit"),
    model = "random", index = c("id", "time"), R = 8,
    method = "bfgs", print.level = 0,
    start = panel_logit_re_start(y1 ~ x + z, dat)
  )
  probit <- pglm::pglm(
    y1 ~ x + z, data = dat, family = binomial("probit"),
    model = "random", index = c("id", "time"), R = 12,
    method = "bfgs", print.level = 0,
    start = panel_logit_re_start(y1 ~ x + z, dat, link = "probit")
  )
  pooling <- pglm::pglm(
    y1 ~ x + z, data = dat, family = binomial("logit"),
    model = "pooling", index = c("id", "time"),
    method = "bfgs", print.level = 0,
    start = panel_logit_re_start(y1 ~ x + z, dat, random = FALSE)
  )

  expect_error(suest(wrong_points, wrong_points), "R = 12", fixed = TRUE)
  expect_s3_class(suest(probit, probit), "suest_model")
  expect_error(suest(pooling, pooling), "individual random-effects binary")
})

test_that("pglm panel logit reproduces the returned Stata benchmark", {
  skip_if_not_installed("pglm")
  skip_if_not_installed("marginaleffects")
  maxLik <- maxLik::maxLik
  fixture <- readRDS(test_path("fixtures", "panel-logit-re-stata.rds"))
  dat <- fixture$data
  fit_model <- function(formula) pglm::pglm(
    formula, data = dat, family = binomial("logit"), model = "random",
    effect = "individual", index = c("id", "time"), R = 12,
    method = "bfgs", print.level = 0, reltol = 1e-12,
    start = panel_logit_re_start(formula, dat)
  )

  y1 <- fit_model(y1 ~ x + z)
  y2 <- fit_model(y2 ~ x + z)
  for (i in seq_along(list(y1, y2))) {
    actual <- panel_logit_re_stata_model(list(y1, y2)[[i]])
    expected <- fixture$balanced[[paste0("model", i)]]
    expect_lt(max(abs(actual$coefficients - expected$coefficients)), 2e-6)
    expect_lt(max(abs(actual$vcov - expected$vcov)), 2e-7)
    expect_lt(abs(actual$logLik - expected$logLik), 2e-6)
  }

  balanced <- suest(y1, y2, model_names = c("Y1", "Y2"),
    observation_id = c("id", "time"))
  predictions <- marginaleffects::avg_predictions(balanced, newdata = dat)
  slopes <- marginaleffects::avg_slopes(
    balanced, variables = "x", newdata = dat)
  expect_lt(max(abs(predictions$estimate -
    fixture$balanced$prediction_mean)), 2e-6)
  expect_lt(max(abs(slopes$estimate -
    fixture$balanced$average_slope_x)), 2e-6)

  partial <- suest(
    fit_model(y1_partial ~ x + z), fit_model(y2_partial ~ x + z),
    model_names = c("Y1", "Y2"), observation_id = c("id", "time"),
    cluster = "higher"
  )
  disjoint <- suest(
    fit_model(y1_left ~ x + z), fit_model(y2_right ~ x + z),
    model_names = c("Y1", "Y2"), observation_id = c("id", "time")
  )

  cases <- list(
    balanced = list(fit = balanced, coefficient = 2e-4, covariance = 3e-4),
    partial_higher = list(
      fit = partial, coefficient = 2e-4, covariance = 4e-4),
    disjoint = list(fit = disjoint, coefficient = 2e-3, covariance = 4e-3)
  )
  for (case in names(cases)) {
    actual <- panel_logit_re_stata_system(cases[[case]]$fit)
    expected <- fixture[[case]]$system
    expect_lt(max(abs(actual$coefficients - expected$coefficients)),
      cases[[case]]$coefficient)
    expect_lt(max(abs(actual$vcov - expected$vcov)),
      cases[[case]]$covariance)
    expect_equal(actual[c("nobs_union", "nobs_overlap", "n_clusters")],
      expected[c("nobs_union", "nobs_overlap", "n_clusters")])
  }
  expect_equal(disjoint$nobs_overlap, 0L)
  expect_equal(unname(vcov(disjoint)[1:4, 5:8]), matrix(0, 4L, 4L),
    tolerance = 1e-12)
})
