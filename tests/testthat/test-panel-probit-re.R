panel_probit_re_test_data <- function(groups = 70L, periods = 6L) {
  set.seed(29107)
  id <- rep(seq_len(groups), each = periods)
  time <- rep(seq_len(periods), groups)
  higher <- ceiling(id / 5)
  x <- rnorm(length(id))
  z <- rnorm(length(id))
  woman <- factor(rbinom(length(id), 1, 0.48), labels = c("Men", "Women"))
  common <- rnorm(groups, sd = 0.85)[id]
  second <- rnorm(groups, sd = 0.45)[id]
  eta1 <- -0.25 + 0.55*x - 0.30*z + 0.20*(woman == "Women") + common
  eta2 <- 0.15 + 0.40*x + 0.25*z - 0.15*(woman == "Women") +
    0.70*common + second
  data.frame(
    id, time, higher, x, z, woman,
    y1 = rbinom(length(id), 1, pnorm(eta1)),
    y2 = rbinom(length(id), 1, pnorm(eta2))
  )
}

panel_probit_re_start <- function(formula, data) {
  c(coef(glm(formula, data = data, family = binomial("probit"))), sigma = 1)
}

panel_probit_re_fit <- function(formula, data) {
  maxLik <- maxLik::maxLik
  pglm::pglm(
    formula, data = data, family = binomial("probit"), model = "random",
    effect = "individual", index = c("id", "time"), R = 12,
    method = "bfgs", print.level = 0, reltol = 1e-12,
    start = panel_probit_re_start(formula, data)
  )
}

panel_probit_re_loglik <- function(parameters, model, panels = NULL) {
  fit_environment <- environment(model$objectiveFn)
  X <- get("X", fit_environment, inherits = FALSE)
  y <- as.numeric(get("y", fit_environment, inherits = FALSE))
  id <- as.character(get("id", fit_environment, inherits = FALSE))
  quadrature <- get("rn", fit_environment, inherits = FALSE)
  beta_names <- setdiff(names(parameters), "sigma")
  eta <- as.numeric(X[, beta_names, drop = FALSE] %*% parameters[beta_names])
  q <- 2*y - 1
  if (is.null(panels)) panels <- unique(id)

  sum(vapply(panels, function(panel) {
    rows <- id == panel
    likelihood <- sum(vapply(seq_along(quadrature$nodes), function(j) {
      quadrature$weights[j] * prod(pnorm(
        q[rows] * (eta[rows] + sqrt(2)*parameters["sigma"]*quadrature$nodes[j])
      ))
    }, numeric(1))) / sqrt(pi)
    log(likelihood)
  }, numeric(1)))
}

panel_probit_re_stata_model <- function(model) {
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

panel_probit_re_stata_system <- function(fit) {
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

test_that("pglm random-effects probit scores and bread match its likelihood", {
  skip_if_not_installed("pglm")
  maxLik <- maxLik::maxLik
  dat <- panel_probit_re_test_data()
  model <- panel_probit_re_fit(y1 ~ x + z, dat)
  adapter <- suest:::.suest_model_adapter(model)
  component <- suest:::.suest_pglm_binary_re_components(
    model, "panel_probit_re")
  parameters <- component$parameters
  panel <- as.character(suest:::.suest_panel_id(model))

  expect_identical(adapter$type, "panel_probit_re")
  expect_equal(panel_probit_re_loglik(parameters, model),
    as.numeric(logLik(model)), tolerance = 1e-9)
  expect_lt(max(abs(colSums(component$score))), 3e-5)
  expect_equal(component$bread/nrow(component$score), {
    sign_matrix <- diag(length(coef(model)))
    sign_matrix[nrow(sign_matrix), ncol(sign_matrix)] <-
      sign(coef(model)["sigma"])
    t(sign_matrix) %*% vcov(model) %*% sign_matrix
  }, tolerance = 1e-10, ignore_attr = TRUE)

  first_panel <- unique(panel)[1L]
  rows <- which(panel == first_panel)
  step <- 1e-5*pmax(1, abs(parameters))
  numerical_score <- vapply(seq_along(parameters), function(j) {
    plus <- minus <- parameters
    plus[j] <- plus[j] + step[j]
    minus[j] <- minus[j] - step[j]
    (panel_probit_re_loglik(plus, model, first_panel) -
       panel_probit_re_loglik(minus, model, first_panel))/(2*step[j])
  }, numeric(1))
  expect_equal(colSums(component$score[rows, , drop = FALSE]),
    numerical_score, tolerance = 2e-6, ignore_attr = TRUE)
  expect_equal(sum(rowSums(abs(component$score)) > 0), length(unique(panel)))
})

test_that("pglm random-effects probit supports panel system alignment", {
  skip_if_not_installed("pglm")
  maxLik <- maxLik::maxLik
  dat <- panel_probit_re_test_data(groups = 80L)
  dat$y1[dat$time == 6 & dat$id %% 4 == 0] <- NA
  dat$y2[dat$time == 1 & dat$id %% 5 == 0] <- NA
  first <- panel_probit_re_fit(y1 ~ x + z, dat)
  second <- panel_probit_re_fit(y2 ~ x + z, dat)
  higher_cluster <- lapply(list(first, second), function(model) {
    frame <- suest:::.suest_model_frame(model, "pglm::pglm")
    dat$higher[match(rownames(frame), rownames(dat))]
  })
  partial <- suest(first, second, model_names = c("First", "Second"),
    observation_id = c("id", "time"), cluster = higher_cluster)

  expect_equal(partial$n_clusters, length(unique(dat$higher)))
  expect_lt(partial$nobs_overlap, partial$nobs_union)
  expect_true(all(c("First::sigma", "Second::sigma") %in% names(coef(partial))))
  expect_true(all(coef(partial)[c("First::sigma", "Second::sigma")] > 0))
  expect_true(all(is.finite(vcov(partial))))
  expect_gt(max(abs(offdiag_vcov(partial))), 0)

  low <- dat[dat$id <= 40, ]
  high <- dat[dat$id > 40, ]
  disjoint <- suest(
    panel_probit_re_fit(y1 ~ x + z, low),
    panel_probit_re_fit(y2 ~ x + z, high),
    model_names = c("Low", "High"), observation_id = c("id", "time")
  )
  expect_equal(disjoint$nobs_overlap, 0L)
  expect_equal(unname(offdiag_vcov(disjoint)), matrix(0, 4L, 4L),
    tolerance = 1e-12)
})

test_that("pglm integrated probit predictions work with marginaleffects", {
  skip_if_not_installed("pglm")
  skip_if_not_installed("marginaleffects")
  maxLik <- maxLik::maxLik
  dat <- panel_probit_re_test_data()
  first <- panel_probit_re_fit(y1 ~ x*woman + z, dat)
  second <- panel_probit_re_fit(y2 ~ x*woman + z, dat)
  fit <- suest(first, second, model_names = c("First", "Second"),
    observation_id = c("id", "time"))
  predictions <- marginaleffects::avg_predictions(fit, newdata = dat)
  slopes <- marginaleffects::avg_slopes(fit, variables = "x", newdata = dat)

  expected <- vapply(list(first, second), function(model) mean(
    suest:::.suest_predict_values(
      model, dat, "response", "pglm::pglm", "panel_probit_re"
    )
  ), numeric(1))
  expect_equal(predictions$estimate, expected, tolerance = 1e-8)
  expect_true(all(is.finite(predictions$std.error)))
  expect_true(all(is.finite(slopes$estimate)))
  expect_true(all(is.finite(slopes$std.error)))

  changed <- coef(fit)
  changed["First::sigma"] <- changed["First::sigma"]*1.5
  altered <- marginaleffects::set_coef(fit, changed)
  original_prediction <- marginaleffects::predictions(
    fit, newdata = dat[1:5, ], vcov = FALSE)
  altered_prediction <- marginaleffects::predictions(
    altered, newdata = dat[1:5, ], vcov = FALSE)
  first_rows <- original_prediction$group == "First"
  expect_gt(max(abs(original_prediction$estimate[first_rows] -
    altered_prediction$estimate[first_rows])), 1e-5)
})

test_that("random-effects panel logit and probit cannot be mixed yet", {
  skip_if_not_installed("pglm")
  maxLik <- maxLik::maxLik
  dat <- panel_probit_re_test_data(groups = 45L)
  probit <- panel_probit_re_fit(y1 ~ x + z, dat)
  logit <- pglm::pglm(
    y2 ~ x + z, data = dat, family = binomial("logit"), model = "random",
    effect = "individual", index = c("id", "time"), R = 12,
    method = "bfgs", print.level = 0,
    start = c(coef(glm(y2 ~ x + z, dat, family = binomial("logit"))), sigma = 1)
  )
  expect_error(suest(probit, logit), "same type", fixed = TRUE)
})

test_that("pglm panel probit reproduces the returned Stata benchmark", {
  skip_if_not_installed("pglm")
  skip_if_not_installed("marginaleffects")
  maxLik <- maxLik::maxLik
  fixture <- readRDS(test_path("fixtures", "panel-probit-re-stata.rds"))
  dat <- fixture$data
  fit_model <- function(formula) pglm::pglm(
    formula, data = dat, family = binomial("probit"), model = "random",
    effect = "individual", index = c("id", "time"), R = 12,
    method = "bfgs", print.level = 0, reltol = 1e-12,
    start = panel_probit_re_start(formula, dat)
  )

  y1 <- fit_model(y1 ~ x + z)
  y2 <- fit_model(y2 ~ x + z)
  for (i in seq_along(list(y1, y2))) {
    actual <- panel_probit_re_stata_model(list(y1, y2)[[i]])
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
    fixture$balanced$prediction_mean)), 1e-6)
  expect_lt(max(abs(slopes$estimate -
    fixture$balanced$average_slope_x)), 1e-6)

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
    balanced = list(fit = balanced, coefficient = 7e-4, covariance = 7e-4),
    partial_higher = list(
      fit = partial, coefficient = 6e-4, covariance = 6e-4),
    disjoint = list(fit = disjoint, coefficient = 2e-3, covariance = 2e-3)
  )
  for (case in names(cases)) {
    actual <- panel_probit_re_stata_system(cases[[case]]$fit)
    expected <- fixture[[case]]$system
    expect_lt(max(abs(actual$coefficients - expected$coefficients)),
      cases[[case]]$coefficient)
    expect_lt(max(abs(actual$vcov - expected$vcov)),
      cases[[case]]$covariance)
    expect_lt(max(abs(diag(actual$vcov)/diag(expected$vcov) - 1)), 0.005)
    expect_equal(actual[c("nobs_union", "nobs_overlap", "n_clusters")],
      expected[c("nobs_union", "nobs_overlap", "n_clusters")])
  }
  expect_equal(unname(vcov(disjoint)[1:4, 5:8]), matrix(0, 4L, 4L),
    tolerance = 1e-12)
})
