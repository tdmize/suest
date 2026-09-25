glmmtmb_nbinom2_ri_data <- function(groups = 70L, periods = 6L) {
  set.seed(229)
  id <- rep(seq_len(groups), each = periods)
  x <- rnorm(length(id))
  z <- rep(rnorm(groups), each = periods)
  woman <- rep(rep(0:1, length.out = groups), each = periods)
  common <- rep(rnorm(groups, sd = 0.6), each = periods)
  second <- rep(rnorm(groups, sd = 0.35), each = periods)
  data.frame(id = factor(id), time = rep(seq_len(periods), groups),
    higher = factor(ceiling(id/5)), x = x, z = z, woman = factor(woman),
    y1 = rnbinom(length(id), size = 1.7,
      mu = exp(0.4 + 0.35*x - 0.2*z + 0.15*woman + common)),
    y2 = rnbinom(length(id), size = 2.5,
      mu = exp(0.2 + 0.2*x + 0.25*z - 0.1*woman + 0.7*common + second)))
}

glmmtmb_nbinom2_ri_fit <- function(formula, data, ...) {
  glmmTMB::glmmTMB(formula, data = data,
    family = glmmTMB::nbinom2("log"), ...)
}

# Independent scalar Laplace calculation, including the NB2 dispersion score.
glmmtmb_nbinom2_ri_loglik <- function(parameters, model, panels = NULL) {
  panel <- as.character(model$frame$id)
  if (is.null(panels)) panels <- unique(panel)
  X <- model.matrix(model, component = "cond")
  eta <- drop(X %*% parameters[colnames(X)])
  y <- model.response(model.frame(model))
  phi <- exp(parameters["log_phi"])
  sigma2 <- exp(2*parameters["log_sigma"])
  sum(vapply(panels, function(value) {
    rows <- panel == value
    score <- function(b) {
      mu <- exp(eta[rows] + b)
      sum(phi*(y[rows] - mu)/(phi + mu)) - b/sigma2
    }
    bhat <- uniroot(score, c(-30, 30), tol = 1e-12)$root
    mu <- exp(eta[rows] + bhat)
    curvature <- sum(phi*mu*(phi + y[rows])/(phi + mu)^2) + 1/sigma2
    sum(dnbinom(y[rows], size = phi, mu = mu, log = TRUE)) -
      bhat^2/(2*sigma2) - parameters["log_sigma"] - 0.5*log(curvature)
  }, numeric(1)))
}

test_that("NB2 rejects the returned Stata boundary case despite nominal convergence", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  dat <- readRDS(test_path("fixtures", "glmmtmb-nbinom2-ri-boundary.rds"))
  model <- suppressWarnings(glmmtmb_nbinom2_ri_fit(y2 ~ x + z + (1 | id), dat))
  expect_lt(exp(2*model$fit$par["theta"]), 1e-6)
  expect_error(suest(model, model, observation_id = c("id", "time")),
    "zero boundary|positive-definite Hessian")
})

test_that("NB2 scores retain dispersion and reproduce the Laplace likelihood", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  dat <- glmmtmb_nbinom2_ri_data()
  model <- glmmtmb_nbinom2_ri_fit(y1 ~ x + z + (1 | id), dat)
  adapter <- suest:::.suest_model_adapter(model)
  component <- suest:::.suest_model_components(model, adapter$type, adapter$engine)
  parameters <- component$parameters
  expect_identical(names(parameters), c("(Intercept)", "x", "z", "log_phi", "log_sigma"))
  expect_equal(unname(parameters["log_phi"]), log(sigma(model)), tolerance = 1e-12)
  expect_equal(glmmtmb_nbinom2_ri_loglik(parameters, model),
    as.numeric(logLik(model)), tolerance = 1e-8)
  expect_equal(unname(component$bread/nrow(dat)),
    unname(as.matrix(vcov(model, full = TRUE))), tolerance = 1e-12)
  expect_lt(max(abs(colSums(component$score))), 5e-3)
  expect_equal(sum(rowSums(abs(component$score)) > 0), nlevels(dat$id))

  step <- .Machine$double.eps^(1/3)*pmax(1, abs(parameters))
  for (panel in c("1", "23", "70")) {
    numerical <- vapply(seq_along(parameters), function(j) {
      plus <- minus <- parameters
      plus[j] <- plus[j] + step[j]
      minus[j] <- minus[j] - step[j]
      (glmmtmb_nbinom2_ri_loglik(plus, model, panel) -
        glmmtmb_nbinom2_ri_loglik(minus, model, panel))/(2*step[j])
    }, numeric(1))
    expect_equal(unname(colSums(component$score[dat$id == panel, , drop = FALSE])),
      numerical, tolerance = 3e-6)
  }
})

test_that("NB2 reproduces the five completed Stata component fits", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  fixture <- readRDS(test_path("fixtures", "glmmtmb-nbinom2-ri-stata-components.rds"))
  dat <- fixture$data
  dat$id <- factor(dat$id)
  for (outcome in names(fixture$models)) {
    model <- glmmtmb_nbinom2_ri_fit(reformulate(c("x", "z", "(1 | id)"), outcome), dat)
    adapter <- suest:::.suest_model_adapter(model)
    component <- suest:::.suest_model_components(model, adapter$type, adapter$engine)
    expected <- fixture$models[[outcome]]
    expect_lt(max(abs(component$parameters - expected$coefficients)), 5e-5)
    expect_lt(max(abs(component$bread/nobs(model) - expected$vcov)), 1e-4)
    expect_lt(abs(as.numeric(logLik(model)) - expected$logLik), 1e-6)
  }
})

test_that("NB2 reproduces revision-2 Stata components and Laplace systems", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  skip_if_not_installed("marginaleffects")
  fixture <- readRDS(test_path("fixtures", "glmmtmb-nbinom2-ri-stata.rds"))
  dat <- fixture$data
  dat$id <- factor(dat$id)
  models <- lapply(names(fixture$components), function(outcome) {
    glmmtmb_nbinom2_ri_fit(reformulate(c("x", "z", "(1 | id)"), outcome),
      dat[!is.na(dat[[outcome]]), ])
  })
  names(models) <- names(fixture$components)
  for (outcome in names(models)) {
    model <- models[[outcome]]
    expected <- fixture$components[[outcome]]
    adapter <- suest:::.suest_model_adapter(model)
    component <- suest:::.suest_model_components(model, adapter$type, adapter$engine)
    native <- component$bread/nobs(model)
    expect_lt(max(abs(component$parameters - expected$coefficients)), 5e-5)
    expect_lt(max(abs(native - expected$vcov)), 1e-4)
    expect_lt(max(abs(sqrt(diag(expected$vcov)/diag(native)) - 1)), .0025)
    expect_lt(abs(as.numeric(logLik(model)) - expected$logLik), 1e-6)
    expect_equal(nobs(model), expected$nobs)
    if (outcome %in% c("y1", "y2")) {
      X <- model.matrix(model, component = "cond")
      b <- component$parameters
      sigma2 <- exp(2*b["log_sigma"])
      mu <- exp(drop(X %*% b[colnames(X)]) + sigma2/2)
      gradient <- c(colMeans(mu*X), 0, sigma2*mean(mu))
      slope_gradient <- gradient*b["x"]
      slope_gradient[2L] <- slope_gradient[2L] + mean(mu)
      G <- rbind(gradient, slope_gradient)
      expect_lt(max(abs(c(mean(mu), b["x"]*mean(mu)) -
        expected$native_margins$estimates)), 2e-5)
      expect_lt(max(abs(sqrt(diag(G %*% native %*% t(G))) -
        sqrt(expected$native_margins$variances))), 1e-4)
    }
  }
  systems <- list(
    balanced = suest(models$y1, models$y2, model_names = c("Y1", "Y2"),
      observation_id = c("id", "time")),
    partial_higher = suest(models$y1_partial, models$y2_partial,
      model_names = c("Y1", "Y2"), observation_id = c("id", "time"), cluster = "higher"),
    disjoint = suest(models$y1_left, models$y2_right, model_names = c("Y1", "Y2"),
      observation_id = c("id", "time")))
  covariance_tolerance <- c(balanced = 1e-4, partial_higher = 1e-4, disjoint = 2e-4)
  se_tolerance <- c(balanced = .002, partial_higher = .002, disjoint = .005)
  for (case in names(systems)) {
    actual <- systems[[case]]
    expected <- fixture$gsem[[case]]
    # GSEM counts the union's 80 groups in the disjoint case; the specialized
    # suest correction counts each model's 40 groups. Preserve the raw fixture.
    corrected <- expected$vcov * fixture$metadata$suest_correction[case] /
      fixture$metadata$gsem_correction[case]
    expect_lt(max(abs(coef(actual) - expected$coefficients)), 5e-5)
    expect_lt(max(abs(vcov(actual) - corrected)), covariance_tolerance[case])
    expect_lt(max(abs(sqrt(diag(corrected)/diag(vcov(actual))) - 1)), se_tolerance[case])
    expect_equal(actual$nobs_union, expected$nobs)
    expect_equal(actual$nobs_overlap, unname(fixture$metadata$nobs_overlap[case]))
    expect_identical(actual$n_clusters, expected$n_clusters)
    indices <- switch(case, balanced = 1:2, partial_higher = 3:4, disjoint = 5:6)
    expect_lt(abs(sum(vapply(models[indices], function(m) as.numeric(logLik(m)),
      numeric(1))) - expected$logLik), 1e-6)
  }
  # Check the disjoint diagonal blocks against independent group influences,
  # including the per-model correction rather than only zero cross-covariance.
  for (i in 1:2) {
    model <- models[[i+4L]]
    influence <- sandwich::estfun(model, full = TRUE) %*% as.matrix(vcov(model, full = TRUE))
    index <- ((i-1L)*5L+1L):(i*5L)
    expect_equal(unname(vcov(systems$disjoint)[index, index]),
      unname(crossprod(influence)*40/39), tolerance = 1e-9)
  }
  expect_equal(unname(offdiag_vcov(systems$disjoint)), matrix(0, 5L, 5L), tolerance = 1e-12)
  expect_equal(unname(fixture$adaptive$disjoint$vcov[1:5, 6:10]), matrix(0, 5L, 5L),
    tolerance = 1e-12)
  expect_lt(max(abs(fixture$gsem$disjoint$vcov[1:5, 6:10])), 1e-4)
  expect_identical(fixture$metadata$laplace_suest2_return_code, 322L)
  predictions <- marginaleffects::avg_predictions(systems$balanced, newdata = dat)
  slopes <- marginaleffects::avg_slopes(systems$balanced, variables = "x", newdata = dat)
  means <- vapply(fixture$components[1:2], function(x) x$native_margins$estimates["mean"], numeric(1))
  slopes_x <- vapply(fixture$components[1:2], function(x) x$native_margins$estimates["slope_x"], numeric(1))
  expect_lt(max(abs(predictions$estimate - means)), 2e-5)
  expect_lt(max(abs(slopes$estimate - slopes_x)), 2e-5)
})

test_that("NB2 joint covariance uses complete group influences and overlap", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  dat <- glmmtmb_nbinom2_ri_data()
  first <- glmmtmb_nbinom2_ri_fit(y1 ~ x + z + (1 | id), dat)
  second <- glmmtmb_nbinom2_ri_fit(y2 ~ x + z + (1 | id), dat)
  fit <- suest(first, second, observation_id = c("id", "time"))
  influence <- function(model) sandwich::estfun(model, full = TRUE) %*%
    as.matrix(vcov(model, full = TRUE))
  joint <- cbind(influence(first), influence(second))
  expect_equal(unname(vcov(fit)), unname(crossprod(joint)*70/69), tolerance = 1e-9)
  expect_identical(fit$n_clusters, 70L)
  expect_gt(max(abs(offdiag_vcov(fit))), 0)
  higher <- suest(first, second, observation_id = c("id", "time"), cluster = "higher")
  grouped <- rowsum(joint, ceiling(as.integer(rownames(joint))/5))
  expect_equal(unname(vcov(higher)), unname(crossprod(grouped)*14/13), tolerance = 1e-9)

  partial1 <- dat[!(dat$time == 6 & as.integer(dat$id) %% 4 == 0), ]
  partial2 <- dat[!(dat$time == 1 & as.integer(dat$id) %% 5 == 0), ]
  partial <- suest(glmmtmb_nbinom2_ri_fit(y1 ~ x + z + (1 | id), partial1),
    glmmtmb_nbinom2_ri_fit(y2 ~ x + z + (1 | id), partial2),
    observation_id = c("id", "time"), cluster = "higher")
  expect_equal(partial$nobs_union, nrow(dat))
  expect_equal(partial$nobs_overlap, nrow(dat) - 31L)
  expect_true(all(is.finite(vcov(partial))))
  expect_gt(max(abs(offdiag_vcov(partial))), 0)

  low <- glmmtmb_nbinom2_ri_fit(y1 ~ x + z + (1 | id),
    dat[as.integer(dat$id) <= 35, ])
  high <- glmmtmb_nbinom2_ri_fit(y2 ~ x + z + (1 | id),
    dat[as.integer(dat$id) > 35, ])
  disjoint <- suest(low, high, observation_id = c("id", "time"))
  expect_equal(unname(offdiag_vcov(disjoint)), matrix(0, 5L, 5L), tolerance = 1e-12)
  expect_error(suest(first, second, observation_id = c("id", "time"), cluster = "time"),
    "Panel IDs must be nested within clusters", fixed = TRUE)
})

test_that("NB2 population means, slopes, and contrasts retain variance uncertainty", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  dat <- glmmtmb_nbinom2_ri_data()
  first <- glmmtmb_nbinom2_ri_fit(y1 ~ x*woman + z + (1 | id), dat)
  second <- glmmtmb_nbinom2_ri_fit(y2 ~ x*woman + z + (1 | id), dat)
  fit <- suest(first, second, model_names = c("First", "Second"),
    observation_id = c("id", "time"))
  means <- slopes <- numeric(2)
  gradients <- matrix(0, 2, length(coef(fit)))
  slope_gradients <- comparison_gradients <- gradients
  comparisons <- numeric(2)
  for (i in 1:2) {
    model <- list(first, second)[[i]]
    X <- model.matrix(model, component = "cond")
    beta <- glmmTMB::fixef(model)$cond
    sigma2 <- as.numeric(glmmTMB::VarCorr(model)$cond$id[1, 1])
    mu <- exp(drop(X %*% beta) + sigma2/2)
    means[i] <- mean(mu)
    multiplier <- beta["x"] + beta["x:woman1"]*(dat$woman == "1")
    slopes[i] <- mean(mu*multiplier)
    indices <- ((i - 1)*7 + 1):(i*7)
    gradients[i, indices] <- c(colMeans(mu*X), 0, sigma2*mean(mu))
    direct <- rep(0, length(beta))
    direct[match("x", names(beta))] <- mean(mu)
    direct[match("x:woman1", names(beta))] <- mean(mu*(dat$woman == "1"))
    slope_gradients[i, indices] <- c(colMeans(mu*multiplier*X) + direct,
      0, sigma2*slopes[i])
    X0 <- X1 <- X
    X0[, c("woman1", "x:woman1")] <- 0
    X1[, "woman1"] <- 1
    X1[, "x:woman1"] <- dat$x
    mu0 <- exp(drop(X0 %*% beta) + sigma2/2)
    mu1 <- exp(drop(X1 %*% beta) + sigma2/2)
    comparisons[i] <- mean(mu1 - mu0)
    comparison_gradients[i, indices] <- c(colMeans(mu1*X1 - mu0*X0),
      0, sigma2*comparisons[i])
  }
  predictions <- marginaleffects::avg_predictions(fit, newdata = dat)
  effects <- marginaleffects::avg_slopes(fit, variables = "x", newdata = dat)
  differences <- marginaleffects::avg_comparisons(fit, variables = "woman", newdata = dat)
  expected_vcov <- gradients %*% vcov(fit) %*% t(gradients)
  expect_equal(predictions$estimate, means, tolerance = 1e-8)
  expect_equal(predictions$std.error, sqrt(diag(expected_vcov)), tolerance = 2e-6)
  expect_equal(effects$estimate, slopes, tolerance = 1e-6)
  expect_equal(effects$std.error,
    sqrt(diag(slope_gradients %*% vcov(fit) %*% t(slope_gradients))), tolerance = 2e-5)
  expect_equal(differences$estimate, comparisons, tolerance = 1e-8)
  expect_equal(differences$std.error,
    sqrt(diag(comparison_gradients %*% vcov(fit) %*% t(comparison_gradients))),
    tolerance = 2e-6)
  expect_identical(as.character(predictions$group), c("First", "Second"))
  contrast <- suppressWarnings(marginaleffects::hypotheses(predictions, "b1 - b2 = 0"))
  expect_equal(as.numeric(contrast$estimate), means[1] - means[2], tolerance = 1e-8)
  expect_equal(contrast$std.error,
    sqrt(expected_vcov[1, 1] + expected_vcov[2, 2] - 2*expected_vcov[1, 2]),
    tolerance = 2e-6)

  changed <- coef(fit)
  changed["First::log_phi"] <- changed["First::log_phi"] + log(2)
  altered <- marginaleffects::set_coef(fit, changed)
  expect_equal(marginaleffects::avg_predictions(altered, newdata = dat, vcov = FALSE)$estimate,
    means, tolerance = 1e-8)
  changed["First::log_sigma"] <- changed["First::log_sigma"] + log(1.5)
  altered <- marginaleffects::set_coef(fit, changed)
  expect_gt(marginaleffects::avg_predictions(altered, newdata = dat, vcov = FALSE)$estimate[1],
    means[1])
  expect_equal(marginaleffects::avg_predictions(altered, newdata = dat, vcov = FALSE)$estimate[2],
    means[2], tolerance = 1e-8)
})

test_that("NB2 excludes unvalidated dispersion and random-effect specifications", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  dat <- glmmtmb_nbinom2_ri_data()
  base <- glmmtmb_nbinom2_ri_fit(y1 ~ x + z + (1 | id), dat)
  dispersion <- glmmtmb_nbinom2_ri_fit(y1 ~ x + z + (1 | id), dat, dispformula = ~x)
  dispersion_offset <- glmmtmb_nbinom2_ri_fit(y1 ~ x + z + (1 | id), dat,
    dispformula = ~1 + offset(x))
  fixed_dispersion <- glmmtmb_nbinom2_ri_fit(y1 ~ x + z + (1 | id), dat,
    start = list(betadisp = log(1.7)), map = list(betadisp = factor(NA)))
  for (model in list(dispersion, dispersion_offset, fixed_dispersion))
    expect_error(suest(model, base), "dispersion", ignore.case = TRUE)
  constrained <- glmmtmb_nbinom2_ri_fit(y1 ~ x + z + (1 | id), dat,
    map = list(beta = factor(c(1, 2, 2))))
  expect_error(suest(constrained, base), "Mapped or constrained", fixed = TRUE)
  permuted <- glmmtmb_nbinom2_ri_fit(y1 ~ x + z + (1 | id), dat,
    map = list(beta = factor(c(3, 2, 1))))
  expect_error(suest(permuted, base), "Mapped or constrained", fixed = TRUE)
  dat$log_phi <- dat$x
  reserved <- glmmtmb_nbinom2_ri_fit(y1 ~ log_phi + z + (1 | id), dat)
  expect_error(suest(reserved, base), "reserves coefficient names", fixed = TRUE)
  zi <- glmmtmb_nbinom2_ri_fit(y1 ~ x + z + (1 | id), dat, ziformula = ~1)
  weighted <- glmmtmb_nbinom2_ri_fit(y1 ~ x + z + (1 | id), dat, weights = rep(2, nrow(dat)))
  offset <- glmmtmb_nbinom2_ri_fit(y1 ~ x + z + offset(x) + (1 | id), dat)
  dat$log_exposure <- 0
  zero_offset <- glmmtmb_nbinom2_ri_fit(
    y1 ~ x + z + offset(log_exposure) + (1 | id), dat)
  zero_call_offset <- glmmTMB::glmmTMB(y1 ~ x + z + (1 | id), data = dat,
    family = glmmTMB::nbinom2("log"), offset = log_exposure)
  slope <- suppressWarnings(glmmtmb_nbinom2_ri_fit(y1 ~ x + z + (1 + x + z | id), dat))
  poisson <- glmmTMB::glmmTMB(y1 ~ x + z + (1 | id), data = dat, family = poisson())
  expect_error(suest(zi, base), "Zero-inflation", fixed = TRUE)
  expect_error(suest(weighted, base), "Weighted glmmTMB", fixed = TRUE)
  expect_error(suest(offset, base), "Offsets are not supported", fixed = TRUE)
  expect_error(suest(zero_offset, base), "Offsets are not supported", fixed = TRUE)
  expect_error(suest(zero_call_offset, base), "Offsets are not supported", fixed = TRUE)
  expect_error(suest(slope, base), "one conditional random intercept", fixed = TRUE)
  expect_error(suest(base, poisson), "same type", fixed = TRUE)
})
