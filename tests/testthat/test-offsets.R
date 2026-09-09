offset_test_data <- function(n = 700L) {
  set.seed(4802)
  x <- rnorm(n)
  z <- rnorm(n)
  offset <- 0.3*rnorm(n)
  eta <- -0.2 + 0.5*x - 0.25*z + offset
  latent <- eta + rlogis(n)
  data.frame(
    x = x,
    z = z,
    offset = offset,
    exposure = exp(offset),
    weight = exp(0.15*z),
    y_lm = eta + rnorm(n),
    y_bin = rbinom(n, 1, plogis(eta)),
    y_count = rpois(n, exp(0.4 + 0.2*x - 0.1*z + offset)),
    y_ord = ordered(
      cut(latent, c(-Inf, -0.8, 0, 0.8, Inf), labels = FALSE),
      levels = 1:4,
      labels = c("A", "B", "C", "D")
    )
  )
}

test_that("offsets work across scalar model families", {
  dat <- offset_test_data()
  linear <- lm(y_lm ~ x + z + offset(offset), data = dat)
  logit <- glm(
    y_bin ~ x + z + offset(offset),
    family = binomial(),
    data = dat
  )
  poisson <- glm(
    y_count ~ x + z + offset(log(exposure)),
    family = poisson(),
    data = dat
  )
  fit <- suest(
    linear,
    logit,
    poisson,
    model_names = c("Linear", "Logit", "Poisson")
  )
  effects <- marginaleffects::avg_comparisons(
    fit,
    variables = "x",
    newdata = dat
  )

  for (i in seq_along(fit$models)) {
    block <- fit$index[[i]]
    if (fit$model_types[i] == "lm")
      block <- block[fit$local_names[[i]] != "lnvar"]
    expect_equal(
      unname(fit$vcov[block, block, drop = FALSE]),
      unname(robust_vcov_direct(fit$models[[i]])),
      tolerance = 1e-8
    )
  }
  expect_equal(nrow(effects), 3L)
  expect_true(all(is.finite(effects$std.error)))
})

test_that("pweighted offset models are scale invariant", {
  dat <- offset_test_data()
  linear <- lm(
    y_lm ~ x + z + offset(offset),
    weights = weight,
    data = dat
  )
  logit <- suppressWarnings(glm(
    y_bin ~ x + z + offset(offset),
    family = binomial(),
    weights = weight,
    data = dat
  ))
  poisson <- glm(
    y_count ~ x + z + offset(log(exposure)),
    family = poisson(),
    weights = weight,
    data = dat
  )
  models <- list(linear, logit, poisson)
  scaled_models <- lapply(models, function(model) {
    model$model[["(weights)"]] <- 10*model$model[["(weights)"]]
    if (!is.null(model$prior.weights))
      model$prior.weights <- 10*model$prior.weights
    model
  })
  fit <- suest(
    linear,
    logit,
    poisson,
    model_names = c("Linear", "Logit", "Poisson"),
    weight_type = "pweight"
  )
  scaled <- suest(
    scaled_models[[1L]],
    scaled_models[[2L]],
    scaled_models[[3L]],
    model_names = c("Linear", "Logit", "Poisson"),
    weight_type = "pweight"
  )

  mean_parameters <- !endsWith(names(coef(fit)), "::lnvar")
  expect_equal(
    coef(scaled)[mean_parameters],
    coef(fit)[mean_parameters],
    tolerance = 1e-12
  )
  expect_equal(
    vcov(scaled)[mean_parameters, mean_parameters],
    vcov(fit)[mean_parameters, mean_parameters],
    tolerance = 1e-12
  )
})

test_that("ordered-model offsets use the correct analytic scores", {
  dat <- offset_test_data()
  polr <- MASS::polr(
    y_ord ~ x + z + offset(offset),
    data = dat,
    method = "logistic",
    Hess = TRUE
  )
  clm <- ordinal::clm(
    y_ord ~ x + z + offset(offset),
    data = dat,
    link = "logit",
    Hess = TRUE,
    model = TRUE
  )
  components <- suest:::.suest_polr_components(polr)
  fit <- suest(polr, clm, model_names = c("polr", "clm"))
  effects <- marginaleffects::avg_comparisons(
    fit,
    variables = "x",
    newdata = dat
  )

  model_frame <- model.frame(polr)
  X <- model.matrix(polr)[, -1L, drop = FALSE]
  y <- as.integer(model.response(model_frame))
  model_offset <- model.offset(model_frame)
  parameters <- c(polr$coefficients, polr$zeta)
  loglik <- function(parameters) {
    beta <- parameters[seq_len(ncol(X))]
    zeta <- parameters[ncol(X) + seq_along(polr$zeta)]
    thresholds <- c(-Inf, zeta, Inf)
    eta <- as.numeric(X %*% beta) + model_offset
    sum(log(
      plogis(thresholds[y + 1L] - eta) -
        plogis(thresholds[y] - eta)
    ))
  }
  epsilon <- 1e-5*pmax(1, abs(parameters))
  numerical_score <- vapply(seq_along(parameters), function(i) {
    upper <- lower <- parameters
    upper[i] <- upper[i] + epsilon[i]
    lower[i] <- lower[i] - epsilon[i]
    (loglik(upper) - loglik(lower))/(2*epsilon[i])
  }, numeric(1))

  expect_lt(max(abs(
    unname(colSums(components$score)) - numerical_score
  )), 1e-6)
  expect_equal(nrow(effects), 8L)
  expect_true(all(is.finite(effects$std.error)))
})
