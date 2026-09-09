pweight_lm_vcov <- function(model, n_union = stats::nobs(model)) {
  mf <- stats::model.frame(model)
  X <- stats::model.matrix(model)
  weights <- stats::model.weights(mf)
  residual <- stats::model.response(mf) - as.numeric(X %*% stats::coef(model))
  inverse <- solve(crossprod(X, X * weights))

  inverse %*%
    crossprod(X * (weights * residual)) %*%
    inverse * n_union / (n_union - 1)
}

pweight_lm_cross_vcov <- function(model1, model2) {
  component <- function(model) {
    mf <- stats::model.frame(model)
    X <- stats::model.matrix(model)
    weights <- stats::model.weights(mf)
    residual <- stats::model.response(mf) -
      as.numeric(X %*% stats::coef(model))

    list(
      rows = rownames(mf),
      score = X * (weights * residual),
      inverse = solve(crossprod(X, X * weights)),
      n = nrow(X)
    )
  }

  first <- component(model1)
  second <- component(model2)
  overlap <- intersect(first$rows, second$rows)
  first_rows <- match(overlap, first$rows)
  second_rows <- match(overlap, second$rows)

  n_union <- length(union(first$rows, second$rows))
  correction <- n_union / (n_union - 1)

  first$inverse %*%
    crossprod(
      first$score[first_rows, , drop = FALSE],
      second$score[second_rows, , drop = FALSE]
    ) %*%
    second$inverse * correction
}

pweight_lm_mean_index <- function(object, model) {
  object$index[[model]][object$local_names[[model]] != "lnvar"]
}


pweight_binary_component <- function(model) {
  mf <- stats::model.frame(model)
  X <- stats::model.matrix(model)
  y <- stats::model.response(mf)
  y <- if (is.factor(y)) as.integer(y) - 1L else as.numeric(y)
  weights <- stats::model.weights(mf)
  eta <- as.numeric(X %*% stats::coef(model))
  link <- model$family$link

  if (link == "logit") {
    probability <- stats::plogis(eta)
    score_eta <- y - probability
    information_eta <- probability * (1 - probability)
  } else {
    probability <- stats::pnorm(eta)
    density <- stats::dnorm(eta)
    variance <- probability * (1 - probability)
    residual <- y - probability
    score_eta <- residual * density / variance
    information_eta <- density * (
      ((eta * residual + density) * variance +
         residual * density * (1 - 2 * probability)) / variance^2
    )
  }

  information <- crossprod(X, X * (weights * information_eta))
  list(
    rows = rownames(mf),
    score = X * (weights * score_eta),
    inverse = solve(information)
  )
}

pweight_binary_vcov <- function(model, n_union = stats::nobs(model)) {
  component <- pweight_binary_component(model)
  component$inverse %*% crossprod(component$score) %*%
    component$inverse * n_union / (n_union - 1)
}

pweight_poisson_component <- function(model) {
  mf <- stats::model.frame(model)
  X <- stats::model.matrix(model)
  y <- as.numeric(stats::model.response(mf))
  weights <- stats::model.weights(mf)
  mu <- exp(as.numeric(X %*% stats::coef(model)))

  list(
    rows = rownames(mf),
    score = X * (weights * (y - mu)),
    inverse = solve(crossprod(X, X * (weights * mu)))
  )
}

pweight_poisson_vcov <- function(model, n_union = stats::nobs(model)) {
  component <- pweight_poisson_component(model)
  component$inverse %*% crossprod(component$score) %*%
    component$inverse * n_union / (n_union - 1)
}

pweight_cross_vcov <- function(model1, model2, component1, component2) {
  first <- component1(model1)
  second <- component2(model2)
  overlap <- intersect(first$rows, second$rows)
  first_rows <- match(overlap, first$rows)
  second_rows <- match(overlap, second$rows)
  n_union <- length(union(first$rows, second$rows))

  first$inverse %*%
    crossprod(
      first$score[first_rows, , drop = FALSE],
      second$score[second_rows, , drop = FALSE]
    ) %*%
    second$inverse * n_union / (n_union - 1)
}

pweight_lm_component <- function(model) {
  mf <- stats::model.frame(model)
  X <- stats::model.matrix(model)
  weights <- stats::model.weights(mf)
  residual <- stats::model.response(mf) -
    as.numeric(X %*% stats::coef(model))

  list(
    rows = rownames(mf),
    score = X * (weights * residual),
    inverse = solve(crossprod(X, X * weights))
  )
}

pweight_test_data <- function(n = 400L) {
  set.seed(3762026)
  dat <- data.frame(
    id = seq_len(n),
    x = rnorm(n),
    z = rnorm(n)
  )
  dat$mediator <- 0.5 * dat$x - 0.3 * dat$z + rnorm(n)
  dat$w <- exp(0.25 * dat$z + 0.15 * rnorm(n))
  dat$w_copy <- dat$w
  dat$w10 <- 10 * dat$w
  dat$y <- 1 + 0.7 * dat$x - 0.4 * dat$z +
    0.6 * dat$mediator + rnorm(n)
  eta <- -0.3 + 0.8 * dat$x - 0.35 * dat$z + 0.5 * dat$mediator
  dat$binary <- stats::rbinom(n, 1, stats::plogis(eta))
  mu <- exp(0.2 + 0.35 * dat$x - 0.2 * dat$z + 0.25 * dat$mediator)
  dat$count <- stats::rpois(n, mu)
  dat
}

testthat::test_that("pweight lm covariance matches the direct sandwich", {
  dat <- pweight_test_data()
  model1 <- lm(y ~ x + z, weights = w, data = dat)
  model2 <- lm(y ~ x + z + mediator, weights = w, data = dat)
  fit <- suest(
    model1,
    model2,
    model_names = c("Base", "Adjusted"),
    weight_type = "pweight"
  )

  testthat::expect_identical(fit$weight_type, "pweight")
  first <- pweight_lm_mean_index(fit, 1L)
  second <- pweight_lm_mean_index(fit, 2L)
  testthat::expect_equal(
    unname(vcov(fit)[first, first]),
    unname(pweight_lm_vcov(model1)),
    tolerance = 1e-10
  )
  testthat::expect_equal(
    unname(vcov(fit)[second, second]),
    unname(pweight_lm_vcov(model2)),
    tolerance = 1e-10
  )
  testthat::expect_equal(
    unname(vcov(fit)[first, second]),
    unname(pweight_lm_cross_vcov(model1, model2)),
    tolerance = 1e-10
  )
})

testthat::test_that("pweight lm mean results are invariant to common scaling", {
  dat <- pweight_test_data()

  model1 <- lm(y ~ x + z, weights = w, data = dat)
  model2 <- lm(y ~ x + z + mediator, weights = w, data = dat)
  fit <- suest(
    model1,
    model2,
    model_names = c("Base", "Adjusted"),
    weight_type = "pweight"
  )

  model1_scaled <- lm(y ~ x + z, weights = w10, data = dat)
  model2_scaled <- lm(y ~ x + z + mediator, weights = w10, data = dat)
  fit_scaled <- suest(
    model1_scaled,
    model2_scaled,
    model_names = c("Base", "Adjusted"),
    weight_type = "pweight"
  )

  mean_parameters <- !endsWith(names(coef(fit)), "::lnvar")
  testthat::expect_equal(
    coef(fit_scaled)[mean_parameters],
    coef(fit)[mean_parameters],
    tolerance = 1e-12
  )
  testthat::expect_equal(
    vcov(fit_scaled)[mean_parameters, mean_parameters],
    vcov(fit)[mean_parameters, mean_parameters],
    tolerance = 1e-10
  )

  weighted_lnvar <- function(model) {
    mf <- stats::model.frame(model)
    X <- stats::model.matrix(model)
    weights <- stats::model.weights(mf)
    residual <- stats::model.response(mf) - as.numeric(X %*% stats::coef(model))
    log(sum(weights * residual^2) / (sum(weights) - ncol(X)))
  }
  ancillary <- endsWith(names(coef(fit)), "::lnvar")
  testthat::expect_equal(
    unname(coef(fit)[ancillary]),
    c(weighted_lnvar(model1), weighted_lnvar(model2)),
    tolerance = 1e-12
  )
  testthat::expect_equal(
    unname(coef(fit_scaled)[ancillary]),
    c(weighted_lnvar(model1_scaled), weighted_lnvar(model2_scaled)),
    tolerance = 1e-12
  )
  testthat::expect_false(isTRUE(all.equal(
    coef(fit_scaled)[ancillary],
    coef(fit)[ancillary],
    tolerance = 1e-12
  )))
})

testthat::test_that("pweights may differ outside a partially overlapping sample", {
  dat <- pweight_test_data()
  dat$w_b <- dat$w
  dat$w_b[dat$id > 300] <- dat$w[dat$id > 300] *
    (1.25 + 0.1 * abs(dat$x[dat$id > 300]))

  model1 <- lm(
    y ~ x + z,
    weights = w,
    data = dat,
    subset = id <= 300
  )
  model2 <- lm(
    y ~ x + z + mediator,
    weights = w_b,
    data = dat,
    subset = id >= 101
  )
  fit <- suest(model1, model2, weight_type = "pweight")

  testthat::expect_equal(fit$nobs_overlap, 200)
  testthat::expect_equal(fit$nobs_union, 400)
  first <- pweight_lm_mean_index(fit, 1L)
  second <- pweight_lm_mean_index(fit, 2L)
  testthat::expect_equal(
    unname(vcov(fit)[first, first]),
    unname(pweight_lm_vcov(model1, n_union = fit$nobs_union)),
    tolerance = 1e-10
  )
  testthat::expect_equal(
    unname(vcov(fit)[second, second]),
    unname(pweight_lm_vcov(model2, n_union = fit$nobs_union)),
    tolerance = 1e-10
  )
  testthat::expect_equal(
    unname(vcov(fit)[first, second]),
    unname(pweight_lm_cross_vcov(model1, model2)),
    tolerance = 1e-10
  )
})

testthat::test_that("different pweights are accepted for disjoint samples", {
  dat <- pweight_test_data()
  dat$w_left <- exp(0.2 * dat$x - 0.1 * dat$z)
  dat$w_right <- exp(-0.15 * dat$x + 0.3 * dat$z)

  model1 <- lm(
    y ~ x + z,
    weights = w_left,
    data = dat,
    subset = id <= 200
  )
  model2 <- lm(
    y ~ x + z + mediator,
    weights = w_right,
    data = dat,
    subset = id > 200
  )
  fit <- suest(model1, model2, weight_type = "pweight")

  testthat::expect_equal(fit$nobs_overlap, 0)
  testthat::expect_true(all(
    unname(vcov(fit)[fit$index[[1]], fit$index[[2]]]) == 0
  ))
})

testthat::test_that("pweight validation compares values on overlapping rows", {
  dat <- pweight_test_data()

  same1 <- lm(y ~ x + z, weights = w, data = dat)
  same2 <- lm(y ~ x + z + mediator, weights = w_copy, data = dat)
  same_fit <- suest(same1, same2, weight_type = "pweight")
  testthat::expect_s3_class(same_fit, "suest_model")

  dat$w_bad <- dat$w
  dat$w_bad[150] <- 1.1 * dat$w_bad[150]
  bad <- lm(y ~ x + z + mediator, weights = w_bad, data = dat)

  testthat::expect_error(
    suest(same1, bad, weight_type = "pweight"),
    "must agree for observations included in both models"
  )

  scaled <- lm(y ~ x + z + mediator, weights = w10, data = dat)
  testthat::expect_error(
    suest(same1, scaled, weight_type = "pweight"),
    "must agree for observations included in both models"
  )
})

testthat::test_that("pweight requests are explicit", {
  dat <- pweight_test_data()
  model1 <- lm(y ~ x, weights = w, data = dat)
  model2 <- lm(y ~ x + z, weights = w, data = dat)

  testthat::expect_error(
    suest(model1, model2),
    'weight_type = "pweight"'
  )
  testthat::expect_error(
    suest(model1, model2, weight_type = "sampling"),
    "must be NULL"
  )

  poisson1 <- glm(
    count ~ x,
    family = poisson(),
    weights = w,
    data = dat
  )
  poisson2 <- glm(
    count ~ x + z,
    family = poisson(),
    weights = w,
    data = dat
  )
  testthat::expect_s3_class(
    suest(poisson1, poisson2, weight_type = "pweight"),
    "suest_model"
  )

  testthat::skip_if_not_installed("MASS")
  nb1 <- suppressWarnings(MASS::glm.nb(count ~ x, weights = w, data = dat))
  nb2 <- suppressWarnings(MASS::glm.nb(count ~ x + z, weights = w, data = dat))
  testthat::expect_s3_class(
    suest(nb1, nb2, weight_type = "pweight"),
    "suest_model"
  )
})

testthat::test_that("pweights must be strictly positive", {
  dat <- pweight_test_data()
  dat$w_zero <- dat$w
  dat$w_zero[1] <- 0

  model1 <- lm(y ~ x, weights = w_zero, data = dat)
  model2 <- lm(y ~ x + z, weights = w_zero, data = dat)

  testthat::expect_error(
    suest(model1, model2, weight_type = "pweight"),
    "strictly positive"
  )
})


testthat::test_that("pweight logit and probit covariance matches direct OIM sandwiches", {
  dat <- pweight_test_data()

  for (link in c("logit", "probit")) {
    model1 <- suppressWarnings(glm(
      binary ~ x + z,
      family = binomial(link = link),
      weights = w,
      data = dat
    ))
    model2 <- suppressWarnings(glm(
      binary ~ x + z + mediator,
      family = binomial(link = link),
      weights = w,
      data = dat
    ))
    fit <- suest(
      model1,
      model2,
      model_names = c("Base", "Adjusted"),
      weight_type = "pweight"
    )

    testthat::expect_equal(
      unname(vcov(fit)[fit$index[[1]], fit$index[[1]]]),
      unname(pweight_binary_vcov(model1)),
      tolerance = 1e-9
    )
    testthat::expect_equal(
      unname(vcov(fit)[fit$index[[2]], fit$index[[2]]]),
      unname(pweight_binary_vcov(model2)),
      tolerance = 1e-9
    )
    testthat::expect_equal(
      unname(vcov(fit)[fit$index[[1]], fit$index[[2]]]),
      unname(pweight_cross_vcov(
        model1,
        model2,
        pweight_binary_component,
        pweight_binary_component
      )),
      tolerance = 1e-9
    )
  }
})



testthat::test_that("pweighted binary models support two-level factor outcomes", {
  dat <- pweight_test_data()
  dat$binary_factor <- factor(dat$binary, levels = c(0, 1))
  numeric1 <- suppressWarnings(glm(
    binary ~ x + z,
    family = binomial(),
    weights = w,
    data = dat
  ))
  numeric2 <- suppressWarnings(glm(
    binary ~ x + z + mediator,
    family = binomial(),
    weights = w,
    data = dat
  ))
  factor1 <- suppressWarnings(glm(
    binary_factor ~ x + z,
    family = binomial(),
    weights = w,
    data = dat
  ))
  factor2 <- suppressWarnings(glm(
    binary_factor ~ x + z + mediator,
    family = binomial(),
    weights = w,
    data = dat
  ))
  numeric_fit <- suest(
    numeric1,
    numeric2,
    model_names = c("Base", "Adjusted"),
    weight_type = "pweight"
  )
  factor_fit <- suest(
    factor1,
    factor2,
    model_names = c("Base", "Adjusted"),
    weight_type = "pweight"
  )

  testthat::expect_equal(coef(factor_fit), coef(numeric_fit), tolerance = 1e-8)
  testthat::expect_equal(vcov(factor_fit), vcov(numeric_fit), tolerance = 1e-8)
})

testthat::test_that("glm2 pweighted binary models match stats glm", {
  testthat::skip_if_not_installed("glm2")
  dat <- pweight_test_data()
  stats1 <- suppressWarnings(glm(
    binary ~ x + z,
    family = binomial(link = "probit"),
    weights = w,
    data = dat
  ))
  stats2 <- suppressWarnings(glm(
    binary ~ x + z + mediator,
    family = binomial(link = "probit"),
    weights = w,
    data = dat
  ))
  glm21 <- suppressWarnings(glm2::glm2(
    binary ~ x + z,
    family = binomial(link = "probit"),
    weights = w,
    data = dat
  ))
  glm22 <- suppressWarnings(glm2::glm2(
    binary ~ x + z + mediator,
    family = binomial(link = "probit"),
    weights = w,
    data = dat
  ))
  stats_fit <- suest(
    stats1,
    stats2,
    model_names = c("Base", "Adjusted"),
    weight_type = "pweight"
  )
  glm2_fit <- suest(
    glm21,
    glm22,
    model_names = c("Base", "Adjusted"),
    weight_type = "pweight"
  )

  testthat::expect_equal(coef(glm2_fit), coef(stats_fit), tolerance = 1e-7)
  testthat::expect_equal(vcov(glm2_fit), vcov(stats_fit), tolerance = 1e-7)
})

testthat::test_that("pweighted probit uses observed rather than Fisher information", {
  dat <- pweight_test_data()
  model <- suppressWarnings(glm(
    binary ~ x + z + mediator,
    family = binomial(link = "probit"),
    weights = w,
    data = dat
  ))
  fit <- suest(
    model,
    model,
    model_names = c("First", "Second"),
    weight_type = "pweight"
  )
  observed <- pweight_binary_vcov(model)
  fisher <- sandwich::sandwich(model) * nrow(dat) / (nrow(dat) - 1)

  testthat::expect_equal(
    unname(vcov(fit)[fit$index[[1]], fit$index[[1]]]),
    unname(observed),
    tolerance = 1e-9
  )
  testthat::expect_gt(max(abs(observed - fisher)), 1e-6)
})

testthat::test_that("pweights support binary overlap, disjoint samples, and mixed families", {
  dat <- pweight_test_data()
  dat$w_b <- dat$w
  dat$w_b[dat$id > 300] <- dat$w[dat$id > 300] *
    (1.25 + 0.1 * abs(dat$x[dat$id > 300]))
  dat$w_left <- exp(0.2 * dat$x - 0.1 * dat$z)
  dat$w_right <- exp(-0.15 * dat$x + 0.3 * dat$z)

  partial1 <- suppressWarnings(glm(
    binary ~ x + z,
    family = binomial(),
    weights = w,
    data = dat,
    subset = id <= 300
  ))
  partial2 <- suppressWarnings(glm(
    binary ~ x + z + mediator,
    family = binomial(),
    weights = w_b,
    data = dat,
    subset = id >= 101
  ))
  partial_fit <- suest(partial1, partial2, weight_type = "pweight")
  testthat::expect_equal(partial_fit$nobs_overlap, 200)

  disjoint1 <- suppressWarnings(glm(
    binary ~ x + z,
    family = binomial(link = "probit"),
    weights = w_left,
    data = dat,
    subset = id <= 200
  ))
  disjoint2 <- suppressWarnings(glm(
    binary ~ x + z + mediator,
    family = binomial(link = "probit"),
    weights = w_right,
    data = dat,
    subset = id > 200
  ))
  disjoint_fit <- suest(disjoint1, disjoint2, weight_type = "pweight")
  testthat::expect_equal(disjoint_fit$nobs_overlap, 0)
  testthat::expect_true(all(
    unname(vcov(disjoint_fit)[disjoint_fit$index[[1]],
                               disjoint_fit$index[[2]]]) == 0
  ))

  linear <- lm(y ~ x + z, weights = w, data = dat)
  logit <- suppressWarnings(glm(
    binary ~ x + z,
    family = binomial(),
    weights = w,
    data = dat
  ))
  mixed <- suest(
    linear,
    logit,
    model_names = c("Linear", "Logit"),
    weight_type = "pweight"
  )
  linear_mean <- pweight_lm_mean_index(mixed, 1L)
  testthat::expect_equal(
    unname(vcov(mixed)[linear_mean, mixed$index[[2]]]),
    unname(pweight_cross_vcov(
      linear,
      logit,
      pweight_lm_component,
      pweight_binary_component
    )),
    tolerance = 1e-9
  )
})

testthat::test_that("pweight binary results are invariant to common scaling", {
  dat <- pweight_test_data()
  model1 <- suppressWarnings(glm(
    binary ~ x + z,
    family = binomial(link = "probit"),
    weights = w,
    data = dat
  ))
  model2 <- suppressWarnings(glm(
    binary ~ x + z + mediator,
    family = binomial(link = "probit"),
    weights = w,
    data = dat
  ))

  scale_stored_weights <- function(model, multiplier) {
    scaled <- model
    scaled$prior.weights <- multiplier * model$prior.weights
    scaled$model[["(weights)"]] <- multiplier * model$model[["(weights)"]]
    scaled
  }

  scaled1 <- scale_stored_weights(model1, 10)
  scaled2 <- scale_stored_weights(model2, 10)
  fit <- suest(
    model1,
    model2,
    model_names = c("Base", "Adjusted"),
    weight_type = "pweight"
  )
  scaled_fit <- suest(
    scaled1,
    scaled2,
    model_names = c("Base", "Adjusted"),
    weight_type = "pweight"
  )

  testthat::expect_equal(
    stats::model.weights(stats::model.frame(scaled1)),
    10 * stats::model.weights(stats::model.frame(model1))
  )
  testthat::expect_equal(coef(scaled_fit), coef(fit), tolerance = 1e-12)
  testthat::expect_equal(vcov(scaled_fit), vcov(fit), tolerance = 1e-12)
})

testthat::test_that("pweight Poisson covariance matches the direct sandwich", {
  dat <- pweight_test_data()
  model1 <- glm(
    count ~ x + z,
    family = poisson(),
    weights = w,
    data = dat
  )
  model2 <- glm(
    count ~ x + z + mediator,
    family = poisson(),
    weights = w,
    data = dat
  )
  fit <- suest(
    model1,
    model2,
    model_names = c("Base", "Adjusted"),
    weight_type = "pweight"
  )

  testthat::expect_equal(
    unname(vcov(fit)[fit$index[[1]], fit$index[[1]]]),
    unname(pweight_poisson_vcov(model1)),
    tolerance = 1e-9
  )
  testthat::expect_equal(
    unname(vcov(fit)[fit$index[[2]], fit$index[[2]]]),
    unname(pweight_poisson_vcov(model2)),
    tolerance = 1e-9
  )
  testthat::expect_equal(
    unname(vcov(fit)[fit$index[[1]], fit$index[[2]]]),
    unname(pweight_cross_vcov(
      model1,
      model2,
      pweight_poisson_component,
      pweight_poisson_component
    )),
    tolerance = 1e-9
  )
})


testthat::test_that("pweights support Poisson overlap and disjoint samples", {
  dat <- pweight_test_data()
  dat$w_b <- dat$w
  dat$w_b[dat$id > 300] <- dat$w[dat$id > 300] *
    (1.25 + 0.1 * abs(dat$x[dat$id > 300]))
  dat$w_left <- exp(0.2 * dat$x - 0.1 * dat$z)
  dat$w_right <- exp(-0.15 * dat$x + 0.3 * dat$z)

  partial1 <- glm(
    count ~ x + z,
    family = poisson(),
    weights = w,
    data = dat,
    subset = id <= 300
  )
  partial2 <- glm(
    count ~ x + z + mediator,
    family = poisson(),
    weights = w_b,
    data = dat,
    subset = id >= 101
  )
  partial_fit <- suest(partial1, partial2, weight_type = "pweight")
  testthat::expect_equal(partial_fit$nobs_overlap, 200)

  disjoint1 <- glm(
    count ~ x + z,
    family = poisson(),
    weights = w_left,
    data = dat,
    subset = id <= 200
  )
  disjoint2 <- glm(
    count ~ x + z + mediator,
    family = poisson(),
    weights = w_right,
    data = dat,
    subset = id > 200
  )
  disjoint_fit <- suest(disjoint1, disjoint2, weight_type = "pweight")
  testthat::expect_equal(disjoint_fit$nobs_overlap, 0)
  testthat::expect_true(all(
    unname(vcov(disjoint_fit)[disjoint_fit$index[[1]],
                               disjoint_fit$index[[2]]]) == 0
  ))
})


testthat::test_that("glm2 pweighted Poisson models match stats glm", {
  testthat::skip_if_not_installed("glm2")
  dat <- pweight_test_data()
  stats1 <- glm(count ~ x + z, family = poisson(), weights = w, data = dat)
  stats2 <- glm(
    count ~ x + z + mediator,
    family = poisson(),
    weights = w,
    data = dat
  )
  glm21 <- glm2::glm2(
    count ~ x + z,
    family = poisson(),
    weights = w,
    data = dat
  )
  glm22 <- glm2::glm2(
    count ~ x + z + mediator,
    family = poisson(),
    weights = w,
    data = dat
  )
  stats_fit <- suest(
    stats1,
    stats2,
    model_names = c("Base", "Adjusted"),
    weight_type = "pweight"
  )
  glm2_fit <- suest(
    glm21,
    glm22,
    model_names = c("Base", "Adjusted"),
    weight_type = "pweight"
  )

  testthat::expect_equal(coef(glm2_fit), coef(stats_fit), tolerance = 1e-8)
  testthat::expect_equal(vcov(glm2_fit), vcov(stats_fit), tolerance = 1e-8)
})


testthat::test_that("pweight Poisson results are invariant to common scaling", {
  dat <- pweight_test_data()
  model1 <- glm(count ~ x + z, family = poisson(), weights = w, data = dat)
  model2 <- glm(
    count ~ x + z + mediator,
    family = poisson(),
    weights = w,
    data = dat
  )

  scale_stored_weights <- function(model, multiplier) {
    scaled <- model
    scaled$prior.weights <- multiplier * model$prior.weights
    scaled$model[["(weights)"]] <- multiplier * model$model[["(weights)"]]
    scaled
  }

  scaled1 <- scale_stored_weights(model1, 10)
  scaled2 <- scale_stored_weights(model2, 10)
  fit <- suest(
    model1,
    model2,
    model_names = c("Base", "Adjusted"),
    weight_type = "pweight"
  )
  scaled_fit <- suest(
    scaled1,
    scaled2,
    model_names = c("Base", "Adjusted"),
    weight_type = "pweight"
  )

  testthat::expect_equal(coef(scaled_fit), coef(fit), tolerance = 1e-12)
  testthat::expect_equal(vcov(scaled_fit), vcov(fit), tolerance = 1e-12)
})


testthat::test_that("suest_newdata exposes model-specific pweights", {
  dat <- pweight_test_data()
  model1 <- lm(y ~ x, weights = w, data = dat, subset = id <= 300)
  model2 <- lm(
    y ~ x + z,
    weights = w,
    data = dat,
    subset = id >= 101
  )
  fit <- suest(
    model1,
    model2,
    model_names = c("First", "Second"),
    weight_type = "pweight"
  )
  nd <- suest_newdata(fit)

  testthat::expect_true(".suest_weight" %in% names(nd))
  testthat::expect_equal(
    nd$.suest_weight[nd$.suest_model == "First"],
    stats::model.weights(stats::model.frame(model1))
  )
  testthat::expect_equal(
    nd$.suest_weight[nd$.suest_model == "Second"],
    stats::model.weights(stats::model.frame(model2))
  )
})


testthat::test_that("marginaleffects can average pweighted lm comparisons", {
  dat <- pweight_test_data()
  model1 <- lm(y ~ x + I(x^2) + z, weights = w, data = dat)
  model2 <- lm(
    y ~ x + I(x^2) + z + mediator,
    weights = w,
    data = dat
  )
  fit <- suest(
    model1,
    model2,
    model_names = c("Base", "Adjusted"),
    weight_type = "pweight"
  )

  effects <- marginaleffects::avg_comparisons(
    fit,
    variables = list(x = 1),
    newdata = dat,
    wts = "w"
  )
  effects <- effects[
    match(c("Base", "Adjusted"), as.character(effects$group)),
    ,
    drop = FALSE
  ]

  change <- 2 * dat$x + 1
  expected <- c(
    stats::coef(model1)["x"] +
      stats::coef(model1)["I(x^2)"] *
        stats::weighted.mean(change, dat$w),
    stats::coef(model2)["x"] +
      stats::coef(model2)["I(x^2)"] *
        stats::weighted.mean(change, dat$w)
  )

  testthat::expect_equal(effects$estimate, unname(expected), tolerance = 1e-8)
  testthat::expect_true(all(is.finite(effects$std.error)))

  nd <- suest_newdata(fit)
  effects_nd <- marginaleffects::avg_comparisons(
    fit,
    variables = list(x = 1),
    newdata = nd,
    wts = ".suest_weight"
  )
  effects_nd <- effects_nd[
    match(c("Base", "Adjusted"), as.character(effects_nd$group)),
    ,
    drop = FALSE
  ]
  testthat::expect_equal(effects_nd$estimate, effects$estimate, tolerance = 1e-8)
  testthat::expect_equal(effects_nd$std.error, effects$std.error, tolerance = 1e-6)
})


testthat::test_that("marginaleffects can average pweighted binary comparisons", {
  dat <- pweight_test_data()
  model1 <- suppressWarnings(glm(
    binary ~ x + z,
    family = binomial(),
    weights = w,
    data = dat
  ))
  model2 <- suppressWarnings(glm(
    binary ~ x + z + mediator,
    family = binomial(),
    weights = w,
    data = dat
  ))
  fit <- suest(
    model1,
    model2,
    model_names = c("Base", "Adjusted"),
    weight_type = "pweight"
  )
  effects <- marginaleffects::avg_comparisons(
    fit,
    variables = "x",
    newdata = dat,
    wts = "w"
  )
  nd <- suest_newdata(fit)
  effects_nd <- marginaleffects::avg_comparisons(
    fit,
    variables = "x",
    newdata = nd,
    wts = ".suest_weight"
  )
  effects <- effects[match(c("Base", "Adjusted"), effects$group), ]
  effects_nd <- effects_nd[match(c("Base", "Adjusted"), effects_nd$group), ]

  testthat::expect_equal(
    effects_nd$estimate,
    effects$estimate,
    tolerance = 1e-8
  )
  testthat::expect_equal(
    effects_nd$std.error,
    effects$std.error,
    tolerance = 1e-6
  )
})


testthat::test_that("marginaleffects can average pweighted Poisson comparisons", {
  dat <- pweight_test_data()
  model1 <- glm(count ~ x + z, family = poisson(), weights = w, data = dat)
  model2 <- glm(
    count ~ x + z + mediator,
    family = poisson(),
    weights = w,
    data = dat
  )
  fit <- suest(
    model1,
    model2,
    model_names = c("Base", "Adjusted"),
    weight_type = "pweight"
  )
  effects <- marginaleffects::avg_comparisons(
    fit,
    variables = "x",
    newdata = dat,
    wts = "w"
  )
  nd <- suest_newdata(fit)
  effects_nd <- marginaleffects::avg_comparisons(
    fit,
    variables = "x",
    newdata = nd,
    wts = ".suest_weight"
  )
  effects <- effects[match(c("Base", "Adjusted"), effects$group), ]
  effects_nd <- effects_nd[match(c("Base", "Adjusted"), effects_nd$group), ]

  testthat::expect_equal(effects_nd$estimate, effects$estimate, tolerance = 1e-8)
  testthat::expect_equal(effects_nd$std.error, effects$std.error, tolerance = 1e-6)
})

testthat::test_that("pweighted negative-binomial results are scale invariant", {
  dat <- pweight_test_data(800)
  mu <- exp(0.3 + 0.3*dat$x - 0.2*dat$z + 0.2*dat$mediator)
  dat$count_nb <- rnbinom(nrow(dat), mu = mu, size = 2)
  control <- glm.control(maxit = 100)
  model1 <- MASS::glm.nb(
    count_nb ~ x + z,
    weights = w,
    data = dat,
    control = control
  )
  model2 <- MASS::glm.nb(
    count_nb ~ x + z + mediator,
    weights = w,
    data = dat,
    control = control
  )
  scaled1 <- model1
  scaled2 <- model2
  scaled1$model[["(weights)"]] <- 10*model1$model[["(weights)"]]
  scaled2$model[["(weights)"]] <- 10*model2$model[["(weights)"]]
  fit <- suest(
    model1,
    model2,
    model_names = c("Base", "Adjusted"),
    weight_type = "pweight"
  )
  scaled <- suest(
    scaled1,
    scaled2,
    model_names = c("Base", "Adjusted"),
    weight_type = "pweight"
  )

  testthat::expect_equal(coef(scaled), coef(fit), tolerance = 1e-12)
  testthat::expect_equal(vcov(scaled), vcov(fit), tolerance = 1e-12)
})

testthat::test_that("pweighted categorical systems work with marginaleffects", {
  dat <- pweight_test_data(800)
  eta <- 0.5*dat$x - 0.25*dat$z + 0.3*dat$mediator + rlogis(nrow(dat))
  dat$ordinal <- ordered(
    cut(eta, c(-Inf, -0.8, 0, 0.8, Inf), labels = FALSE),
    levels = 1:4,
    labels = c("A", "B", "C", "D")
  )
  dat$nominal <- factor(dat$ordinal, levels = levels(dat$ordinal))
  ologit <- suppressWarnings(MASS::polr(
    ordinal ~ x + z,
    data = dat,
    weights = w,
    method = "logistic",
    Hess = TRUE
  ))
  oprobit <- suppressWarnings(MASS::polr(
    ordinal ~ x + z,
    data = dat,
    weights = w,
    method = "probit",
    Hess = TRUE
  ))
  multinom <- nnet::multinom(
    nominal ~ x + z,
    data = dat,
    weights = w,
    Hess = TRUE,
    trace = FALSE
  )
  fit <- suest(
    ologit,
    oprobit,
    multinom,
    model_names = c("Ologit", "Oprobit", "Multinom"),
    weight_type = "pweight"
  )
  effects <- marginaleffects::avg_comparisons(
    fit,
    variables = "x",
    newdata = suest_newdata(fit),
    wts = ".suest_weight"
  )

  testthat::expect_equal(nrow(effects), 12L)
  testthat::expect_true(all(is.finite(effects$estimate)))
  testthat::expect_true(all(is.finite(effects$std.error)))
})

testthat::test_that("pweighted categorical results are scale invariant", {
  dat <- pweight_test_data(800)
  eta <- 0.4*dat$x - 0.2*dat$z + rlogis(nrow(dat))
  dat$ordinal <- ordered(
    cut(eta, c(-Inf, -0.7, 0, 0.7, Inf), labels = FALSE),
    levels = 1:4,
    labels = c("A", "B", "C", "D")
  )
  dat$nominal <- factor(dat$ordinal, levels = levels(dat$ordinal))
  models <- list(
    suppressWarnings(MASS::polr(
      ordinal ~ x + z,
      data = dat,
      weights = w,
      method = "logistic",
      Hess = TRUE
    )),
    nnet::multinom(
      nominal ~ x + z,
      data = dat,
      weights = w,
      Hess = TRUE,
      trace = FALSE
    )
  )
  scaled_models <- models
  scaled_models[[1L]]$model[["(weights)"]] <-
    10*models[[1L]]$model[["(weights)"]]
  scaled_models[[2L]]$weights <- 10*models[[2L]]$weights
  fit <- suest(
    models[[1L]],
    models[[2L]],
    model_names = c("Ordered", "Multinomial"),
    weight_type = "pweight"
  )
  scaled <- suest(
    scaled_models[[1L]],
    scaled_models[[2L]],
    model_names = c("Ordered", "Multinomial"),
    weight_type = "pweight"
  )

  testthat::expect_equal(coef(scaled), coef(fit), tolerance = 1e-12)
  testthat::expect_equal(vcov(scaled), vcov(fit), tolerance = 1e-12)
})
