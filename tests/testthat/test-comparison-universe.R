comparison_universe_data <- function(n = 500L) {
  set.seed(6109)
  x <- stats::rnorm(n)
  z <- stats::rnorm(n)
  mediator <- 0.6 * x + 0.2 * z + stats::rnorm(n)
  eta <- -0.2 + 0.4 * x + 0.25 * mediator - 0.15 * z

  latent <- eta + stats::rlogis(n)
  ordered_outcome <- ordered(
    cut(
      latent,
      breaks = c(-Inf, -0.8, 0, 0.8, Inf),
      labels = FALSE
    ),
    levels = 1:4,
    labels = c("A", "B", "C", "D")
  )

  nominal_outcome <- factor(
    ordered_outcome,
    levels = levels(ordered_outcome)
  )

  data.frame(
    x = x,
    z = z,
    mediator = mediator,
    woman = rep(c(0, 1), length.out = n),
    y_lm = eta + stats::rnorm(n),
    y_bin = stats::rbinom(n, 1, stats::plogis(eta)),
    y_pois = stats::rpois(n, exp(0.3 + 0.15 * eta)),
    y_nb = stats::rnbinom(n, mu = exp(0.3 + 0.15 * eta), size = 2),
    y_ord = ordered_outcome,
    y_nom = nominal_outcome
  )
}

test_that("nested comparisons work for every supported family", {
  dat <- comparison_universe_data()

  models <- list(
    lm = list(
      stats::lm(y_lm ~ x + z, data = dat),
      stats::lm(y_lm ~ x + mediator + z, data = dat)
    ),
    logit = list(
      stats::glm(y_bin ~ x + z, family = binomial(), data = dat),
      stats::glm(
        y_bin ~ x + mediator + z,
        family = binomial(),
        data = dat
      )
    ),
    probit = list(
      stats::glm(
        y_bin ~ x + z,
        family = binomial("probit"),
        data = dat
      ),
      stats::glm(
        y_bin ~ x + mediator + z,
        family = binomial("probit"),
        data = dat
      )
    ),
    poisson = list(
      stats::glm(
        y_pois ~ x + z,
        family = poisson(),
        data = dat
      ),
      stats::glm(
        y_pois ~ x + mediator + z,
        family = poisson(),
        data = dat
      )
    ),
    negbin = list(
      MASS::glm.nb(y_nb ~ x + z, data = dat),
      MASS::glm.nb(y_nb ~ x + mediator + z, data = dat)
    ),
    ologit = list(
      MASS::polr(
        y_ord ~ x + z,
        data = dat,
        method = "logistic",
        Hess = TRUE
      ),
      MASS::polr(
        y_ord ~ x + mediator + z,
        data = dat,
        method = "logistic",
        Hess = TRUE
      )
    ),
    oprobit = list(
      MASS::polr(
        y_ord ~ x + z,
        data = dat,
        method = "probit",
        Hess = TRUE
      ),
      MASS::polr(
        y_ord ~ x + mediator + z,
        data = dat,
        method = "probit",
        Hess = TRUE
      )
    ),
    multinom = list(
      nnet::multinom(
        y_nom ~ x + z,
        data = dat,
        Hess = TRUE,
        trace = FALSE
      ),
      nnet::multinom(
        y_nom ~ x + mediator + z,
        data = dat,
        Hess = TRUE,
        trace = FALSE
      )
    )
  )

  for (family_name in names(models)) {
    fit <- suest(
      models[[family_name]][[1L]],
      models[[family_name]][[2L]],
      model_names = c("Base", "Mediator")
    )

    effects <- marginaleffects::avg_comparisons(
      fit,
      variables = "x",
      newdata = dat
    )

    expect_true(all(is.finite(effects$estimate)))
    expect_true(all(is.finite(effects$std.error)))
    expect_gt(max(abs(fit$vcov[
      fit$index[[1L]],
      fit$index[[2L]],
      drop = FALSE
    ])), 0)
  }
})

test_that("sex-stratified subset comparisons are disjoint", {
  dat <- comparison_universe_data()

  men <- stats::glm(
    y_bin ~ x + mediator + z,
    family = binomial(),
    data = dat,
    subset = woman == 0
  )
  women <- stats::glm(
    y_bin ~ x + mediator + z,
    family = binomial(),
    data = dat,
    subset = woman == 1
  )

  fit <- suest(men, women, model_names = c("Men", "Women"))

  expect_equal(fit$nobs_overlap, 0)
  expect_equal(
    fit$vcov[
      fit$index[[1L]],
      fit$index[[2L]],
      drop = FALSE
    ],
    matrix(
      0,
      length(fit$index[[1L]]),
      length(fit$index[[2L]])
    )
  )

  effects <- marginaleffects::avg_comparisons(
    fit,
    variables = "x",
    newdata = suest_newdata(fit)
  )

  expect_true(all(is.finite(effects$std.error)))
})
