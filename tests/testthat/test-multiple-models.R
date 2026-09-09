multiple_model_data <- function(n = 600L) {
  set.seed(9403)
  x <- rnorm(n)
  z <- rnorm(n)
  eta <- -0.2 + 0.5*x - 0.25*z
  latent <- eta + rlogis(n)
  y_ord <- ordered(
    cut(latent, c(-Inf, -0.8, 0, 0.8, Inf), labels = FALSE),
    levels = 1:4,
    labels = c("A", "B", "C", "D")
  )
  data.frame(
    id = seq_len(n),
    x = x,
    z = z,
    y_lm = eta + rnorm(n),
    y_bin = rbinom(n, 1, plogis(eta)),
    y_pois = rpois(n, exp(0.4 + 0.2*eta)),
    y_ord = y_ord,
    y_nom = factor(y_ord, levels = levels(y_ord))
  )
}

test_that("three heterogeneous scalar models work with marginaleffects", {
  dat <- multiple_model_data()
  linear <- lm(y_lm ~ x + z, data = dat)
  logit <- glm(y_bin ~ x + z, family = binomial(), data = dat)
  poisson <- glm(y_pois ~ x + z, family = poisson(), data = dat)
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
  differences <- marginaleffects::hypotheses(
    effects,
    hypothesis = difference ~ pairwise
  )

  expect_equal(length(fit$models), 3L)
  expect_equal(fit$nobs_union, nrow(dat))
  expect_equal(fit$nobs_overlap, nrow(dat))
  expect_equal(fit$nobs_overlap_pairwise, matrix(
    nrow(dat),
    3L,
    3L,
    dimnames = list(fit$model_names, fit$model_names)
  ))
  expect_equal(nrow(effects), 3L)
  expect_equal(nrow(differences), 3L)
  expect_true(all(is.finite(effects$std.error)))
  expect_true(all(is.finite(differences$std.error)))
})

test_that("three categorical models work with marginaleffects", {
  dat <- multiple_model_data()
  ologit <- MASS::polr(
    y_ord ~ x + z,
    data = dat,
    method = "logistic",
    Hess = TRUE
  )
  oprobit <- MASS::polr(
    y_ord ~ x + z,
    data = dat,
    method = "probit",
    Hess = TRUE
  )
  multinom <- nnet::multinom(
    y_nom ~ x + z,
    data = dat,
    Hess = TRUE,
    trace = FALSE
  )
  fit <- suest(
    ologit,
    oprobit,
    multinom,
    model_names = c("Ologit", "Oprobit", "Multinom")
  )
  effects <- marginaleffects::avg_comparisons(
    fit,
    variables = "x",
    newdata = dat
  )

  expect_equal(nrow(effects), 12L)
  expect_true(all(is.finite(effects$estimate)))
  expect_true(all(is.finite(effects$std.error)))
  totals <- tapply(effects$estimate, sub("::.*", "", effects$group), sum)
  expect_equal(as.numeric(totals), rep(0, 3L), tolerance = 1e-10)
})

test_that("pairwise overlap is tracked for three separate data objects", {
  dat <- multiple_model_data()
  data1 <- dat[dat$id <= 400, ]
  data2 <- dat[dat$id >= 201, ]
  data3 <- dat[dat$id >= 101 & dat$id <= 500, ]
  model1 <- lm(y_lm ~ x, data = data1)
  model2 <- lm(y_lm ~ x + z, data = data2)
  model3 <- lm(y_lm ~ z, data = data3)
  fit <- suest(
    model1,
    model2,
    model3,
    model_names = c("A", "B", "C"),
    observation_id = "id"
  )

  expected <- matrix(
    c(400, 200, 300, 200, 400, 300, 300, 300, 400),
    3L,
    dimnames = list(fit$model_names, fit$model_names)
  )
  expect_equal(fit$nobs_union, 600L)
  expect_equal(fit$nobs_overlap, 200L)
  expect_equal(fit$nobs_overlap_pairwise, expected)
})

test_that("model names must cover every model", {
  dat <- multiple_model_data(200)
  model1 <- lm(y_lm ~ x, data = dat)
  model2 <- lm(y_lm ~ z, data = dat)
  model3 <- lm(y_lm ~ x + z, data = dat)

  expect_error(
    suest(model1, model2, model3, model_names = c("A", "B")),
    "one nonempty name per model"
  )
})

test_that("mixed scalar and categorical systems work with marginaleffects", {
  dat <- multiple_model_data()
  linear <- lm(y_lm ~ x + z, data = dat)
  ordered <- MASS::polr(
    y_ord ~ x + z,
    data = dat,
    method = "logistic",
    Hess = TRUE
  )
  nominal <- nnet::multinom(
    y_nom ~ x + z,
    data = dat,
    Hess = TRUE,
    trace = FALSE
  )
  fit <- suest(
    linear,
    ordered,
    nominal,
    model_names = c("Linear", "Ordered", "Nominal")
  )
  effects <- marginaleffects::avg_comparisons(
    fit,
    variables = "x",
    newdata = dat
  )

  expect_equal(nrow(effects), 9L)
  expect_true(all(is.finite(effects$estimate)))
  expect_true(all(is.finite(effects$std.error)))
})

test_that("three-model pweight systems validate every overlap", {
  dat <- multiple_model_data()
  dat$w <- exp(0.15*dat$x)
  dat$w_bad <- dat$w
  dat$w_bad[1L] <- 2*dat$w_bad[1L]
  linear <- lm(y_lm ~ x + z, data = dat, weights = w)
  logit <- suppressWarnings(glm(
    y_bin ~ x + z,
    family = binomial(),
    data = dat,
    weights = w
  ))
  poisson <- glm(
    y_pois ~ x + z,
    family = poisson(),
    data = dat,
    weights = w
  )
  poisson_bad <- glm(
    y_pois ~ x + z,
    family = poisson(),
    data = dat,
    weights = w_bad
  )
  fit <- suest(
    linear,
    logit,
    poisson,
    model_names = c("Linear", "Logit", "Poisson"),
    weight_type = "pweight"
  )
  effects <- marginaleffects::avg_comparisons(
    fit,
    variables = "x",
    newdata = suest_newdata(fit),
    wts = ".suest_weight"
  )

  expect_equal(nrow(effects), 3L)
  expect_true(all(is.finite(effects$std.error)))
  expect_error(
    suest(
      linear,
      logit,
      poisson_bad,
      model_names = c("Linear", "Logit", "Poisson"),
      weight_type = "pweight"
    ),
    "models 'Linear' and 'Poisson'"
  )
})
