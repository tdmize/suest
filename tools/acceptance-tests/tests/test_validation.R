cat("\n\nVALIDATION AND FAILURE TESTS\n")

test_case("Validation: Gaussian glm is supported", {
  data <- datasets::mtcars
  g1 <- stats::glm(mpg ~ wt, family = stats::gaussian(), data = data)
  g2 <- stats::glm(mpg ~ wt + hp, family = stats::gaussian(), data = data)
  combined <- suest(g1, g2)
  expect_true(all(is.finite(stats::coef(combined))))
  expect_true(all(is.finite(stats::vcov(combined))))
  expect_true(combined$nobs_overlap == nrow(data))
  combined
})

test_case("Validation: scalar-categorical pair is supported", {
  data <- datasets::mtcars
  data$am <- factor(data$am)
  data$cyl_ord <- ordered(data$cyl)
  logit <- stats::glm(
    am ~ wt,
    family = stats::binomial(),
    data = data
  )
  ordered <- MASS::polr(
    cyl_ord ~ wt,
    data = data,
    method = "logistic",
    Hess = TRUE
  )
  combined <- suest(logit, ordered, model_names = c("Logit", "Ordered"))
  effects <- marginaleffects::avg_comparisons(
    combined,
    variables = "wt",
    newdata = data
  )
  expect_true(nrow(effects) == 4L)
  expect_true(all(is.finite(effects$std.error)))
  effects
})

test_case("Validation: ordered logit-ordered probit is supported", {
  data <- MASS::housing[
    rep(seq_len(nrow(MASS::housing)), MASS::housing$Freq),
    c("Sat", "Infl", "Type", "Cont")
  ]
  data$Sat <- ordered(data$Sat, levels = c("Low", "Medium", "High"))
  ol <- MASS::polr(
    Sat ~ Cont,
    data = data,
    method = "logistic",
    Hess = TRUE
  )
  op <- MASS::polr(
    Sat ~ Cont,
    data = data,
    method = "probit",
    Hess = TRUE
  )
  combined <- suest(ol, op)
  expect_true(all(is.finite(stats::coef(combined))))
  expect_true(all(is.finite(stats::vcov(combined))))
  expect_true(combined$nobs_overlap == nrow(data))
})

test_case("Validation: different categorical outcome levels are supported", {
  data <- datasets::iris
  data$Species2 <- factor(
    data$Species,
    levels = rev(levels(data$Species))
  )
  m1 <- nnet::multinom(
    Species ~ Sepal.Length,
    data = data,
    Hess = TRUE,
    trace = FALSE
  )
  m2 <- nnet::multinom(
    Species2 ~ Sepal.Length,
    data = data,
    Hess = TRUE,
    trace = FALSE
  )
  combined <- suest(m1, m2, model_names = c("Forward", "Reverse"))
  predictions <- marginaleffects::avg_predictions(combined, newdata = data)
  expect_true(nrow(predictions) == 6L)
  expect_true(all(is.finite(predictions$std.error)))
  predictions
})

test_case("Validation: nonunit weights", {
  data <- datasets::mtcars
  m1 <- stats::lm(mpg ~ wt, data = data, weights = rep(2, nrow(data)))
  m2 <- stats::lm(mpg ~ wt + hp, data = data)
  expect_error(
    suest(m1, m2),
    "Weights"
  )
})

test_case("Validation: offsets are supported", {
  data <- datasets::warpbreaks
  p1 <- stats::glm(
    breaks ~ wool + offset(rep(0.5, nrow(data))),
    family = stats::poisson(),
    data = data
  )
  p2 <- stats::glm(
    breaks ~ wool + tension,
    family = stats::poisson(),
    data = data
  )
  combined <- suest(p1, p2)
  expect_true(all(is.finite(stats::coef(combined))))
  expect_true(all(is.finite(stats::vcov(combined))))
  expect_true(combined$nobs_overlap == nrow(data))
  combined
})

test_case("Validation: aliased parameters", {
  data <- datasets::mtcars
  data$wt_copy <- data$wt
  m1 <- stats::lm(mpg ~ wt + wt_copy, data = data)
  m2 <- stats::lm(mpg ~ wt + hp, data = data)
  expect_error(
    suest(m1, m2),
    "aliased|missing parameters"
  )
})

test_case("Validation: cross-type predictions are response scale only", {
  data <- datasets::mtcars
  data$am <- factor(data$am)
  logit <- stats::glm(
    am ~ wt,
    family = stats::binomial("logit"),
    data = data
  )
  probit <- stats::glm(
    am ~ wt,
    family = stats::binomial("probit"),
    data = data
  )
  object <- suest(logit, probit)
  expect_error(
    marginaleffects::predictions(
      object,
      newdata = data,
      type = "link"
    ),
    "response scale"
  )
})
