cat("\n\nVALIDATION AND FAILURE TESTS\n")

test_case("Validation: unsupported Gaussian glm", {
  data <- datasets::mtcars
  g1 <- stats::glm(mpg ~ wt, family = stats::gaussian(), data = data)
  g2 <- stats::glm(mpg ~ wt + hp, family = stats::gaussian(), data = data)
  expect_error(
    suest(g1, g2),
    "supported"
  )
})

test_case("Validation: unsupported binary logit-Poisson pair", {
  data <- datasets::mtcars
  data$am <- factor(data$am)
  logit <- stats::glm(
    am ~ wt,
    family = stats::binomial(),
    data = data
  )
  poisson <- stats::glm(
    cyl ~ wt,
    family = stats::poisson(),
    data = data
  )
  expect_error(
    suest(logit, poisson),
    "not supported"
  )
})

test_case("Validation: ordered logit-ordered probit is rejected", {
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
  expect_error(
    suest(ol, op),
    "not supported"
  )
})

test_case("Validation: incompatible categorical outcome levels", {
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
  expect_error(
    suest(m1, m2),
    "same outcome categories"
  )
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

test_case("Validation: offsets", {
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
  expect_error(
    suest(p1, p2),
    "Offsets"
  )
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
