test_that("marginaleffects computes cross-model comparisons", {
  dat <- mtcars
  dat$am <- factor(dat$am)

  model1 <- glm(am ~ wt, family = binomial(), data = dat)
  model2 <- glm(am ~ wt + hp, family = binomial(), data = dat)

  fit <- suest(
    model1,
    model2,
    model_names = c("Base", "Adjusted")
  )

  effects <- marginaleffects::avg_comparisons(
    fit,
    variables = "wt",
    newdata = dat
  )

  expect_equal(nrow(effects), 2)
  expect_true(all(is.finite(effects$estimate)))
  expect_true(all(is.finite(effects$std.error)))

  difference <- marginaleffects::hypotheses(
    effects,
    hypothesis = difference ~ revpairwise
  )

  expect_equal(nrow(difference), 1)
  expect_true(is.finite(difference$estimate))
  expect_true(is.finite(difference$std.error))
})

test_that("suest_newdata averages within separate samples", {
  set.seed(376)
  dat <- data.frame(
    group = rep(c("A", "B"), each = 200),
    x = rnorm(400),
    z = rnorm(400)
  )
  dat$y <- rbinom(
    400,
    1,
    plogis(-0.3 + 0.5 * dat$x - 0.2 * dat$z)
  )

  model1 <- glm(
    y ~ x + z,
    family = binomial(),
    data = dat,
    subset = group == "A"
  )
  model2 <- glm(
    y ~ x + z,
    family = binomial(),
    data = dat,
    subset = group == "B"
  )

  fit <- suest(model1, model2, model_names = c("A", "B"))
  nd <- suest_newdata(fit)

  effects <- marginaleffects::avg_comparisons(
    fit,
    variables = "x",
    newdata = nd
  )

  expect_equal(nrow(effects), 2)
  expect_true(all(is.finite(effects$std.error)))
})
