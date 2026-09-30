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

test_that("disjoint samples retain every multi-category factor contrast", {
  set.seed(377)
  n <- 700L
  reltrad_levels <- paste0("R", seq_len(7L))
  dat <- data.frame(
    sample = rep(c("Pre", "Post"), each = n / 2L),
    reltrad = factor(
      rep(reltrad_levels, length.out = n),
      levels = reltrad_levels
    ),
    age = stats::rnorm(n, mean = 45, sd = 15),
    woman = factor(rep(c("Men", "Women"), length.out = n))
  )
  reltrad_effect <- c(-0.45, -0.25, -0.05, 0.10, 0.25, 0.40, 0.15)
  eta <- -0.3 +
    reltrad_effect[as.integer(dat$reltrad)] +
    0.01 * (dat$age - 45) +
    0.2 * (dat$sample == "Post")
  dat$y <- stats::rbinom(n, 1, stats::plogis(eta))

  premod <- stats::glm(
    y ~ reltrad + age + woman,
    family = stats::binomial(),
    data = dat,
    subset = sample == "Pre"
  )
  postmod <- stats::glm(
    y ~ reltrad + age + woman,
    family = stats::binomial(),
    data = dat,
    subset = sample == "Post"
  )

  fit <- suest(premod, postmod, model_names = c("Pre", "Post"))
  nd <- suest_newdata(fit)

  expect_true(is.factor(nd$reltrad))
  expect_equal(levels(nd$reltrad), reltrad_levels)

  joint_pw <- marginaleffects::avg_comparisons(
    fit,
    variables = list(reltrad = "pairwise"),
    newdata = nd
  )
  joint_ref <- marginaleffects::avg_comparisons(
    fit,
    variables = list(reltrad = "reference"),
    newdata = nd
  )

  expect_equal(as.integer(table(joint_pw$group)), c(21L, 21L))
  expect_equal(as.integer(table(joint_ref$group)), c(6L, 6L))

  sep_pre <- marginaleffects::avg_comparisons(
    premod,
    variables = list(reltrad = "pairwise")
  )
  sep_post <- marginaleffects::avg_comparisons(
    postmod,
    variables = list(reltrad = "pairwise")
  )

  compare_estimates <- function(joint, separate, group) {
    j <- joint[joint$group == group, c("contrast", "estimate")]
    s <- separate[, c("contrast", "estimate")]
    idx <- match(as.character(s$contrast), as.character(j$contrast))
    expect_false(anyNA(idx))
    expect_equal(j$estimate[idx], s$estimate, tolerance = 1e-8)
  }

  compare_estimates(joint_pw, sep_pre, "Pre")
  compare_estimates(joint_pw, sep_post, "Post")

  binary <- marginaleffects::avg_comparisons(
    fit,
    variables = "woman",
    newdata = nd
  )
  binary_pre <- marginaleffects::avg_comparisons(premod, variables = "woman")
  binary_post <- marginaleffects::avg_comparisons(postmod, variables = "woman")

  expect_equal(nrow(binary), 2L)
  expect_equal(binary$estimate[binary$group == "Pre"], binary_pre$estimate,
               tolerance = 1e-8)
  expect_equal(binary$estimate[binary$group == "Post"], binary_post$estimate,
               tolerance = 1e-8)

  predictions <- marginaleffects::avg_predictions(
    fit,
    newdata = nd,
    by = "reltrad"
  )
  pred_pre <- marginaleffects::avg_predictions(premod, by = "reltrad")
  pred_post <- marginaleffects::avg_predictions(postmod, by = "reltrad")

  expect_equal(nrow(predictions), 14L)

  compare_predictions <- function(joint, separate, group) {
    j <- joint[joint$group == group, c("reltrad", "estimate")]
    idx <- match(as.character(separate$reltrad), as.character(j$reltrad))
    expect_false(anyNA(idx))
    expect_equal(j$estimate[idx], separate$estimate, tolerance = 1e-8)
  }

  compare_predictions(predictions, pred_pre, "Pre")
  compare_predictions(predictions, pred_post, "Post")
})
