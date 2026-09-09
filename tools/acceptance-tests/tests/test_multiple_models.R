cat("\n\nMULTIPLE-MODEL TESTS\n")

test_case("Multiple models: heterogeneous scalar system", {
  linear <- universe_fit("lm", "y_lm", "x + mediator + z")
  logit <- universe_fit("logit", "y_bin", "x + mediator + z")
  poisson <- universe_fit("poisson", "y_pois", "x + mediator + z")
  combined <- suest(
    linear,
    logit,
    poisson,
    model_names = c("Linear", "Logit", "Poisson")
  )
  effects <- marginaleffects::avg_comparisons(
    combined,
    variables = "x",
    newdata = universe_data
  )
  differences <- marginaleffects::hypotheses(
    effects,
    hypothesis = difference ~ pairwise
  )

  expect_true(length(combined$models) == 3L)
  expect_true(combined$nobs_overlap == nrow(universe_data))
  expect_true(nrow(effects) == 3L)
  expect_true(nrow(differences) == 3L)
  expect_true(all(is.finite(effects$std.error)))
  expect_true(all(is.finite(differences$std.error)))

  list(object = combined, effects = effects, differences = differences)
})

test_case("Multiple models: heterogeneous categorical system", {
  universe_data$y_ord_nom <- factor(
    universe_data$y_ord,
    levels = levels(universe_data$y_ord)
  )
  assign("universe_data", universe_data, envir = globalenv())
  ologit <- universe_fit("ologit", "y_ord", "x + mediator + z")
  oprobit <- universe_fit("oprobit", "y_ord", "x + mediator + z")
  multinom <- universe_fit("multinom", "y_ord_nom", "x + mediator + z")
  combined <- suest(
    ologit,
    oprobit,
    multinom,
    model_names = c("Ologit", "Oprobit", "Multinom")
  )
  effects <- marginaleffects::avg_comparisons(
    combined,
    variables = "x",
    newdata = universe_data
  )

  expect_true(nrow(effects) == 12L)
  expect_true(all(is.finite(effects$estimate)))
  expect_true(all(is.finite(effects$std.error)))
  totals <- tapply(effects$estimate, sub("::.*", "", effects$group), sum)
  expect_near(
    as.numeric(totals),
    rep(0, 3L),
    tolerance = 1e-8,
    label = "three-model category-effect sums"
  )

  list(object = combined, effects = effects)
})
