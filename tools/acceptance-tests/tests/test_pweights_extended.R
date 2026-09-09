cat("\n\nPWEIGHT NEGATIVE-BINOMIAL AND CATEGORICAL TESTS\n")

test_case("Pweights: negative-binomial scale invariance", {
  pwdat$y_nb <- pwdat$y_count + 3L*(pwdat$id %% 11L == 0L) +
    5L*(pwdat$id %% 37L == 0L)
  control <- stats::glm.control(maxit = 100)
  model1 <- MASS::glm.nb(
    y_nb ~ x1 + x2 + z,
    weights = pw,
    data = pwdat,
    control = control
  )
  model2 <- MASS::glm.nb(
    y_nb ~ x1 + x2 + z + mediator,
    weights = pw,
    data = pwdat,
    control = control
  )
  scaled1 <- model1
  scaled2 <- model2
  scaled1$model[["(weights)"]] <- 10*model1$model[["(weights)"]]
  scaled2$model[["(weights)"]] <- 10*model2$model[["(weights)"]]
  combined <- suest(
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

  expected <- c(
    `Base::(Intercept)` = 0.654173840317753,
    `Base::x1` = 0.358477396490002,
    `Base::x2` = -0.221957616055218,
    `Base::z1` = 0.319138291940327,
    `Base::ln_theta` = 1.39232124289059,
    `Adjusted::(Intercept)` = 0.626377507464159,
    `Adjusted::x1` = 0.2553295308383,
    `Adjusted::x2` = -0.152706719332316,
    `Adjusted::z1` = 0.26207790366685,
    `Adjusted::mediator` = 0.202172869964593,
    `Adjusted::ln_theta` = 1.63239377715187
  )

  expect_near(
    stats::coef(combined)[names(expected)],
    expected,
    tolerance = 1e-7,
    label = "Stata suest2 pweight negative-binomial coefficients"
  )

  V <- stats::vcov(combined)
  base <- combined$index[[1L]]
  adjusted <- combined$index[[2L]]
  names(base) <- combined$local_names[[1L]]
  names(adjusted) <- combined$local_names[[2L]]

  # Stata reports ln(alpha), whereas MASS and the R adapter report
  # ln(theta) = -ln(alpha). Signs below are transformed to the R scale.
  expect_near(
    c(
      V[base["x1"], base["x1"]],
      V[base["ln_theta"], base["ln_theta"]],
      V[adjusted["x1"], adjusted["x1"]],
      V[adjusted["ln_theta"], adjusted["ln_theta"]],
      V[base["x1"], adjusted["x1"]],
      V[base["ln_theta"], adjusted["ln_theta"]],
      V[base["x1"], base["ln_theta"]],
      V[base["ln_theta"], adjusted["x1"]],
      V[adjusted["mediator"], adjusted["ln_theta"]]
    ),
    c(
      0.000389715612719831,
      0.00800523809287427,
      0.000434346191920883,
      0.0123448886056109,
      0.000361911085284041,
      0.009196175372499,
      0.000396039333490631,
      0.000557374755216374,
      0.000366908554008439
    ),
    tolerance = 2e-9,
    label = "Stata suest2 pweight negative-binomial covariance"
  )

  expect_near(
    stats::coef(scaled),
    stats::coef(combined),
    tolerance = 1e-12,
    label = "negative-binomial pweight scaled coefficients"
  )
  expect_near(
    stats::vcov(scaled),
    stats::vcov(combined),
    tolerance = 1e-12,
    label = "negative-binomial pweight scaled covariance"
  )

  list(object = combined, theta = combined$theta)
})

test_case("Pweights: three-model categorical system", {
  ordered_rank <- rank(pwdat$y_linear, ties.method = "first")
  category <- cut(
    ordered_rank,
    breaks = c(0, 600, 1200, 1800, 2400),
    labels = c("A", "B", "C", "D")
  )
  pwdat$y_ord <- ordered(category, levels = c("A", "B", "C", "D"))
  pwdat$y_nom <- factor(pwdat$y_ord, levels = levels(pwdat$y_ord))
  ologit <- suppressWarnings(MASS::polr(
    y_ord ~ x1 + x2 + z,
    weights = pw,
    data = pwdat,
    method = "logistic",
    Hess = TRUE
  ))
  oprobit <- suppressWarnings(MASS::polr(
    y_ord ~ x1 + x2 + z,
    weights = pw,
    data = pwdat,
    method = "probit",
    Hess = TRUE
  ))
  multinom <- nnet::multinom(
    y_nom ~ x1 + x2 + z,
    weights = pw,
    data = pwdat,
    Hess = TRUE,
    trace = FALSE
  )
  combined <- suest(
    ologit,
    oprobit,
    multinom,
    model_names = c("Ologit", "Oprobit", "Multinom"),
    weight_type = "pweight"
  )

  expected <- c(
    `Ologit::x1` = 1.54866398033679,
    `Ologit::x2` = -0.88515202885839,
    `Ologit::z1` = 1.11708493241771,
    `Ologit::A|B` = -1.14326435793516,
    `Ologit::B|C` = 0.609018181943836,
    `Ologit::C|D` = 2.29600482060063,
    `Oprobit::x1` = 0.869518710723655,
    `Oprobit::x2` = -0.514810473193384,
    `Oprobit::z1` = 0.644857447524903,
    `Oprobit::A|B` = -0.640679630178486,
    `Oprobit::B|C` = 0.368541658219638,
    `Oprobit::C|D` = 1.34129044306089,
    `Multinom::B:(Intercept)` = 0.466098718919969,
    `Multinom::B:x1` = 1.31010916218924,
    `Multinom::B:x2` = -0.603582977394302,
    `Multinom::B:z1` = 0.724610977371939,
    `Multinom::C:(Intercept)` = 0.133242216943626,
    `Multinom::C:x1` = 2.08956505911389,
    `Multinom::C:x2` = -1.1164430671209,
    `Multinom::C:z1` = 1.29362447308339,
    `Multinom::D:(Intercept)` = -0.955177723675979,
    `Multinom::D:x1` = 2.90112932317841,
    `Multinom::D:x2` = -1.70012016832426,
    `Multinom::D:z1` = 2.15138728587451
  )

  expect_near(
    stats::coef(combined)[names(expected)],
    expected,
    tolerance = 6e-5,
    label = "Stata suest2 pweight categorical coefficients"
  )

  V <- stats::vcov(combined)
  ordered_logit <- combined$index[[1L]]
  ordered_probit <- combined$index[[2L]]
  multinomial <- combined$index[[3L]]
  names(ordered_logit) <- combined$local_names[[1L]]
  names(ordered_probit) <- combined$local_names[[2L]]
  names(multinomial) <- combined$local_names[[3L]]

  expect_near(
    c(
      V[ordered_logit["x1"], ordered_logit["x1"]],
      V[ordered_logit["A|B"], ordered_logit["A|B"]],
      V[ordered_probit["x1"], ordered_probit["x1"]],
      V[ordered_probit["C|D"], ordered_probit["C|D"]],
      V[multinomial["B:x1"], multinomial["B:x1"]],
      V[multinomial["C:x1"], multinomial["C:x1"]],
      V[multinomial["D:z1"], multinomial["D:z1"]],
      V[ordered_logit["x1"], ordered_probit["x1"]],
      V[ordered_logit["x1"], multinomial["B:x1"]],
      V[ordered_logit["B|C"], multinomial["C:z1"]],
      V[ordered_probit["B|C"], multinomial["B:x1"]],
      V[multinomial["B:x1"], multinomial["C:x1"]],
      V[multinomial["C:x1"], multinomial["D:x1"]]
    ),
    c(
      0.00382824980178405,
      0.00536480349083837,
      0.00122846895986791,
      0.00215472320874687,
      0.0121049707735076,
      0.0168857847756299,
      0.033349068041068,
      0.00210068065862767,
      0.00389361830522395,
      0.00535188615150891,
      0.000639662105433944,
      0.0110670732868086,
      0.0166866135427096
    ),
    tolerance = 4e-7,
    label = "Stata suest2 pweight categorical covariance"
  )

  effects <- marginaleffects::avg_comparisons(
    combined,
    variables = "x1",
    newdata = suest_newdata(combined),
    wts = ".suest_weight"
  )

  expect_true(nrow(effects) == 12L)
  expect_true(all(is.finite(effects$estimate)))
  expect_true(all(is.finite(effects$std.error)))
  expect_true(max(abs(offdiag_vcov(combined))) > 0)

  list(object = combined, effects = effects)
})
