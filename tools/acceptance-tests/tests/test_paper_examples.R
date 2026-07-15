cat("\n\nPAPER REPLICATION TESTS\n")

############################################################
# Example 6.1: Polynomial effect and mediation in lm()
############################################################

test_case("Paper 6.1: load and prepare Add Health data", {
  ah <- load_stata_data("ah4_cme.dta")
  vars61 <- c(
    "depsympB", "income", "inc10", "age",
    "woman", "race", "college", "jobsat"
  )
  ah61 <- complete_data(ah, vars61)
  ah61 <- factorize(ah61, c("woman", "race", "jobsat"))
  expect_true(nrow(ah61) == 4307, "Example 6.1 should have N=4307.")
  nrow(ah61)
})

test_case("Paper 6.1: fit and combine nested linear models", {
  lm61_base <- stats::lm(
    depsympB ~ income + I(income^2) + age + woman + race,
    data = ah61
  )
  lm61_med <- stats::lm(
    depsympB ~ income + I(income^2) + age + woman + race + jobsat,
    data = ah61
  )
  s61 <- suest(
    lm61_base,
    lm61_med,
    model_names = c("Base", "Mediator")
  )
  check_basic_invariants(s61)
  s61
})

test_case("Paper 6.1: income +1 SD benchmark", {
  amount <- stats::sd(ah61$income)
  ame61 <- marginaleffects::avg_comparisons(
    s61,
    variables = list(income = forward_change(amount)),
    newdata = ah61
  )
  diff61 <- marginaleffects::hypotheses(
    ame61,
    hypothesis = difference ~ revpairwise
  )

  expect_near(
    ame61$estimate,
    c(-0.816, -0.648),
    tolerance = 0.005,
    label = "Example 6.1 estimates"
  )
  expect_near(
    ame61$std.error,
    c(0.076, 0.074),
    tolerance = 0.005,
    label = "Example 6.1 standard errors"
  )
  expect_near(
    c(diff61$estimate, diff61$std.error),
    c(-0.168, 0.022),
    tolerance = 0.005,
    label = "Example 6.1 difference"
  )

  list(effects = ame61, difference = diff61)
})

############################################################
# Load GSS once for Examples 6.2 to 6.6
############################################################

test_case("Paper examples: load GSS data", {
  gss <- load_stata_data("gss_cme.dta")
  expect_true(nrow(gss) > 0, "GSS data did not load.")
  dim(gss)
})

############################################################
# Example 6.2: Nested logit models
############################################################

test_case("Paper 6.2: prepare data", {
  vars62 <- c(
    "vhappy", "college", "wages", "occprest", "age",
    "married", "parent", "woman", "conserv", "reltrad", "year",
    "employed"
  )
  gss62 <- subset(gss, year >= 2000 & employed == 1)
  gss62 <- complete_data(gss62, vars62)
  gss62 <- factorize(
    gss62,
    c(
      "college", "married", "parent", "woman",
      "conserv", "reltrad", "year"
    )
  )
  expect_true(nrow(gss62) == 9216, "Example 6.2 should have N=9216.")
  nrow(gss62)
})

test_case("Paper 6.2: nested logit benchmark", {
  logit62_1 <- stats::glm(
    vhappy ~ college,
    family = stats::binomial("logit"),
    data = gss62
  )
  logit62_2 <- stats::glm(
    vhappy ~ college + married + parent + woman + conserv +
      reltrad + year + age + I(age^2),
    family = stats::binomial("logit"),
    data = gss62
  )
  s62 <- suest(
    logit62_1,
    logit62_2,
    model_names = c("Model 1", "Model 2")
  )
  ame62 <- marginaleffects::avg_comparisons(
    s62,
    variables = "college",
    newdata = gss62
  )
  diff62 <- marginaleffects::hypotheses(
    ame62,
    hypothesis = difference ~ revpairwise
  )

  expect_near(
    ame62$estimate,
    c(0.0718060, 0.0599196),
    tolerance = 1e-4,
    label = "Example 6.2 estimates"
  )
  expect_near(
    ame62$std.error,
    c(0.010334, 0.010530),
    tolerance = 2e-4,
    label = "Example 6.2 standard errors"
  )
  expect_near(
    c(diff62$estimate, diff62$std.error),
    c(0.0118864, 0.004023),
    tolerance = 2e-4,
    label = "Example 6.2 difference"
  )

  list(object = s62, effects = ame62, difference = diff62)
})

############################################################
# Example 6.3: Alternative predictor operationalizations
############################################################

test_case("Paper 6.3: prepare data", {
  vars63 <- c(
    "samesexB", "sexident", "sexbehav", "college",
    "woman", "race", "age", "year"
  )
  gss63 <- complete_data(gss, vars63)
  gss63 <- factorize(
    gss63,
    c("sexident", "sexbehav", "college", "woman", "race", "year")
  )
  expect_true(nrow(gss63) == 4921, "Example 6.3 should have N=4921.")
  nrow(gss63)
})

test_case("Paper 6.3: alternative predictors benchmark", {
  logit63_behavior <- stats::glm(
    samesexB ~ sexbehav + woman + college + age + race + year,
    family = stats::binomial("logit"),
    data = gss63
  )
  logit63_identity <- stats::glm(
    samesexB ~ sexident + woman + college + age + race + year,
    family = stats::binomial("logit"),
    data = gss63
  )
  s63 <- suest(
    logit63_behavior,
    logit63_identity,
    model_names = c("Behavior", "Identity")
  )

  ame63 <- marginaleffects::avg_comparisons(
    s63,
    variables = list(
      sexbehav = "reference",
      sexident = "reference"
    ),
    newdata = gss63
  )

  behavior <- ame63[
    as.character(ame63$term) == "sexbehav" &
      as.character(ame63$group) == "Behavior",
  ]
  identity <- ame63[
    as.character(ame63$term) == "sexident" &
      as.character(ame63$group) == "Identity",
  ]

  expect_near(
    behavior$estimate,
    c(-0.097, -0.362),
    tolerance = 0.003,
    label = "Example 6.3 behavior estimates"
  )
  expect_near(
    behavior$std.error,
    c(0.023, 0.038),
    tolerance = 0.003,
    label = "Example 6.3 behavior SEs"
  )
  expect_near(
    identity$estimate,
    c(-0.274, -0.428),
    tolerance = 0.003,
    label = "Example 6.3 identity estimates"
  )
  expect_near(
    identity$std.error,
    c(0.041, 0.034),
    tolerance = 0.003,
    label = "Example 6.3 identity SEs"
  )

  contrast63 <- marginaleffects::avg_comparisons(
    s63,
    variables = list(
      sexbehav = "reference",
      sexident = "reference"
    ),
    newdata = gss63,
    hypothesis = alternative_predictor_hypothesis
  )

  expect_near(
    contrast63$estimate,
    c(0.177, 0.066, -0.111),
    tolerance = 0.004,
    label = "Example 6.3 cross-model contrasts"
  )
  expect_near(
    contrast63$std.error,
    c(0.043, 0.041, 0.059),
    tolerance = 0.004,
    label = "Example 6.3 cross-model contrast SEs"
  )

  list(effects = ame63, contrasts = contrast63)
})

############################################################
# Example 6.4: Negative binomial models for different outcomes
############################################################

test_case("Paper 6.4: prepare data", {
  vars64 <- c(
    "mntlhlth", "physhlth", "woman", "married", "age",
    "faminc", "race", "college", "parent", "reltrad", "year"
  )
  gss64 <- complete_data(gss, vars64)
  gss64 <- factorize(
    gss64,
    c("woman", "married", "parent", "college", "race", "year")
  )
  expect_true(nrow(gss64) == 5062, "Example 6.4 should have N=5062.")
  nrow(gss64)
})

test_case("Paper 6.4: negative-binomial outcome benchmark", {
  nb64_mental <- MASS::glm.nb(
    mntlhlth ~ woman + married + parent + college + age +
      faminc + race + year,
    data = gss64
  )
  nb64_physical <- MASS::glm.nb(
    physhlth ~ woman + married + parent + college + age +
      faminc + race + year,
    data = gss64
  )
  s64 <- suest(
    nb64_mental,
    nb64_physical,
    model_names = c("Mental", "Physical")
  )

  ame64 <- marginaleffects::avg_comparisons(
    s64,
    variables = list(
      woman = "reference",
      married = "reference",
      parent = "reference",
      college = "reference",
      age = forward_change(stats::sd(gss64$age)),
      faminc = forward_change(stats::sd(gss64$faminc)),
      race = "reference",
      year = "reference"
    ),
    newdata = gss64
  )

  diff64 <- marginaleffects::avg_comparisons(
    s64,
    variables = list(
      woman = "reference",
      married = "reference",
      parent = "reference",
      college = "reference",
      age = forward_change(stats::sd(gss64$age)),
      faminc = forward_change(stats::sd(gss64$faminc)),
      race = "reference",
      year = "reference"
    ),
    newdata = gss64,
    hypothesis = pair_by_term_hypothesis("Mental", "Physical")
  )

  print(ame64[
    as.character(ame64$term) %in%
      c("woman", "married", "parent", "college", "age", "faminc"),
  ])
  print(diff64[
    grepl(
      "^(woman|married|parent|college|age|faminc):",
      as.character(diff64$term)
    ),
  ])

  expected <- list(
    woman = list(
      estimate = c(0.993, 0.773),
      se = c(0.208, 0.174),
      difference = c(0.220, 0.229)
    ),
    married = list(
      estimate = c(-1.010, -0.159),
      se = c(0.230, 0.192),
      difference = c(-0.851, 0.250)
    ),
    parent = list(
      estimate = c(0.274, -0.269),
      se = c(0.248, 0.219),
      difference = c(0.543, 0.274)
    ),
    college = list(
      estimate = c(-0.879, -0.542),
      se = c(0.229, 0.189),
      difference = c(-0.337, 0.254)
    ),
    age = list(
      estimate = c(-0.462, 0.491),
      se = c(0.108, 0.118),
      difference = c(-0.953, 0.132)
    ),
    faminc = list(
      estimate = c(-0.445, -0.381),
      se = c(0.118, 0.102),
      difference = c(-0.064, 0.133)
    )
  )

  for (term in names(expected)) {
    rows <- ame64[as.character(ame64$term) == term, ]
    rows <- rows[match(c("Mental", "Physical"), as.character(rows$group)), ]

    expect_near(
      rows$estimate,
      expected[[term]]$estimate,
      tolerance = 0.012,
      label = paste("Example 6.4", term, "estimates")
    )
    expect_near(
      rows$std.error,
      expected[[term]]$se,
      tolerance = 0.012,
      label = paste("Example 6.4", term, "SEs")
    )

    difference_row <- diff64[
      grepl(paste0("^", term, ":"), as.character(diff64$term)),
    ]

    expect_near(
      c(difference_row$estimate, difference_row$std.error),
      expected[[term]]$difference,
      tolerance = 0.015,
      label = paste("Example 6.4", term, "difference")
    )
  }

  # Full race and year benchmark vectors.
  race_mental <- ame64[
    as.character(ame64$term) == "race" &
      as.character(ame64$group) == "Mental",
  ]
  race_physical <- ame64[
    as.character(ame64$term) == "race" &
      as.character(ame64$group) == "Physical",
  ]
  year_mental <- ame64[
    as.character(ame64$term) == "year" &
      as.character(ame64$group) == "Mental",
  ]
  year_physical <- ame64[
    as.character(ame64$term) == "year" &
      as.character(ame64$group) == "Physical",
  ]

  expect_near(
    race_mental$estimate,
    c(-1.016, -0.438),
    0.015,
    "Example 6.4 race mental estimates"
  )
  expect_near(
    race_physical$estimate,
    c(-0.531, 0.145),
    0.015,
    "Example 6.4 race physical estimates"
  )
  expect_near(
    year_mental$estimate,
    c(-0.944, -0.068, -0.582),
    0.015,
    "Example 6.4 year mental estimates"
  )
  expect_near(
    year_physical$estimate,
    c(-0.292, 0.323, -0.229),
    0.015,
    "Example 6.4 year physical estimates"
  )

  list(effects = ame64, differences = diff64, theta = s64$theta)
})

############################################################
# Example 6.5: Ordered logit versus multinomial logit
############################################################

test_case("Paper 6.5: prepare data", {
  vars65 <- c(
    "partyid5", "woman", "edyrs", "age", "parent", "married",
    "faminc", "employed", "region4", "year", "race"
  )
  gss65 <- subset(gss, year >= 2010)
  gss65 <- complete_data(gss65, vars65)

  factor_vars65 <- c(
    "woman", "parent", "married", "race",
    "employed", "region4", "year"
  )
  gss65 <- factorize(gss65, factor_vars65)
  levels_party <- sort(unique(gss65$partyid5))
  gss65$party_ord <- ordered(gss65$partyid5, levels = levels_party)
  gss65$party_nom <- factor(gss65$partyid5, levels = levels_party)

  expect_true(nrow(gss65) == 8179, "Example 6.5 should have N=8179.")
  nrow(gss65)
})

test_case("Paper 6.5: ordered versus nominal benchmark at covariate means", {
  ologit65 <- MASS::polr(
    party_ord ~ age + I(age^2) + woman + edyrs + parent +
      married + race + faminc + employed + region4 + year,
    data = gss65,
    method = "logistic",
    Hess = TRUE
  )
  mlogit65 <- nnet::multinom(
    party_nom ~ age + I(age^2) + woman + edyrs + parent +
      married + race + faminc + employed + region4 + year,
    data = gss65,
    Hess = TRUE,
    trace = FALSE
  )
  s65 <- suest(
    ologit65,
    mlogit65,
    model_names = c("Ordered", "Multinomial")
  )

  results65 <- paper65_results(s65, 20, 30)
  effects65 <- results65$effects
  differences65 <- results65$differences

  ordered <- effects65$estimate[1:5]
  nominal <- effects65$estimate[6:10]
  se65 <- effects65$std.error
  differences <- differences65$estimate
  se_diff <- differences65$std.error

  J_ordered65 <- results65$jacobian[
    1:5,
    s65$index[[1]],
    drop = FALSE
  ]
  J_nominal65 <- results65$jacobian[
    6:10,
    s65$index[[2]],
    drop = FALSE
  ]
  ordered_model_se65 <- sqrt(diag(
    J_ordered65 %*% stats::vcov(ologit65) %*% t(J_ordered65)
  ))
  nominal_model_se65 <- sqrt(diag(
    J_nominal65 %*% stats::vcov(mlogit65) %*% t(J_nominal65)
  ))

  print(data.frame(
    category = ologit65$lev,
    ordered = ordered,
    ordered_robust_se = se65[1:5],
    ordered_model_se = ordered_model_se65,
    nominal = nominal,
    nominal_robust_se = se65[6:10],
    nominal_model_se = nominal_model_se65,
    difference = differences,
    difference_se = se_diff
  ))

  expect_near(
    ordered,
    c(0.020, 0.023, -0.003, -0.024, -0.016),
    tolerance = 0.003,
    label = "Example 6.5 ordered effects"
  )
  expect_near(
    se65[1:5],
    c(0.004, 0.005, 0.000, 0.005, 0.004),
    tolerance = 0.003,
    label = "Example 6.5 ordered SEs"
  )
  expect_near(
    nominal,
    c(0.032, -0.010, 0.000, -0.031, 0.010),
    tolerance = 0.004,
    label = "Example 6.5 nominal effects"
  )
  expect_near(
    se65[6:10],
    c(0.003, 0.011, 0.010, 0.011, 0.003),
    tolerance = 0.003,
    label = "Example 6.5 nominal SEs"
  )
  expect_near(
    differences,
    c(-0.012, 0.034, -0.003, 0.007, -0.026),
    tolerance = 0.004,
    label = "Example 6.5 differences"
  )
  expect_near(
    se_diff,
    c(0.003, 0.009, 0.010, 0.009, 0.004),
    tolerance = 0.003,
    label = "Example 6.5 difference SEs"
  )

  data.frame(
    category = ologit65$lev,
    ordered = ordered,
    ordered_se = se65[1:5],
    nominal = nominal,
    nominal_se = se65[6:10],
    difference = differences,
    difference_se = se_diff
  )
})

############################################################
# Example 6.6: Different samples
############################################################

test_case("Paper 6.6: prepare data and fit disjoint samples", {
  vars66 <- c(
    "helpsickB", "polviews", "conserv", "faminc", "employed",
    "woman", "age", "college", "married", "parent", "race", "year"
  )
  gss66 <- complete_data(gss, vars66)
  gss66 <- factorize(
    gss66,
    c(
      "conserv", "employed", "woman", "college",
      "married", "parent", "race"
    )
  )

  logit66_1986 <- stats::glm(
    helpsickB ~ conserv + faminc + employed + woman + age +
      college + married + parent + race,
    family = stats::binomial("logit"),
    data = gss66,
    subset = year == 1986
  )
  logit66_2016 <- stats::glm(
    helpsickB ~ conserv + faminc + employed + woman + age +
      college + married + parent + race,
    family = stats::binomial("logit"),
    data = gss66,
    subset = year == 2016
  )
  s66 <- suest(
    logit66_1986,
    logit66_2016,
    model_names = c("1986", "2016")
  )

  expect_true(s66$nobs_models[1] == 1254, "1986 N should be 1254.")
  expect_true(s66$nobs_models[2] == 1670, "2016 N should be 1670.")
  expect_true(s66$nobs_overlap == 0, "Samples should be disjoint.")
  expect_near(
    offdiag_vcov(s66),
    matrix(0, nrow(offdiag_vcov(s66)), ncol(offdiag_vcov(s66))),
    tolerance = 0,
    label = "Example 6.6 cross covariance"
  )

  s66
})

test_case("Paper 6.6: different-sample benchmark", {
  nd66 <- suest_newdata(s66)
  ame66 <- marginaleffects::avg_comparisons(
    s66,
    variables = "conserv",
    newdata = nd66
  )
  diff66 <- marginaleffects::hypotheses(
    ame66,
    hypothesis = difference ~ revpairwise
  )

  expect_near(
    ame66$estimate,
    c(-0.092, -0.258),
    tolerance = 0.003,
    label = "Example 6.6 estimates"
  )
  expect_near(
    ame66$std.error,
    c(0.030, 0.025),
    tolerance = 0.003,
    label = "Example 6.6 SEs"
  )
  expect_near(
    c(diff66$estimate, diff66$std.error),
    c(0.166, 0.039),
    tolerance = 0.004,
    label = "Example 6.6 difference"
  )

  list(effects = ame66, difference = diff66)
})
