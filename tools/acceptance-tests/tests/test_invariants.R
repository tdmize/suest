cat("\n\nCORE AND MODEL-FAMILY INVARIANT TESTS\n")

############################################################
# Linear models
############################################################

test_case("Invariants: two linear models", {
  data <- datasets::mtcars
  m1 <- stats::lm(mpg ~ wt, data = data)
  m2 <- stats::lm(mpg ~ wt + hp, data = data)
  s <- suest(m1, m2, model_names = c("L1", "L2"))

  check_basic_invariants(s)
  expect_near(
    vcov_block(s, 1),
    robust_vcov_direct(m1),
    1e-10,
    "linear block 1"
  )
  expect_near(
    vcov_block(s, 2),
    robust_vcov_direct(m2),
    1e-10,
    "linear block 2"
  )

  combined <- marginaleffects::avg_slopes(
    s,
    variables = "wt",
    newdata = data
  )
  separate <- c(
    marginaleffects::avg_slopes(
      m1,
      variables = "wt",
      newdata = data,
      vcov = vcov_block(s, 1)
    )$estimate,
    marginaleffects::avg_slopes(
      m2,
      variables = "wt",
      newdata = data,
      vcov = vcov_block(s, 2)
    )$estimate
  )
  expect_near(combined$estimate, separate, 1e-10, "linear slopes")
  s
})

############################################################
# Binary logit, probit, and cross-link comparisons
############################################################

test_case("Invariants: logit, probit, and mixed logit-probit", {
  data <- datasets::mtcars
  data$am <- factor(data$am)
  logit1 <- stats::glm(
    am ~ wt,
    family = stats::binomial("logit"),
    data = data
  )
  logit2 <- stats::glm(
    am ~ wt + hp,
    family = stats::binomial("logit"),
    data = data
  )
  probit2 <- stats::glm(
    am ~ wt + hp,
    family = stats::binomial("probit"),
    data = data
  )

  s_logit <- suest(logit1, logit2, c("G1", "G2"))
  s_mixed <- suest(logit1, probit2, c("Logit", "Probit"))
  check_basic_invariants(s_logit)
  check_basic_invariants(s_mixed)

  effects <- marginaleffects::avg_comparisons(
    s_mixed,
    variables = "wt",
    newdata = data
  )
  expect_true(nrow(effects) == 2L)
  expect_true(all(is.finite(effects$estimate)))
  list(logit = s_logit, mixed = effects)
})

############################################################
# Poisson and negative binomial
############################################################

test_case("Invariants: Poisson, negative binomial, and mixed count models", {
  data <- datasets::warpbreaks
  p1 <- stats::glm(
    breaks ~ wool,
    family = stats::poisson("log"),
    data = data
  )
  p2 <- stats::glm(
    breaks ~ wool + tension,
    family = stats::poisson("log"),
    data = data
  )
  n2 <- MASS::glm.nb(
    breaks ~ wool + tension,
    data = data
  )

  s_pois <- suest(p1, p2, c("P1", "P2"))
  s_count <- suest(p1, n2, c("Poisson", "NB"))
  check_basic_invariants(s_pois)
  check_basic_invariants(s_count)

  effects <- marginaleffects::avg_comparisons(
    s_count,
    variables = "wool",
    newdata = data
  )
  expect_true(nrow(effects) == 2L)
  expect_true(all(is.finite(effects$estimate)))
  list(poisson = s_pois, mixed = effects)
})

test_case("Invariants: full negative-binomial score and dispersion", {
  data <- datasets::warpbreaks
  model <- MASS::glm.nb(
    breaks ~ wool + tension,
    data = data
  )
  components <- suest:::.suest_negbin_components(model)

  expect_true(
    ncol(components$score) == length(stats::coef(model)) + 1L,
    "Negative-binomial score must include log(theta)."
  )
  expect_true(
    identical(
      colnames(components$score),
      names(components$parameters)
    ),
    "Negative-binomial score and parameter names do not match."
  )
  expect_near(
    colSums(components$score),
    rep(0, ncol(components$score)),
    tolerance = 2e-5,
    label = "negative-binomial score sums"
  )
  expect_true(
    all(is.finite(components$bread)),
    "Negative-binomial bread contains nonfinite values."
  )

  list(
    theta = model$theta,
    score_sums = colSums(components$score)
  )
})

test_case("Invariants: analytic ordered and multinomial components", {
  housing <- MASS::housing[
    rep(seq_len(nrow(MASS::housing)), MASS::housing$Freq),
    c("Sat", "Infl", "Type", "Cont")
  ]
  rownames(housing) <- NULL
  housing$Sat_ord <- ordered(
    housing$Sat,
    levels = c("Low", "Medium", "High")
  )
  housing$Sat_nom <- factor(
    housing$Sat,
    levels = c("Low", "Medium", "High")
  )

  ordered <- MASS::polr(
    Sat_ord ~ Cont + Infl + Type,
    data = housing,
    method = "logistic",
    Hess = TRUE
  )
  nominal <- nnet::multinom(
    Sat_nom ~ Cont + Infl + Type,
    data = housing,
    Hess = TRUE,
    trace = FALSE
  )

  ordered_components <- suest:::.suest_polr_components(ordered)
  nominal_components <- suest:::.suest_multinom_components(nominal)

  expect_near(
    ordered_components$score,
    sandwich::estfun(ordered),
    tolerance = 1e-10,
    label = "ordered analytic scores"
  )
  nominal_sandwich_score <- tryCatch(
    sandwich::estfun(nominal),
    error = function(e) NULL
  )

  if (!is.null(nominal_sandwich_score)) {
    expect_near(
      nominal_components$score,
      nominal_sandwich_score,
      tolerance = 1e-10,
      label = "multinomial analytic scores"
    )
  }
  expect_near(
    colMeans(ordered_components$score),
    rep(0, ncol(ordered_components$score)),
    tolerance = 1e-7,
    label = "ordered mean scores"
  )
  expect_near(
    colMeans(nominal_components$score),
    rep(0, ncol(nominal_components$score)),
    tolerance = 1e-6,
    label = "multinomial mean scores"
  )
  expect_true(all(is.finite(ordered_components$bread)))
  expect_true(all(is.finite(nominal_components$bread)))

  list(
    ordered_mean_scores = colMeans(ordered_components$score),
    multinomial_mean_scores = colMeans(nominal_components$score),
    sandwich_multinom_score_available =
      !is.null(nominal_sandwich_score)
  )
})

############################################################
# Ordered logit and ordered probit
############################################################

test_case("Invariants: ordered logit and ordered probit", {
  data <- MASS::housing[
    rep(seq_len(nrow(MASS::housing)), MASS::housing$Freq),
    c("Sat", "Infl", "Type", "Cont")
  ]
  rownames(data) <- NULL
  data$Sat <- ordered(data$Sat, levels = c("Low", "Medium", "High"))

  ol1 <- MASS::polr(
    Sat ~ Cont,
    data = data,
    method = "logistic",
    Hess = TRUE
  )
  ol2 <- MASS::polr(
    Sat ~ Cont + Infl + Type,
    data = data,
    method = "logistic",
    Hess = TRUE
  )
  op1 <- MASS::polr(
    Sat ~ Cont,
    data = data,
    method = "probit",
    Hess = TRUE
  )
  op2 <- MASS::polr(
    Sat ~ Cont + Infl + Type,
    data = data,
    method = "probit",
    Hess = TRUE
  )

  s_ol <- suest(ol1, ol2, c("OL1", "OL2"))
  s_op <- suest(op1, op2, c("OP1", "OP2"))
  check_basic_invariants(s_ol)
  check_basic_invariants(s_op)

  for (object in list(s_ol, s_op)) {
    predictions <- marginaleffects::avg_predictions(
      object,
      newdata = data
    )
    totals <- aggregate(
      predictions$estimate,
      list(model = sub(
        "::.*$",
        "",
        as.character(predictions$group)
      )),
      sum
    )
    expect_near(totals$x, c(1, 1), 1e-10, "ordinal probabilities")

    effects <- marginaleffects::avg_comparisons(
      object,
      variables = "Cont",
      newdata = data
    )
    sums <- aggregate(
      effects$estimate,
      list(model = sub(
        "::.*$",
        "",
        as.character(effects$group)
      )),
      sum
    )
    expect_near(sums$x, c(0, 0), 1e-10, "ordinal effect sums")
  }

  list(ordered_logit = s_ol, ordered_probit = s_op)
})

############################################################
# Multinomial logit
############################################################

test_case("Invariants: multinomial logit", {
  data <- datasets::iris
  data$wide <- factor(
    ifelse(
      data$Sepal.Width > stats::median(data$Sepal.Width),
      "Wide",
      "Narrow"
    ),
    levels = c("Narrow", "Wide")
  )

  m1 <- nnet::multinom(
    Species ~ wide,
    data = data,
    Hess = TRUE,
    trace = FALSE
  )
  m2 <- nnet::multinom(
    Species ~ wide + Sepal.Length,
    data = data,
    Hess = TRUE,
    trace = FALSE
  )
  s <- suest(m1, m2, c("MN1", "MN2"))
  check_basic_invariants(s)

  predictions <- marginaleffects::avg_predictions(s, newdata = data)
  totals <- aggregate(
    predictions$estimate,
    list(model = sub(
      "::.*$",
      "",
      as.character(predictions$group)
    )),
    sum
  )
  expect_near(totals$x, c(1, 1), 1e-10, "multinomial probabilities")

  effects <- marginaleffects::avg_comparisons(
    s,
    variables = "wide",
    newdata = data
  )
  sums <- aggregate(
    effects$estimate,
    list(model = sub(
      "::.*$",
      "",
      as.character(effects$group)
    )),
    sum
  )
  expect_near(sums$x, c(0, 0), 1e-10, "multinomial effect sums")
  s
})

############################################################
# Sample alignment
############################################################

test_case("Invariants: identical, partial, and disjoint samples", {
  set.seed(376)
  data <- data.frame(
    id = seq_len(600),
    year = rep(c(1986, 2016), each = 300),
    x = stats::rnorm(600),
    z = stats::rnorm(600)
  )
  data$y <- stats::rbinom(
    600,
    1,
    stats::plogis(-0.2 + 0.4 * data$x - 0.2 * data$z)
  )

  same1 <- stats::glm(
    y ~ x,
    family = stats::binomial(),
    data = data
  )
  same2 <- stats::glm(
    y ~ x + z,
    family = stats::binomial(),
    data = data
  )
  partial1 <- stats::glm(
    y ~ x,
    family = stats::binomial(),
    data = data,
    subset = id <= 450
  )
  partial2 <- stats::glm(
    y ~ x + z,
    family = stats::binomial(),
    data = data,
    subset = id >= 151
  )
  disjoint1 <- stats::glm(
    y ~ x + z,
    family = stats::binomial(),
    data = data,
    subset = year == 1986
  )
  disjoint2 <- stats::glm(
    y ~ x + z,
    family = stats::binomial(),
    data = data,
    subset = year == 2016
  )

  s_same <- suest(same1, same2)
  s_partial <- suest(partial1, partial2)
  s_disjoint <- suest(disjoint1, disjoint2)

  expect_true(s_same$nobs_overlap == 600)
  expect_true(max(abs(offdiag_vcov(s_same))) > 0)
  expect_true(s_partial$nobs_overlap == 300)
  expect_true(s_partial$nobs_union == 600)
  expect_true(max(abs(offdiag_vcov(s_partial))) > 0)
  expect_true(s_disjoint$nobs_overlap == 0)
  expect_near(
    offdiag_vcov(s_disjoint),
    matrix(
      0,
      nrow(offdiag_vcov(s_disjoint)),
      ncol(offdiag_vcov(s_disjoint))
    ),
    0,
    "disjoint covariance"
  )

  expect_near(
    vcov_block(s_disjoint, 1),
    robust_vcov_direct(disjoint1),
    1e-10,
    "disjoint block 1"
  )
  expect_near(
    vcov_block(s_disjoint, 2),
    robust_vcov_direct(disjoint2),
    1e-10,
    "disjoint block 2"
  )

  list(
    same = s_same$nobs_overlap,
    partial = s_partial$nobs_overlap,
    disjoint = s_disjoint$nobs_overlap
  )
})

############################################################
# Reordering and renaming
############################################################

test_case("Invariants: model reordering and renaming", {
  data <- datasets::mtcars
  data$am <- factor(data$am)
  m1 <- stats::glm(
    am ~ wt,
    family = stats::binomial(),
    data = data
  )
  m2 <- stats::glm(
    am ~ wt + hp,
    family = stats::binomial(),
    data = data
  )

  s12 <- suest(m1, m2, c("A", "B"))
  s21 <- suest(m2, m1, c("B", "A"))
  renamed <- suest(m1, m2, c("First", "Second"))

  effects12 <- marginaleffects::avg_comparisons(
    s12,
    variables = "wt",
    newdata = data
  )
  effects21 <- marginaleffects::avg_comparisons(
    s21,
    variables = "wt",
    newdata = data
  )

  order21 <- match(
    as.character(effects12$group),
    as.character(effects21$group)
  )

  expect_true(
    !anyNA(order21),
    "Reordered model labels could not be matched."
  )
  expect_near(
    effects12$estimate,
    effects21$estimate[order21],
    1e-10,
    "reordered estimates"
  )
  expect_near(
    effects12$std.error,
    effects21$std.error[order21],
    1e-8,
    "reordered SEs"
  )
  expect_equal(
    unname(stats::coef(s12)),
    unname(stats::coef(renamed)),
    tolerance = 0
  )

  TRUE
})
