adapter_housing_data <- function() {
  data <- MASS::housing[
    rep(seq_len(nrow(MASS::housing)), MASS::housing$Freq),
    c("Sat", "Infl", "Type", "Cont")
  ]
  rownames(data) <- NULL
  data$Sat <- ordered(
    data$Sat,
    levels = c("Low", "Medium", "High")
  )
  data$woman <- rep(c(0, 1), length.out = nrow(data))
  data
}

test_that("glm2 logit, probit, and Poisson models are supported", {
  skip_if_not_installed("glm2")

  dat <- mtcars
  dat$am <- factor(dat$am)

  fits <- list(
    logit = list(
      glm2::glm2(
        am ~ wt,
        family = binomial("logit"),
        data = dat
      ),
      glm2::glm2(
        am ~ wt + hp,
        family = binomial("logit"),
        data = dat
      )
    ),
    probit = list(
      glm2::glm2(
        am ~ wt,
        family = binomial("probit"),
        data = dat
      ),
      glm2::glm2(
        am ~ wt + hp,
        family = binomial("probit"),
        data = dat
      )
    ),
    poisson = list(
      glm2::glm2(
        cyl ~ wt,
        family = poisson("log"),
        data = dat
      ),
      glm2::glm2(
        cyl ~ wt + hp,
        family = poisson("log"),
        data = dat
      )
    )
  )

  for (family_name in names(fits)) {
    combined <- suest(
      fits[[family_name]][[1L]],
      fits[[family_name]][[2L]],
      model_names = c("Base", "Adjusted")
    )

    expect_equal(
      unname(combined$model_engines),
      rep("glm2::glm2", 2L)
    )

    effects <- marginaleffects::avg_comparisons(
      combined,
      variables = "wt",
      newdata = dat
    )

    expect_true(all(is.finite(effects$estimate)))
    expect_true(all(is.finite(effects$std.error)))
  }
})

test_that("glm2 and stats glm can be compared directly", {
  skip_if_not_installed("glm2")

  dat <- mtcars
  dat$am <- factor(dat$am)

  standard <- stats::glm(
    am ~ wt + hp,
    family = binomial("logit"),
    data = dat
  )
  alternative <- glm2::glm2(
    am ~ wt + hp,
    family = binomial("logit"),
    data = dat
  )

  combined <- suest(
    standard,
    alternative,
    model_names = c("glm", "glm2")
  )

  expect_equal(
    unname(combined$model_engines),
    c("stats::glm", "glm2::glm2")
  )

  effects <- marginaleffects::avg_comparisons(
    combined,
    variables = "wt",
    newdata = dat
  )

  expect_equal(
    effects$estimate[1L],
    effects$estimate[2L],
    tolerance = 1e-8
  )
})

test_that("ordinal clm logit and probit models are supported", {
  skip_if_not_installed("ordinal")

  dat <- adapter_housing_data()

  for (link in c("logit", "probit")) {
    base <- ordinal::clm(
      Sat ~ Cont + Infl,
      data = dat,
      link = link
    )
    adjusted <- ordinal::clm(
      Sat ~ Cont + Infl + Type,
      data = dat,
      link = link
    )

    combined <- suest(
      base,
      adjusted,
      model_names = c("Base", "Adjusted")
    )

    expect_equal(
      unname(combined$model_engines),
      rep("ordinal::clm", 2L)
    )
    expect_equal(
      combined$category_levels[[1L]],
      levels(dat$Sat)
    )
    expect_gt(max(abs(offdiag_vcov(combined))), 0)

    predict_method <- getFromNamespace(
      "get_predict.suest_model",
      "suest"
    )
    set_method <- getFromNamespace(
      "set_coef.suest_model",
      "suest"
    )

    raw_predictions <- predict_method(
      combined,
      newdata = utils::head(dat),
      type = "response"
    )
    expect_equal(
      nrow(raw_predictions),
      2L * nlevels(dat$Sat) * nrow(utils::head(dat))
    )

    perturbed <- combined$coefficients
    perturbed[1L] <- perturbed[1L] + 1e-4
    perturbed_model <- set_method(combined, perturbed)
    perturbed_predictions <- predict_method(
      perturbed_model,
      newdata = utils::head(dat),
      type = "response"
    )
    expect_true(any(
      abs(
        perturbed_predictions$estimate -
          raw_predictions$estimate
      ) > 0
    ))

    effects <- marginaleffects::avg_comparisons(
      combined,
      variables = "Cont",
      newdata = dat
    )

    expect_equal(nrow(effects), 2L * nlevels(dat$Sat))
    expect_true(all(is.finite(effects$estimate)))
    expect_true(all(is.finite(effects$std.error)))
  }
})

test_that("clm and polr models can be compared directly", {
  skip_if_not_installed("ordinal")

  dat <- adapter_housing_data()

  clm_fit <- ordinal::clm(
    Sat ~ Cont + Infl + Type,
    data = dat,
    link = "logit"
  )
  polr_fit <- MASS::polr(
    Sat ~ Cont + Infl + Type,
    data = dat,
    method = "logistic",
    Hess = TRUE
  )

  combined <- suest(
    clm_fit,
    polr_fit,
    model_names = c("clm", "polr")
  )

  expect_equal(
    unname(combined$model_types),
    c("ologit", "ologit")
  )
  expect_equal(
    unname(combined$model_engines),
    c("ordinal::clm", "MASS::polr")
  )

  effects <- marginaleffects::avg_comparisons(
    combined,
    variables = "Cont",
    newdata = dat
  )

  expect_equal(nrow(effects), 2L * nlevels(dat$Sat))
  expect_true(all(is.finite(effects$std.error)))
})

test_that("clm sample subsets are aligned correctly", {
  skip_if_not_installed("ordinal")

  dat <- adapter_housing_data()

  men <- ordinal::clm(
    Sat ~ Cont + Infl + Type,
    data = dat,
    subset = woman == 0,
    link = "logit"
  )
  women <- ordinal::clm(
    Sat ~ Cont + Infl + Type,
    data = dat,
    subset = woman == 1,
    link = "logit"
  )

  combined <- suest(
    men,
    women,
    model_names = c("Men", "Women")
  )

  expect_equal(combined$nobs_overlap, 0L)
  expect_true(all(offdiag_vcov(combined) == 0))

  effects <- marginaleffects::avg_comparisons(
    combined,
    variables = "Cont",
    newdata = suest_newdata(combined)
  )

  expect_true(all(is.finite(effects$std.error)))
})

test_that("clm covariance blocks match direct HC1 calculations", {
  skip_if_not_installed("ordinal")

  dat <- adapter_housing_data()
  model1 <- ordinal::clm(
    Sat ~ Cont + Infl,
    data = dat,
    link = "logit"
  )
  model2 <- ordinal::clm(
    Sat ~ Cont + Infl + Type,
    data = dat,
    link = "logit"
  )

  combined <- suest(model1, model2)
  V <- vcov(combined)

  expect_equal(
    unname(V[combined$index[[1L]], combined$index[[1L]]]),
    unname(robust_vcov_direct(model1)),
    tolerance = 1e-8
  )
  expect_equal(
    unname(V[combined$index[[2L]], combined$index[[2L]]]),
    unname(robust_vcov_direct(model2)),
    tolerance = 1e-8
  )
})

test_that("unsupported clm extensions are rejected clearly", {
  skip_if_not_installed("ordinal")

  dat <- adapter_housing_data()
  ordinary <- ordinal::clm(
    Sat ~ Cont + Infl,
    data = dat,
    link = "logit"
  )

  scale_model <- ordinal::clm(
    Sat ~ Cont,
    scale = ~ Infl,
    data = dat,
    link = "logit"
  )
  nominal_model <- ordinal::clm(
    Sat ~ Cont,
    nominal = ~ Infl,
    data = dat,
    link = "logit"
  )
  structured_threshold <- ordinal::clm(
    Sat ~ Cont + Infl,
    data = dat,
    link = "logit",
    threshold = "equidistant"
  )

  expect_error(
    suest(ordinary, scale_model),
    "scale models"
  )
  expect_error(
    suest(ordinary, nominal_model),
    "nominal-effects"
  )
  expect_error(
    suest(ordinary, structured_threshold),
    "threshold = 'flexible'"
  )
})

test_that("adjusted-score GLM fits are rejected", {
  skip_if_not_installed("brglm2")

  dat <- mtcars
  dat$am <- factor(dat$am)

  ordinary <- stats::glm(
    am ~ wt,
    family = binomial(),
    data = dat
  )
  adjusted <- stats::glm(
    am ~ wt,
    family = binomial(),
    data = dat,
    method = brglm2::brglmFit
  )

  expect_error(
    suest(ordinary, adjusted),
    "Bias-reduced|adjusted-score|Firth|penalized"
  )
})
