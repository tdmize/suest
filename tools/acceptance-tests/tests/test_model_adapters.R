############################################################
# Alternative model-fitting engines
############################################################

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
  data$id <- seq_len(nrow(data))
  data
}

test_case("Model adapters: prepare common data", {
  adapter_data <- adapter_housing_data()
  assign("adapter_data", adapter_data, envir = globalenv())

  expect_true(nrow(adapter_data) == 1681L)
  expect_true(all(table(adapter_data$Sat) > 0L))

  table(adapter_data$Sat)
})

for (link in c("logit", "probit")) {
  local({
    link <- link

    test_case(
      paste0("Model adapters: glm2 binary ", link),
      {
        dat <- mtcars
        dat$am <- factor(dat$am)
        assign("glm2_dat", dat, envir = globalenv())

        base <- glm2::glm2(
          am ~ wt,
          family = stats::binomial(link),
          data = glm2_dat
        )
        adjusted <- glm2::glm2(
          am ~ wt + hp,
          family = stats::binomial(link),
          data = glm2_dat
        )

        combined <- suest(
          base,
          adjusted,
          model_names = c("Base", "Adjusted")
        )

        expect_true(all(
          unname(combined$model_engines) == "glm2::glm2"
        ))
        expect_true(combined$nobs_overlap == nrow(glm2_dat))
        expect_true(max(abs(offdiag_vcov(combined))) > 0)

        effects <- marginaleffects::avg_comparisons(
          combined,
          variables = "wt",
          newdata = glm2_dat
        )
        difference <- marginaleffects::hypotheses(
          effects,
          hypothesis = difference ~ revpairwise
        )

        expect_true(all(is.finite(effects$estimate)))
        expect_true(all(is.finite(effects$std.error)))
        expect_true(all(is.finite(difference$estimate)))
        expect_true(all(is.finite(difference$std.error)))

        list(
          link = link,
          engines = combined$model_engines,
          difference = difference
        )
      }
    )
  })
}

test_case("Model adapters: glm2 Poisson", {
  dat <- mtcars
  assign("glm2_count_dat", dat, envir = globalenv())

  base <- glm2::glm2(
    cyl ~ wt,
    family = stats::poisson("log"),
    data = glm2_count_dat
  )
  adjusted <- glm2::glm2(
    cyl ~ wt + hp,
    family = stats::poisson("log"),
    data = glm2_count_dat
  )

  combined <- suest(
    base,
    adjusted,
    model_names = c("Base", "Adjusted")
  )

  effects <- marginaleffects::avg_comparisons(
    combined,
    variables = "wt",
    newdata = glm2_count_dat
  )
  difference <- marginaleffects::hypotheses(
    effects,
    hypothesis = difference ~ revpairwise
  )

  expect_true(all(is.finite(difference$std.error)))

  list(
    engines = combined$model_engines,
    difference = difference
  )
})

test_case("Model adapters: stats glm versus glm2", {
  dat <- mtcars
  dat$am <- factor(dat$am)
  assign("glm_cross_dat", dat, envir = globalenv())

  standard <- stats::glm(
    am ~ wt + hp,
    family = stats::binomial("logit"),
    data = glm_cross_dat
  )
  alternative <- glm2::glm2(
    am ~ wt + hp,
    family = stats::binomial("logit"),
    data = glm_cross_dat
  )

  combined <- suest(
    standard,
    alternative,
    model_names = c("glm", "glm2")
  )

  expect_true(identical(
    unname(combined$model_engines),
    c("stats::glm", "glm2::glm2")
  ))

  effects <- marginaleffects::avg_comparisons(
    combined,
    variables = "wt",
    newdata = glm_cross_dat
  )

  expect_near(
    effects$estimate[1L],
    effects$estimate[2L],
    tolerance = 1e-8,
    label = "glm versus glm2 effects"
  )

  effects
})

for (link in c("logit", "probit")) {
  local({
    link <- link

    test_case(
      paste0("Model adapters: ordinal clm ", link),
      {
        base <- ordinal::clm(
          Sat ~ Cont + Infl,
          data = adapter_data,
          link = link
        )
        adjusted <- ordinal::clm(
          Sat ~ Cont + Infl + Type,
          data = adapter_data,
          link = link
        )

        combined <- suest(
          base,
          adjusted,
          model_names = c("Base", "Adjusted")
        )

        expect_true(all(
          unname(combined$model_engines) == "ordinal::clm"
        ))
        expect_true(identical(
          combined$category_levels[[1L]],
          levels(adapter_data$Sat)
        ))
        expect_true(max(abs(offdiag_vcov(combined))) > 0)

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
          newdata = utils::head(adapter_data),
          type = "response"
        )
        expect_true(
          nrow(raw_predictions) ==
            2L * nlevels(adapter_data$Sat) * nrow(utils::head(adapter_data))
        )

        perturbed <- combined$coefficients
        perturbed[1L] <- perturbed[1L] + 1e-4
        perturbed_model <- set_method(combined, perturbed)
        perturbed_predictions <- predict_method(
          perturbed_model,
          newdata = utils::head(adapter_data),
          type = "response"
        )
        expect_true(any(
          abs(
            perturbed_predictions$estimate -
              raw_predictions$estimate
          ) > 0
        ))

        direct1 <- robust_vcov_direct(base)
        direct2 <- robust_vcov_direct(adjusted)

        expect_near(
          unname(combined$vcov[
            combined$index[[1L]],
            combined$index[[1L]]
          ]),
          unname(direct1),
          tolerance = 1e-8,
          label = paste0("clm ", link, " model 1 covariance")
        )
        expect_near(
          unname(combined$vcov[
            combined$index[[2L]],
            combined$index[[2L]]
          ]),
          unname(direct2),
          tolerance = 1e-8,
          label = paste0("clm ", link, " model 2 covariance")
        )

        effects <- marginaleffects::avg_comparisons(
          combined,
          variables = "Cont",
          newdata = adapter_data
        )

        expect_true(
          nrow(effects) == 2L * nlevels(adapter_data$Sat)
        )
        expect_true(all(is.finite(effects$estimate)))
        expect_true(all(is.finite(effects$std.error)))

        list(
          link = link,
          engines = combined$model_engines,
          effects = effects
        )
      }
    )
  })
}

test_case("Model adapters: clm versus polr", {
  clm_fit <- ordinal::clm(
    Sat ~ Cont + Infl + Type,
    data = adapter_data,
    link = "logit"
  )
  polr_fit <- MASS::polr(
    Sat ~ Cont + Infl + Type,
    data = adapter_data,
    method = "logistic",
    Hess = TRUE
  )

  combined <- suest(
    clm_fit,
    polr_fit,
    model_names = c("clm", "polr")
  )

  expect_true(identical(
    unname(combined$model_types),
    c("ologit", "ologit")
  ))
  expect_true(identical(
    unname(combined$model_engines),
    c("ordinal::clm", "MASS::polr")
  ))

  effects <- marginaleffects::avg_comparisons(
    combined,
    variables = "Cont",
    newdata = adapter_data
  )

  expect_true(all(is.finite(effects$std.error)))

  list(
    engines = combined$model_engines,
    effects = effects
  )
})

test_case("Model adapters: clm disjoint sex subsets", {
  men <- ordinal::clm(
    Sat ~ Cont + Infl + Type,
    data = adapter_data,
    subset = woman == 0,
    link = "logit"
  )
  women <- ordinal::clm(
    Sat ~ Cont + Infl + Type,
    data = adapter_data,
    subset = woman == 1,
    link = "logit"
  )

  combined <- suest(
    men,
    women,
    model_names = c("Men", "Women")
  )

  expect_true(combined$nobs_overlap == 0L)
  expect_true(all(offdiag_vcov(combined) == 0))

  effects <- marginaleffects::avg_comparisons(
    combined,
    variables = "Cont",
    newdata = suest_newdata(combined)
  )

  expect_true(all(is.finite(effects$std.error)))

  list(
    overlap = combined$nobs_overlap,
    effects = effects
  )
})

test_case("Model adapters: clm partial sample overlap", {
  early <- ordinal::clm(
    Sat ~ Cont + Infl + Type,
    data = adapter_data,
    subset = id <= 1200,
    link = "logit"
  )
  late <- ordinal::clm(
    Sat ~ Cont + Infl + Type,
    data = adapter_data,
    subset = id >= 482,
    link = "logit"
  )

  combined <- suest(
    early,
    late,
    model_names = c("Early", "Late")
  )

  expect_true(combined$nobs_models[1L] == 1200L)
  expect_true(combined$nobs_models[2L] == 1200L)
  expect_true(combined$nobs_overlap == 719L)
  expect_true(combined$nobs_union == 1681L)
  expect_true(max(abs(offdiag_vcov(combined))) > 0)

  effects <- marginaleffects::avg_comparisons(
    combined,
    variables = "Cont",
    newdata = suest_newdata(combined)
  )

  expect_true(all(is.finite(effects$std.error)))

  list(
    overlap = combined$nobs_overlap,
    union = combined$nobs_union
  )
})

test_case("Model adapters: unsupported clm extensions", {
  ordinary <- ordinal::clm(
    Sat ~ Cont + Infl,
    data = adapter_data,
    link = "logit"
  )
  scale_model <- ordinal::clm(
    Sat ~ Cont,
    scale = ~ Infl,
    data = adapter_data,
    link = "logit"
  )
  nominal_model <- ordinal::clm(
    Sat ~ Cont,
    nominal = ~ Infl,
    data = adapter_data,
    link = "logit"
  )
  threshold_model <- ordinal::clm(
    Sat ~ Cont + Infl,
    data = adapter_data,
    link = "logit",
    threshold = "equidistant"
  )

  scale_error <- expect_error(
    suest(ordinary, scale_model),
    "scale models"
  )
  nominal_error <- expect_error(
    suest(ordinary, nominal_model),
    "nominal-effects"
  )
  threshold_error <- expect_error(
    suest(ordinary, threshold_model),
    "threshold = 'flexible'"
  )

  list(
    scale = scale_error,
    nominal = nominal_error,
    threshold = threshold_error
  )
})

test_case("Model adapters: adjusted-score GLMs are rejected", {
  dat <- mtcars
  dat$am <- factor(dat$am)
  assign("brglm_dat", dat, envir = globalenv())

  ordinary <- stats::glm(
    am ~ wt,
    family = stats::binomial(),
    data = brglm_dat
  )
  adjusted <- stats::glm(
    am ~ wt,
    family = stats::binomial(),
    data = brglm_dat,
    method = brglm2::brglmFit
  )

  error <- expect_error(
    suest(ordinary, adjusted),
    "Bias-reduced|adjusted-score|Firth|penalized"
  )

  error
})
