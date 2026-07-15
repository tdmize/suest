############################################################
# Cross-model comparison universe
############################################################

universe_simulate_data <- function(n = 1000L, seed = 6109L) {
  set.seed(seed)

  id <- seq_len(n)
  woman <- rep(c(0, 1), length.out = n)
  z <- stats::rnorm(n)
  x <- 0.25 * woman + 0.35 * z + stats::rnorm(n)
  x_alt <- 0.70 * x + stats::rnorm(n, sd = 0.75)
  mediator <- 0.65 * x + 0.25 * z + 0.20 * woman +
    stats::rnorm(n)

  linear_predictor <- -0.30 + 0.45 * x + 0.30 * mediator -
    0.20 * z + 0.15 * woman
  linear_predictor_2 <- 0.10 + 0.25 * x + 0.20 * mediator +
    0.15 * z - 0.10 * woman

  y_lm <- 1 + linear_predictor + stats::rnorm(n, sd = 1.25)
  y_lm_2 <- -0.25 + linear_predictor_2 +
    stats::rnorm(n, sd = 1.10)

  y_bin <- stats::rbinom(
    n,
    1,
    stats::plogis(linear_predictor)
  )
  y_bin_2 <- stats::rbinom(
    n,
    1,
    stats::plogis(linear_predictor_2)
  )

  mu_pois <- exp(0.35 + 0.18 * x + 0.12 * mediator -
    0.08 * z + 0.08 * woman)
  mu_pois_2 <- exp(0.20 - 0.10 * x + 0.16 * mediator +
    0.10 * z + 0.05 * woman)

  y_pois <- stats::rpois(n, mu_pois)
  y_pois_2 <- stats::rpois(n, mu_pois_2)
  y_nb <- stats::rnbinom(n, mu = mu_pois, size = 2.2)
  y_nb_2 <- stats::rnbinom(n, mu = mu_pois_2, size = 1.8)

  latent_ord <- linear_predictor + stats::rlogis(n)
  latent_ord_2 <- linear_predictor_2 + stats::rlogis(n)
  cuts <- c(-Inf, -0.85, 0, 0.85, Inf)
  y_ord <- ordered(
    cut(latent_ord, breaks = cuts, labels = FALSE),
    levels = seq_len(4L),
    labels = c("Low", "Medium", "High", "Very high")
  )
  y_ord_2 <- ordered(
    cut(latent_ord_2, breaks = cuts, labels = FALSE),
    levels = seq_len(4L),
    labels = c("Low", "Medium", "High", "Very high")
  )

  eta_nom <- cbind(
    0,
    -0.20 + 0.30 * x + 0.18 * mediator - 0.10 * z,
    0.15 - 0.20 * x + 0.25 * mediator + 0.12 * woman,
    -0.10 + 0.12 * x - 0.15 * mediator + 0.20 * z
  )
  eta_nom <- eta_nom - apply(eta_nom, 1L, max)
  probability_nom <- exp(eta_nom)
  probability_nom <- probability_nom / rowSums(probability_nom)

  eta_nom_2 <- cbind(
    0,
    0.10 + 0.20 * x + 0.10 * mediator + 0.10 * z,
    -0.15 - 0.12 * x + 0.22 * mediator - 0.08 * woman,
    0.20 + 0.08 * x - 0.10 * mediator + 0.16 * z
  )
  eta_nom_2 <- eta_nom_2 - apply(eta_nom_2, 1L, max)
  probability_nom_2 <- exp(eta_nom_2)
  probability_nom_2 <- probability_nom_2 /
    rowSums(probability_nom_2)

  draw_category <- function(probability) {
    vapply(
      seq_len(nrow(probability)),
      function(i) sample.int(
        ncol(probability),
        size = 1L,
        prob = probability[i, ]
      ),
      integer(1)
    )
  }

  y_nom <- factor(
    draw_category(probability_nom),
    levels = seq_len(4L),
    labels = c("A", "B", "C", "D")
  )
  y_nom_2 <- factor(
    draw_category(probability_nom_2),
    levels = seq_len(4L),
    labels = c("A", "B", "C", "D")
  )

  mediator_missing <- mediator
  mediator_missing[id %% 5L == 0L] <- NA_real_

  data.frame(
    id = id,
    woman = woman,
    x = x,
    x_alt = x_alt,
    mediator = mediator,
    mediator_missing = mediator_missing,
    z = z,
    y_lm = y_lm,
    y_lm_2 = y_lm_2,
    y_bin = y_bin,
    y_bin_2 = y_bin_2,
    y_pois = y_pois,
    y_pois_2 = y_pois_2,
    y_nb = y_nb,
    y_nb_2 = y_nb_2,
    y_ord = y_ord,
    y_ord_2 = y_ord_2,
    y_nom = y_nom,
    y_nom_2 = y_nom_2
  )
}

universe_specs <- list(
  lm = list(
    outcome = "y_lm",
    outcome2 = "y_lm_2",
    fit = "lm"
  ),
  logit = list(
    outcome = "y_bin",
    outcome2 = "y_bin_2",
    fit = "logit"
  ),
  probit = list(
    outcome = "y_bin",
    outcome2 = "y_bin_2",
    fit = "probit"
  ),
  poisson = list(
    outcome = "y_pois",
    outcome2 = "y_pois_2",
    fit = "poisson"
  ),
  negbin = list(
    outcome = "y_nb",
    outcome2 = "y_nb_2",
    fit = "negbin"
  ),
  ologit = list(
    outcome = "y_ord",
    outcome2 = "y_ord_2",
    fit = "ologit"
  ),
  oprobit = list(
    outcome = "y_ord",
    outcome2 = "y_ord_2",
    fit = "oprobit"
  ),
  multinom = list(
    outcome = "y_nom",
    outcome2 = "y_nom_2",
    fit = "multinom"
  )
)

universe_formula <- function(outcome, rhs) {
  formula <- stats::as.formula(
    paste(outcome, "~", rhs),
    env = globalenv()
  )
  environment(formula) <- globalenv()
  formula
}

universe_fit <- function(type, outcome, rhs,
                         subset = NULL,
                         data_name = "universe_data") {
  formula <- universe_formula(outcome, rhs)
  arguments <- list(
    formula = formula,
    data = as.name(data_name)
  )

  if (!is.null(subset))
    arguments$subset <- parse(text = subset)[[1L]]

  call <- switch(
    type,
    lm = as.call(c(
      list(quote(stats::lm)),
      arguments
    )),
    logit = as.call(c(
      list(quote(stats::glm)),
      arguments,
      list(family = quote(stats::binomial("logit")))
    )),
    probit = as.call(c(
      list(quote(stats::glm)),
      arguments,
      list(family = quote(stats::binomial("probit")))
    )),
    poisson = as.call(c(
      list(quote(stats::glm)),
      arguments,
      list(family = quote(stats::poisson("log")))
    )),
    negbin = as.call(c(
      list(quote(MASS::glm.nb)),
      arguments
    )),
    ologit = as.call(c(
      list(quote(MASS::polr)),
      arguments,
      list(
        method = "logistic",
        Hess = TRUE
      )
    )),
    oprobit = as.call(c(
      list(quote(MASS::polr)),
      arguments,
      list(
        method = "probit",
        Hess = TRUE
      )
    )),
    multinom = as.call(c(
      list(quote(nnet::multinom)),
      arguments,
      list(
        Hess = TRUE,
        trace = FALSE
      )
    )),
    stop("Unknown model type: ", type, call. = FALSE)
  )

  eval(call, envir = globalenv())
}

universe_category <- function(group, model) {
  prefix <- paste0(model, "::")
  ifelse(
    startsWith(group, prefix),
    substring(group, nchar(prefix) + 1L),
    NA_character_
  )
}

universe_pair_hypothesis <- function(model1, model2) {
  force(model1)
  force(model2)

  function(x) {
    group <- as.character(x$group)

    scalar1 <- which(group == model1)
    scalar2 <- which(group == model2)

    if (length(scalar1) == 1L && length(scalar2) == 1L) {
      return(data.frame(
        term = paste0(model1, " - ", model2),
        estimate = x$estimate[scalar1] - x$estimate[scalar2]
      ))
    }

    category1 <- universe_category(group, model1)
    category2 <- universe_category(group, model2)
    common <- intersect(
      category1[!is.na(category1)],
      category2[!is.na(category2)]
    )

    if (!length(common))
      stop("No matched scalar or category-specific effects were found.")

    output <- lapply(common, function(category) {
      i1 <- which(category1 == category)
      i2 <- which(category2 == category)

      if (length(i1) != 1L || length(i2) != 1L)
        stop("Effects could not be matched uniquely.")

      data.frame(
        term = paste0(
          model1,
          "::",
          category,
          " - ",
          model2,
          "::",
          category
        ),
        estimate = x$estimate[i1] - x$estimate[i2]
      )
    })

    do.call(rbind, output)
  }
}

universe_alternative_hypothesis <- function(model1, model2,
                                             variable1 = "x",
                                             variable2 = "x_alt") {
  force(model1)
  force(model2)
  force(variable1)
  force(variable2)

  function(x) {
    group <- as.character(x$group)
    term <- as.character(x$term)

    rows1 <- term == variable1
    rows2 <- term == variable2
    scalar1 <- which(rows1 & group == model1)
    scalar2 <- which(rows2 & group == model2)

    if (length(scalar1) == 1L && length(scalar2) == 1L) {
      return(data.frame(
        term = paste0(variable1, " - ", variable2),
        estimate = x$estimate[scalar1] - x$estimate[scalar2]
      ))
    }

    category1 <- universe_category(group, model1)
    category2 <- universe_category(group, model2)
    common <- intersect(
      category1[rows1 & !is.na(category1)],
      category2[rows2 & !is.na(category2)]
    )

    if (!length(common))
      stop("No matched alternative-predictor effects were found.")

    output <- lapply(common, function(category) {
      i1 <- which(rows1 & category1 == category)
      i2 <- which(rows2 & category2 == category)

      if (length(i1) != 1L || length(i2) != 1L)
        stop("Alternative-predictor effects could not be matched uniquely.")

      data.frame(
        term = paste0(
          variable1,
          "::",
          category,
          " - ",
          variable2,
          "::",
          category
        ),
        estimate = x$estimate[i1] - x$estimate[i2]
      )
    })

    do.call(rbind, output)
  }
}

universe_effects <- function(object, variables, newdata,
                             model1, model2,
                             alternative = FALSE) {
  effects <- marginaleffects::avg_comparisons(
    object,
    variables = variables,
    newdata = newdata
  )

  hypothesis <- if (alternative) {
    universe_alternative_hypothesis(model1, model2)
  } else {
    universe_pair_hypothesis(model1, model2)
  }

  differences <- marginaleffects::hypotheses(
    effects,
    hypothesis = hypothesis
  )

  expect_true(
    all(is.finite(differences$estimate)),
    "Cross-model estimates must be finite."
  )
  expect_true(
    all(is.finite(differences$std.error)),
    "Cross-model standard errors must be finite."
  )
  expect_true(
    all(differences$std.error >= 0),
    "Cross-model standard errors must be nonnegative."
  )

  list(
    effects = effects,
    differences = differences
  )
}

test_case("Comparison universe: simulate common data", {
  universe_data <- universe_simulate_data()
  assign("universe_data", universe_data, envir = globalenv())

  expect_true(nrow(universe_data) == 1000L)
  expect_true(all(table(universe_data$woman) == 500L))
  expect_true(all(table(universe_data$y_ord) > 0L))
  expect_true(all(table(universe_data$y_nom) > 0L))

  sapply(
    universe_data[
      c("y_lm", "y_bin", "y_pois", "y_nb", "y_ord", "y_nom")
    ],
    function(x) length(unique(x))
  )
})

############################################################
# Same-family comparison archetypes
############################################################

for (family_name in names(universe_specs)) {
  local({
    family_name <- family_name
    specification <- universe_specs[[family_name]]

    test_case(
      paste0(
        "Comparison universe: nested/mediator - ",
        family_name
      ),
      {
        base <- universe_fit(
          specification$fit,
          specification$outcome,
          "x + z + woman"
        )
        mediator <- universe_fit(
          specification$fit,
          specification$outcome,
          "x + mediator + z + woman"
        )

        combined <- suest(
          base,
          mediator,
          model_names = c("Base", "Mediator")
        )

        expect_true(combined$nobs_overlap == nrow(universe_data))
        expect_true(max(abs(offdiag_vcov(combined))) > 0)

        result <- universe_effects(
          combined,
          variables = "x",
          newdata = universe_data,
          model1 = "Base",
          model2 = "Mediator"
        )

        list(
          model_type = family_name,
          comparison = "nested/mediator",
          contrasts = nrow(result$differences)
        )
      }
    )

    test_case(
      paste0(
        "Comparison universe: alternative predictors - ",
        family_name
      ),
      {
        first <- universe_fit(
          specification$fit,
          specification$outcome,
          "x + z + woman"
        )
        second <- universe_fit(
          specification$fit,
          specification$outcome,
          "x_alt + z + woman"
        )

        combined <- suest(
          first,
          second,
          model_names = c("Predictor x", "Predictor x_alt")
        )

        result <- universe_effects(
          combined,
          variables = c("x", "x_alt"),
          newdata = universe_data,
          model1 = "Predictor x",
          model2 = "Predictor x_alt",
          alternative = TRUE
        )

        list(
          model_type = family_name,
          comparison = "alternative predictors",
          contrasts = nrow(result$differences)
        )
      }
    )

    test_case(
      paste0(
        "Comparison universe: different outcomes - ",
        family_name
      ),
      {
        outcome1 <- universe_fit(
          specification$fit,
          specification$outcome,
          "x + mediator + z + woman"
        )
        outcome2 <- universe_fit(
          specification$fit,
          specification$outcome2,
          "x + mediator + z + woman"
        )

        combined <- suest(
          outcome1,
          outcome2,
          model_names = c("Outcome 1", "Outcome 2")
        )

        result <- universe_effects(
          combined,
          variables = "x",
          newdata = universe_data,
          model1 = "Outcome 1",
          model2 = "Outcome 2"
        )

        list(
          model_type = family_name,
          comparison = "different outcomes",
          contrasts = nrow(result$differences)
        )
      }
    )

    test_case(
      paste0(
        "Comparison universe: sex-stratified subset - ",
        family_name
      ),
      {
        men <- universe_fit(
          specification$fit,
          specification$outcome,
          "x + mediator + z",
          subset = "woman == 0"
        )
        women <- universe_fit(
          specification$fit,
          specification$outcome,
          "x + mediator + z",
          subset = "woman == 1"
        )

        combined <- suest(
          men,
          women,
          model_names = c("Men", "Women")
        )

        expect_true(combined$nobs_overlap == 0L)
        expect_true(combined$nobs_union == nrow(universe_data))
        expect_true(max(abs(offdiag_vcov(combined))) == 0)

        result <- universe_effects(
          combined,
          variables = "x",
          newdata = suest_newdata(combined),
          model1 = "Men",
          model2 = "Women"
        )

        list(
          model_type = family_name,
          comparison = "sex-stratified subset",
          contrasts = nrow(result$differences)
        )
      }
    )

    test_case(
      paste0(
        "Comparison universe: missing mediator overlap - ",
        family_name
      ),
      {
        base <- universe_fit(
          specification$fit,
          specification$outcome,
          "x + z + woman"
        )
        mediator <- universe_fit(
          specification$fit,
          specification$outcome,
          "x + mediator_missing + z + woman"
        )

        combined <- suest(
          base,
          mediator,
          model_names = c("Complete base", "Observed mediator")
        )

        expected_observed <- sum(!is.na(universe_data$mediator_missing))
        expect_true(combined$nobs_models[1L] == nrow(universe_data))
        expect_true(combined$nobs_models[2L] == expected_observed)
        expect_true(combined$nobs_overlap == expected_observed)
        expect_true(combined$nobs_union == nrow(universe_data))
        expect_true(max(abs(offdiag_vcov(combined))) > 0)

        result <- universe_effects(
          combined,
          variables = "x",
          newdata = suest_newdata(combined),
          model1 = "Complete base",
          model2 = "Observed mediator"
        )

        list(
          model_type = family_name,
          comparison = "missing mediator overlap",
          overlap = combined$nobs_overlap,
          contrasts = nrow(result$differences)
        )
      }
    )

    test_case(
      paste0(
        "Comparison universe: partially overlapping subsets - ",
        family_name
      ),
      {
        early <- universe_fit(
          specification$fit,
          specification$outcome,
          "x + mediator + z + woman",
          subset = "id <= 750"
        )
        late <- universe_fit(
          specification$fit,
          specification$outcome,
          "x + mediator + z + woman",
          subset = "id >= 251"
        )

        combined <- suest(
          early,
          late,
          model_names = c("Early subset", "Late subset")
        )

        expect_true(combined$nobs_models[1L] == 750L)
        expect_true(combined$nobs_models[2L] == 750L)
        expect_true(combined$nobs_overlap == 500L)
        expect_true(combined$nobs_union == 1000L)
        expect_true(max(abs(offdiag_vcov(combined))) > 0)

        result <- universe_effects(
          combined,
          variables = "x",
          newdata = suest_newdata(combined),
          model1 = "Early subset",
          model2 = "Late subset"
        )

        list(
          model_type = family_name,
          comparison = "partially overlapping subsets",
          overlap = combined$nobs_overlap,
          contrasts = nrow(result$differences)
        )
      }
    )
  })
}

############################################################
# Allowed cross-family pairs
############################################################

cross_family_specs <- list(
  list(
    name = "logit-probit",
    type1 = "logit",
    type2 = "probit",
    outcome1 = "y_bin",
    outcome2 = "y_bin"
  ),
  list(
    name = "logit-linear",
    type1 = "logit",
    type2 = "lm",
    outcome1 = "y_bin",
    outcome2 = "y_bin"
  ),
  list(
    name = "probit-linear",
    type1 = "probit",
    type2 = "lm",
    outcome1 = "y_bin",
    outcome2 = "y_bin"
  ),
  list(
    name = "Poisson-negative binomial",
    type1 = "poisson",
    type2 = "negbin",
    outcome1 = "y_nb",
    outcome2 = "y_nb"
  ),
  list(
    name = "ordered-multinomial",
    type1 = "ologit",
    type2 = "multinom",
    outcome1 = "y_ord",
    outcome2 = "y_ord"
  )
)

# Multinomial needs an unordered response with the same category labels.
universe_data$y_ord_nom <- factor(
  universe_data$y_ord,
  levels = levels(universe_data$y_ord)
)
assign("universe_data", universe_data, envir = globalenv())
cross_family_specs[[5L]]$outcome2 <- "y_ord_nom"

for (specification in cross_family_specs) {
  local({
    specification <- specification

    test_case(
      paste0(
        "Comparison universe: allowed cross-family same sample - ",
        specification$name
      ),
      {
        model1 <- universe_fit(
          specification$type1,
          specification$outcome1,
          "x + mediator + z + woman"
        )
        model2 <- universe_fit(
          specification$type2,
          specification$outcome2,
          "x + mediator + z + woman"
        )

        combined <- suest(
          model1,
          model2,
          model_names = c("Model 1", "Model 2")
        )

        expect_true(combined$nobs_overlap == nrow(universe_data))
        expect_true(max(abs(offdiag_vcov(combined))) > 0)

        result <- universe_effects(
          combined,
          variables = "x",
          newdata = universe_data,
          model1 = "Model 1",
          model2 = "Model 2"
        )

        list(
          comparison = specification$name,
          sample = "same",
          contrasts = nrow(result$differences)
        )
      }
    )

    test_case(
      paste0(
        "Comparison universe: allowed cross-family disjoint samples - ",
        specification$name
      ),
      {
        model1 <- universe_fit(
          specification$type1,
          specification$outcome1,
          "x + mediator + z",
          subset = "woman == 0"
        )
        model2 <- universe_fit(
          specification$type2,
          specification$outcome2,
          "x + mediator + z",
          subset = "woman == 1"
        )

        combined <- suest(
          model1,
          model2,
          model_names = c("Men model", "Women model")
        )

        expect_true(combined$nobs_overlap == 0L)
        expect_true(max(abs(offdiag_vcov(combined))) == 0)

        result <- universe_effects(
          combined,
          variables = "x",
          newdata = suest_newdata(combined),
          model1 = "Men model",
          model2 = "Women model"
        )

        list(
          comparison = specification$name,
          sample = "disjoint",
          contrasts = nrow(result$differences)
        )
      }
    )
  })
}

############################################################
# Data-object selection behavior
############################################################

test_case("Comparison universe: separately filtered data objects are disjoint", {
  universe_men <- universe_data[universe_data$woman == 0, , drop = FALSE]
  universe_women <- universe_data[universe_data$woman == 1, , drop = FALSE]
  assign("universe_men", universe_men, envir = globalenv())
  assign("universe_women", universe_women, envir = globalenv())

  men <- universe_fit(
    "logit",
    "y_bin",
    "x + mediator + z",
    data_name = "universe_men"
  )
  women <- universe_fit(
    "logit",
    "y_bin",
    "x + mediator + z",
    data_name = "universe_women"
  )

  combined <- suest(
    men,
    women,
    model_names = c("Men object", "Women object")
  )

  expect_true(combined$nobs_overlap == 0L)
  expect_true(max(abs(offdiag_vcov(combined))) == 0)

  result <- universe_effects(
    combined,
    variables = "x",
    newdata = suest_newdata(combined),
    model1 = "Men object",
    model2 = "Women object"
  )

  list(
    behavior = "different data objects treated as disjoint",
    contrasts = nrow(result$differences)
  )
})
