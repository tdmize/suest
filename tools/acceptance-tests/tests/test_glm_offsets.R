cat("\n\nADDITIONAL GLM AND OFFSET TESTS\n")

acceptance_loglog_link <- function() {
  structure(
    list(
      linkfun = function(mu) -log(-log(mu)),
      linkinv = function(eta) exp(-exp(-eta)),
      mu.eta = function(eta) exp(-eta - exp(-eta)),
      valideta = function(eta) TRUE,
      name = "loglog"
    ),
    class = "link-glm"
  )
}

test_case("GLM families: Gaussian, Gamma, cloglog, and fractional", {
  set.seed(8135)
  n <- 1200L
  glm_data <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  eta <- -0.2 + 0.45*glm_data$x - 0.25*glm_data$z
  probability <- 1 - exp(-exp(eta))
  glm_data$y_gaussian <- 1 + eta + stats::rnorm(n)
  glm_data$y_gamma <- stats::rgamma(
    n,
    shape = 3,
    scale = exp(0.4 + 0.15*eta)/3
  )
  glm_data$y_binary <- stats::rbinom(n, 1, probability)
  glm_data$y_fractional <- stats::plogis(
    eta + stats::rnorm(n, sd = 0.8)
  )
  gaussian <- stats::glm(
    y_gaussian ~ x + z,
    family = stats::gaussian("identity"),
    data = glm_data
  )
  gamma <- stats::glm(
    y_gamma ~ x + z,
    family = stats::Gamma("log"),
    data = glm_data
  )
  cloglog <- stats::glm(
    y_binary ~ x + z,
    family = stats::binomial("cloglog"),
    data = glm_data
  )
  fractional <- stats::glm(
    y_fractional ~ x + z,
    family = stats::quasibinomial("probit"),
    data = glm_data
  )
  fractional_loglog <- stats::glm(
    y_fractional ~ x + z,
    family = stats::quasibinomial(acceptance_loglog_link()),
    data = glm_data
  )
  combined <- suest(
    gaussian,
    gamma,
    cloglog,
    fractional,
    fractional_loglog,
    model_names = c(
      "Gaussian", "Gamma", "Cloglog", "Fractional", "Frac log-log"
    )
  )
  effects <- marginaleffects::avg_comparisons(
    combined,
    variables = "x",
    newdata = glm_data
  )
  differences <- marginaleffects::hypotheses(
    effects,
    hypothesis = difference ~ pairwise
  )

  expect_true(nrow(effects) == 5L)
  expect_true(nrow(differences) == 10L)
  expect_true(all(is.finite(effects$std.error)))
  expect_true(all(is.finite(differences$std.error)))

  list(object = combined, effects = effects, differences = differences)
})

test_case("Offsets: scalar and ordered marginaleffects", {
  set.seed(4802)
  n <- 900L
  offset_data <- data.frame(
    x = stats::rnorm(n),
    z = stats::rnorm(n),
    offset = 0.3*stats::rnorm(n)
  )
  eta <- -0.2 + 0.5*offset_data$x - 0.25*offset_data$z +
    offset_data$offset
  offset_data$exposure <- exp(offset_data$offset)
  offset_data$y_lm <- eta + stats::rnorm(n)
  offset_data$y_count <- stats::rpois(
    n,
    exp(0.4 + 0.2*offset_data$x - 0.1*offset_data$z + offset_data$offset)
  )
  offset_data$y_ord <- ordered(
    cut(
      eta + stats::rlogis(n),
      c(-Inf, -0.8, 0, 0.8, Inf),
      labels = FALSE
    ),
    levels = 1:4,
    labels = c("A", "B", "C", "D")
  )
  linear <- stats::lm(
    y_lm ~ x + z + offset(offset),
    data = offset_data
  )
  poisson <- stats::glm(
    y_count ~ x + z + offset(log(exposure)),
    family = stats::poisson(),
    data = offset_data
  )
  scalar <- suest(linear, poisson, model_names = c("Linear", "Poisson"))
  scalar_effects <- marginaleffects::avg_comparisons(
    scalar,
    variables = "x",
    newdata = offset_data
  )
  polr <- MASS::polr(
    y_ord ~ x + z + offset(offset),
    data = offset_data,
    method = "logistic",
    Hess = TRUE
  )
  clm <- ordinal::clm(
    y_ord ~ x + z + offset(offset),
    data = offset_data,
    link = "logit",
    Hess = TRUE,
    model = TRUE
  )
  categorical <- suest(polr, clm, model_names = c("polr", "clm"))
  categorical_effects <- marginaleffects::avg_comparisons(
    categorical,
    variables = "x",
    newdata = offset_data
  )

  expect_true(all(is.finite(scalar_effects$std.error)))
  expect_true(all(is.finite(categorical_effects$std.error)))

  list(scalar = scalar_effects, categorical = categorical_effects)
})
