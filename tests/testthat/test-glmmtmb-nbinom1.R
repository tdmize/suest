nbinom1_data <- function(n = 600L, groups = 40L) {
  set.seed(3301)
  g <- rep(seq_len(groups), length.out = n)
  x <- rnorm(n)
  z <- rbinom(n, 1L, 0.4)
  u <- rnorm(groups, sd = 0.5)[g]
  mu1 <- exp(0.3 + 0.5*x - 0.4*z + u)
  mu2 <- exp(0.1 + 0.3*x + 0.2*z + 0.6*u)
  data.frame(row = seq_len(n), id = factor(g), x = x, z = z,
    y1 = rnbinom(n, size = mu1/1.5, mu = mu1),
    y2 = rnbinom(n, size = 2, mu = mu2))
}

nbinom1_fit <- function(formula, data, ...)
  glmmTMB::glmmTMB(formula, data = data, family = glmmTMB::nbinom1(), ...)

test_that("constant-dispersion negative binomial scores are exact", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  skip_if_not_installed("numDeriv")
  d <- nbinom1_data()
  model <- nbinom1_fit(y1 ~ x + z, d)
  adapter <- suest:::.suest_model_adapter(model)
  expect_identical(adapter$type, "nb1")
  component <- suest:::.suest_model_components(model, adapter$type, adapter$engine)
  parameters <- component$parameters
  expect_identical(names(parameters), c("(Intercept)", "x", "z", "log_phi"))
  expect_equal(unname(parameters["log_phi"]), log(sigma(model)), tolerance = 1e-12)
  expect_equal(unname(component$bread/nrow(d)),
    unname(as.matrix(vcov(model, full = TRUE))), tolerance = 1e-12)
  expect_lt(max(abs(colSums(component$score))), 1e-3)

  # Variance mu*(1 + phi): size mu/phi in R's parameterization.
  X <- model.matrix(~ x + z, d)
  density <- function(theta, i) {
    mu <- exp(sum(X[i, ]*theta[1:3]))
    dnbinom(d$y1[i], size = mu/exp(theta[4]), mu = mu, log = TRUE)
  }
  for (i in c(1L, 77L, 400L)) {
    expect_equal(unname(component$score[i, ]),
      numDeriv::grad(density, unname(parameters), i = i), tolerance = 1e-7)
  }
  information <- -numDeriv::hessian(function(theta)
    sum(vapply(seq_len(nrow(d)), density, numeric(1), theta = theta)), unname(parameters))
  expect_equal(unname(component$bread), unname(nrow(d)*solve(information)), tolerance = 1e-6)

  if (requireNamespace("sandwich", quietly = TRUE)) {
    theirs <- sandwich::estfun(model, full = TRUE, cluster = factor(seq_len(nrow(d))))
    expect_equal(unname(component$score), unname(as.matrix(theirs)), tolerance = 1e-7)
  }
})

test_that("constant- and mean-dispersion models combine with a robust system", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  skip_if_not_installed("MASS")
  d <- nbinom1_data()
  constant <- nbinom1_fit(y1 ~ x + z, d)
  mean_dispersion <- MASS::glm.nb(y1 ~ x + z, data = d)
  fit <- suest(constant, mean_dispersion, model_names = c("Constant", "Mean"),
    observation_id = "row")
  expect_identical(unname(fit$model_types), c("nb1", "negbin"))
  expect_identical(fit$comparison_scale, "expected counts")
  expect_identical(names(coef(fit))[1:4],
    c("Constant::(Intercept)", "Constant::x", "Constant::z", "Constant::log_phi"))

  component <- suest:::.suest_model_components(constant, "nb1", "glmmTMB::glmmTMB")
  n <- nrow(d)
  B <- component$bread
  expected <- B %*% crossprod(component$score) %*% B/n^2*n/(n - 1)
  index <- fit$index[["Constant"]]
  expect_equal(unname(vcov(fit)[index, index]), unname(expected), tolerance = 1e-10)
  expect_equal(unname(coef(fit)[index]), unname(component$parameters), tolerance = 1e-12)
})

test_that("constant-dispersion predictions are expected counts", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  skip_if_not_installed("marginaleffects")
  d <- nbinom1_data()
  first <- nbinom1_fit(y1 ~ x*z, d)
  second <- nbinom1_fit(y2 ~ x + z, d)
  expect_equal(
    suest:::.suest_predict_values(first, d[1:6, ], "response", "glmmTMB::glmmTMB", "nb1"),
    unname(predict(first, newdata = d[1:6, ], type = "response")), tolerance = 1e-10)
  expect_equal(
    suest:::.suest_predict_values(first, d[1:6, ], "link", "glmmTMB::glmmTMB", "nb1"),
    unname(predict(first, newdata = d[1:6, ], type = "link")), tolerance = 1e-10)

  fit <- suest(first, second, model_names = c("First", "Second"), observation_id = "row")
  averages <- marginaleffects::avg_predictions(fit, newdata = d)
  expect_equal(averages$estimate, c(mean(predict(first, type = "response")),
    mean(predict(second, type = "response"))), tolerance = 1e-8)
  expect_true(all(is.finite(averages$std.error)))

  changed <- coef(fit)
  changed["First::x"] <- changed["First::x"] + 0.1
  altered <- marginaleffects::set_coef(fit, changed)
  before <- marginaleffects::predictions(fit, newdata = d[1:3, ], vcov = FALSE)
  after <- marginaleffects::predictions(altered, newdata = d[1:3, ], vcov = FALSE)
  rows <- before$group == "First"
  expect_equal(after$estimate[rows], before$estimate[rows]*exp(0.1*d$x[1:3]),
    tolerance = 1e-10)
})

test_that("constant-dispersion random effects use the multilevel engine", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  d <- nbinom1_data()
  intercept <- suppressWarnings(nbinom1_fit(y1 ~ x + z + (1 | id), d))
  expect_identical(suest:::.suest_model_adapter(intercept)$type, "glmm_nbinom1_ri")
  component <- suest:::.suest_model_components(intercept, "glmm_nbinom1_ri", "glmmTMB::glmmTMB")
  expect_identical(names(component$parameters),
    c("(Intercept)", "x", "z", "log_phi", "log_sigma"))
  expect_equal(unname(component$bread/nrow(d)),
    unname(as.matrix(vcov(intercept, full = TRUE))), tolerance = 1e-10)
  fit <- suest(intercept, nbinom1_fit(y1 ~ x + z, d), model_names = c("Multilevel", "Single"),
    observation_id = "row")
  expect_identical(fit$n_clusters, 40L)
  expect_identical(unname(fit$model_types), c("glmm_nbinom1_ri", "nb1"))
})

test_that("unsupported single-level glmmTMB models are refused with a reason", {
  skip_if_not_installed("glmmTMB", minimum_version = "1.1.14")
  d <- nbinom1_data()
  expect_error(suest(glmmTMB::glmmTMB(y1 ~ x, d, family = poisson), glm(y1 ~ x, poisson, d)),
    "constant-dispersion negative binomial")
  expect_error(suest(glmmTMB::glmmTMB(y1 ~ x, d, family = glmmTMB::nbinom2()),
    glm(y1 ~ x, poisson, d)), "glm.nb")
  expect_error(suest(nbinom1_fit(y1 ~ x, d, dispformula = ~ z), glm(y1 ~ x, poisson, d)),
    "constant dispersion")
  expect_error(suest(suppressWarnings(nbinom1_fit(y1 ~ x, d, ziformula = ~ 1)),
    glm(y1 ~ x, poisson, d)), "Zero-inflation")
  expect_error(suest(nbinom1_fit(y1 ~ x + offset(z), d), glm(y1 ~ x, poisson, d)),
    "Offsets")
})
