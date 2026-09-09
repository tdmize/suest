panel_ml_test_data <- function(groups = 60L, periods = 4L) {
  set.seed(9472)
  id <- rep(seq_len(groups), each = periods)
  time <- rep(seq_len(periods), groups)
  higher <- ceiling(id / 3)
  x <- rnorm(length(id))
  z <- rnorm(length(id))
  random1 <- rnorm(groups, sd = 0.9)[id]
  random2 <- rnorm(groups, sd = 0.7)[id]
  y1 <- 0.8 + 0.6*x - 0.25*z + 0.08*time + random1 +
    rnorm(length(id), sd = 0.65)
  y2 <- -0.4 + 0.3*x + 0.4*z - 0.04*time + random2 +
    0.25*random1 + rnorm(length(id), sd = 0.75)
  data.frame(id = factor(id), time = factor(time), higher, x, z, y1, y2)
}

panel_ml_loglik <- function(parameters, model, frame) {
  beta_names <- names(model$coefficients$fixed)
  beta <- parameters[beta_names]
  sigma_u <- unname(parameters["sigma_u"])
  sigma_e <- unname(parameters["sigma_e"])
  if (sigma_u <= 0 || sigma_e <= 0)
    return(-Inf)

  X <- model.matrix(
    delete.response(terms(model)), frame,
    contrasts.arg = model$contrasts
  )[, beta_names, drop = FALSE]
  residual <- model.response(frame) - as.numeric(X %*% beta)
  panel <- as.character(model$groups[[1L]])
  value <- 0
  for (rows in split(seq_along(panel), panel)) {
    r <- residual[rows]
    Ti <- length(rows)
    sigma_e2 <- sigma_e^2
    sigma_u2 <- sigma_u^2
    denominator <- sigma_e2 + Ti*sigma_u2
    quadratic <- sum(r^2)/sigma_e2 -
      sigma_u2*sum(r)^2/(sigma_e2*denominator)
    value <- value - 0.5*(
      Ti*log(2*pi) + (Ti - 1)*log(sigma_e2) +
        log(denominator) + quadratic
    )
  }
  value
}

test_that("nlme random-intercept ML scores and bread match its likelihood", {
  skip_if_not_installed("nlme")
  dat <- panel_ml_test_data()
  model <- nlme::lme(
    y1 ~ x + z + time, random = ~1 | id, data = dat,
    method = "ML", na.action = na.omit
  )
  component <- suest:::.suest_lme_ml_components(model)
  frame <- suest:::.suest_model_frame(model, "nlme::lme")
  parameters <- component$parameters

  expect_equal(unname(parameters["sigma_u"]),
               sqrt(as.numeric(nlme::getVarCov(model)[1L, 1L])),
               tolerance = 1e-14)
  expect_equal(panel_ml_loglik(parameters, model, frame),
               as.numeric(logLik(model)), tolerance = 1e-10)

  expect_lt(max(abs(colSums(component$score))), 1e-4)

  likelihood <- function(x) -panel_ml_loglik(x, model, frame)
  numerical_information <- optimHess(parameters, likelihood)
  numerical_covariance <- solve(numerical_information)
  expect_lt(
    max(abs(component$bread/nrow(frame) - numerical_covariance)),
    5e-5
  )

  first_panel <- which(as.character(model$groups[[1L]]) ==
                         as.character(model$groups[[1L]][1L]))
  panel_likelihood <- function(x) {
    panel_frame <- frame[first_panel, , drop = FALSE]
    panel_model <- model
    panel_model$groups <- model$groups[first_panel, , drop = FALSE]
    panel_ml_loglik(x, panel_model, panel_frame)
  }
  step <- 1e-6*pmax(1, abs(parameters))
  numerical_score <- vapply(seq_along(parameters), function(j) {
    plus <- minus <- parameters
    plus[j] <- plus[j] + step[j]
    minus[j] <- minus[j] - step[j]
    (panel_likelihood(plus) - panel_likelihood(minus))/(2*step[j])
  }, numeric(1))
  names(numerical_score) <- names(parameters)
  expect_equal(
    colSums(component$score[first_panel, , drop = FALSE]),
    numerical_score,
    tolerance = 2e-6
  )
})

test_that("nlme random-intercept ML systems include both variance components", {
  skip_if_not_installed("nlme")
  dat <- panel_ml_test_data()
  first <- nlme::lme(
    y1 ~ x + z + time, random = ~1 | id, data = dat,
    method = "ML", na.action = na.omit
  )
  second <- nlme::lme(
    y2 ~ x + z + time, random = ~1 | id, data = dat,
    method = "ML", na.action = na.omit
  )
  fit <- suest(first, second, model_names = c("First", "Second"))

  expect_equal(fit$n_clusters, nlevels(dat$id))
  expect_true(all(c(
    "First::sigma_u", "First::sigma_e",
    "Second::sigma_u", "Second::sigma_e"
  ) %in% names(coef(fit))))
  expect_true(all(is.finite(vcov(fit))))
  expect_equal(vcov(fit), t(vcov(fit)), tolerance = 1e-12)

  effects <- marginaleffects::avg_slopes(
    fit, variables = "x", newdata = dat
  )
  expect_equal(
    effects$estimate,
    c(unname(nlme::fixef(first)["x"]), unname(nlme::fixef(second)["x"])),
    tolerance = 1e-8
  )
  expect_true(all(is.finite(effects$std.error)))
})

test_that("nlme random-intercept ML supports overlapping samples and higher clusters", {
  skip_if_not_installed("nlme")
  dat <- panel_ml_test_data()
  dat$y1[dat$time == "4" & as.integer(dat$id) %% 4 == 0] <- NA
  dat$y2[dat$time == "1" & as.integer(dat$id) %% 5 == 0] <- NA
  first <- nlme::lme(
    y1 ~ x + z, random = ~1 | id, data = dat,
    method = "ML", na.action = na.omit
  )
  second <- nlme::lme(
    y2 ~ x + z, random = ~1 | id, data = dat,
    method = "ML", na.action = na.omit
  )
  fit <- suest(
    first, second, model_names = c("First", "Second"), cluster = "higher"
  )

  expect_equal(fit$n_clusters, length(unique(dat$higher)))
  expect_lt(fit$nobs_overlap, fit$nobs_union)
  expect_true(all(is.finite(vcov(fit))))
  expect_error(
    suest(first, second, cluster = "time"),
    "Panel IDs must be nested within clusters"
  )
})

test_that("nlme panel ML validation rejects different likelihoods and covariance models", {
  skip_if_not_installed("nlme")
  dat <- panel_ml_test_data(groups = 35L)
  reml <- nlme::lme(y1 ~ x + z, random = ~1 | id, data = dat)
  correlated <- nlme::lme(
    y1 ~ x + z, random = ~1 | id,
    correlation = nlme::corAR1(form = ~as.integer(time) | id),
    data = dat, method = "ML"
  )

  expect_error(
    suest(reml, reml),
    "method = \"ML\"",
    fixed = TRUE
  )
  expect_error(
    suest(correlated, correlated),
    "Correlated or heteroskedastic"
  )
  slope_only <- nlme::lme(
    y1 ~ x + z, random = ~0 + x | id, data = dat, method = "ML"
  )
  expect_error(suest(slope_only, slope_only), "random-intercept")
})
