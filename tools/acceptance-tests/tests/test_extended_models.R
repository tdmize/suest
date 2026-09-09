cat("\n\nEXTENDED MODEL-FAMILY TESTS\n")

extended_robust_vcov <- function(model) {
  U <- sandwich::estfun(model)
  B <- sandwich::bread(model)
  n <- nrow(U)
  B %*% crossprod(U) %*% B / n^2 * n/(n - 1)
}

set.seed(9628)
extended_n <- 850L
extended_data <- data.frame(
  id = seq_len(extended_n),
  x = stats::rnorm(extended_n),
  z = stats::rnorm(extended_n)
)
extended_eta <- 0.3 + 0.65*extended_data$x - 0.25*extended_data$z
extended_data$latent <- extended_eta + stats::rnorm(extended_n, sd = 0.85)
extended_data$censored <- pmax(0, extended_data$latent)
extended_data$lower <- floor(2*extended_data$latent)/2
extended_data$upper <- extended_data$lower + 0.5
extended_data$beta_y <- stats::plogis(
  extended_eta + stats::rnorm(extended_n, sd = 0.75)
)
extended_data$count <- stats::rpois(
  extended_n,
  exp(0.15 + 0.3*extended_data$x - 0.2*extended_data$z)
)
extra_zero <- stats::rbinom(
  extended_n,
  1,
  stats::plogis(-0.7 + 0.35*extended_data$z)
)
extended_data$zero_count <- ifelse(extra_zero == 1L, 0L, extended_data$count)
extended_data$hetero_binary <- stats::rbinom(
  extended_n,
  1,
  stats::pnorm(
    (0.2 + 0.6*extended_data$x - 0.2*extended_data$z)/
      exp(0.3*extended_data$z)
  )
)

test_case("Censored regression: censReg and Gaussian survreg equations", {
  direct <- censReg::censReg(
    censored ~ x + z,
    left = 0,
    data = extended_data
  )
  survival <- survival::survreg(
    survival::Surv(censored, latent > 0, type = "left") ~ x + z,
    data = extended_data,
    dist = "gaussian",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  interval <- survival::survreg(
    survival::Surv(lower, upper, type = "interval2") ~ x + z,
    data = extended_data,
    dist = "gaussian",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  combined <- suest(direct, survival, interval)

  expect_near(
    stats::coef(direct),
    c(stats::coef(survival), log(survival$scale)),
    tolerance = 1e-6,
    label = "Tobit coefficients"
  )
  expect_near(
    sandwich::estfun(direct),
    sandwich::estfun(survival),
    tolerance = 1e-5,
    label = "Tobit scores"
  )
  expect_true(all(is.finite(stats::vcov(combined))))
  effects <- marginaleffects::avg_comparisons(
    combined,
    variables = "x",
    newdata = extended_data
  )
  expect_true(nrow(effects) == 3L)
  expect_true(all(is.finite(effects$std.error)))

  list(object = combined, effects = effects)
})

test_case("Beta regression: links and modeled precision", {
  constant <- betareg::betareg(
    beta_y ~ x,
    data = extended_data,
    link = "cloglog"
  )
  precision <- betareg::betareg(
    beta_y ~ x + z | z,
    data = extended_data,
    link = "loglog"
  )
  combined <- suest(constant, precision)
  effects <- marginaleffects::avg_comparisons(
    combined,
    variables = "x",
    newdata = extended_data
  )

  expect_near(
    stats::vcov(combined)[combined$index[[1L]], combined$index[[1L]]],
    extended_robust_vcov(constant),
    tolerance = 1e-8,
    label = "beta robust covariance"
  )
  expect_true(any(grepl("\\(phi\\)", names(stats::coef(combined)))))
  expect_true(nrow(effects) == 2L)
  expect_true(all(is.finite(effects$std.error)))

  list(object = combined, effects = effects)
})

test_case("Zero-inflated Poisson and negative binomial", {
  zip <- pscl::zeroinfl(
    zero_count ~ x + z | z,
    data = extended_data,
    dist = "poisson",
    x = TRUE,
    y = TRUE
  )
  zinb <- pscl::zeroinfl(
    zero_count ~ x | z,
    data = extended_data,
    dist = "negbin",
    x = TRUE,
    y = TRUE
  )
  combined <- suest(zip, zinb)
  effects <- marginaleffects::avg_comparisons(
    combined,
    variables = "x",
    newdata = extended_data
  )

  expect_near(
    stats::vcov(combined)[combined$index[[1L]], combined$index[[1L]]],
    extended_robust_vcov(zip),
    tolerance = 1e-8,
    label = "ZIP robust covariance"
  )
  expect_true("zinb::ln_theta" %in% names(stats::coef(combined)))
  expect_true(nrow(effects) == 2L)
  expect_true(all(is.finite(effects$std.error)))

  list(object = combined, effects = effects)
})

test_case("Truncated Gaussian regression: partial samples", {
  left_data <- extended_data[extended_data$latent > 0, ]
  right_data <- extended_data[extended_data$latent < 1, ]
  left <- truncreg::truncreg(
    latent ~ x + z,
    data = left_data,
    point = 0,
    direction = "left",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  right <- truncreg::truncreg(
    latent ~ x + z,
    data = right_data,
    point = 1,
    direction = "right",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  combined <- suest(left, right, observation_id = "id")
  effects <- marginaleffects::avg_comparisons(
    combined,
    variables = "x",
    newdata = extended_data
  )

  expect_true(
    combined$nobs_overlap == length(intersect(left_data$id, right_data$id))
  )
  expect_near(
    stats::vcov(combined)[combined$index[[1L]], combined$index[[1L]]],
    extended_robust_vcov(left),
    tolerance = 1e-8,
    label = "truncated robust covariance"
  )
  expect_true(nrow(effects) == 2L)
  expect_true(all(is.finite(effects$std.error)))

  list(object = combined, effects = effects)
})

test_case("Instrumental variables: 2SLS equations and predictions", {
  set.seed(7735)
  n <- 1000L
  iv_data <- data.frame(
    id = seq_len(n),
    x = stats::rnorm(n),
    z1 = stats::rnorm(n),
    z2 = stats::rnorm(n),
    disturbance = stats::rnorm(n)
  )
  iv_data$endogenous <- 0.7*iv_data$z1 + 0.25*iv_data$z2 +
    0.4*iv_data$disturbance + stats::rnorm(n)
  iv_data$y <- 0.4 + 0.65*iv_data$x + 1.1*iv_data$endogenous +
    iv_data$disturbance
  base <- fixest::feols(y ~ x | endogenous ~ z1, data = iv_data)
  adjusted <- fixest::feols(
    y ~ x | endogenous ~ z1 + z2,
    data = iv_data
  )
  combined <- suest(base, adjusted)

  X <- cbind(1, iv_data$endogenous, iv_data$x)
  Z <- cbind(1, iv_data$z1, iv_data$z2, iv_data$x)
  Xhat <- qr.fitted(qr(Z), X)
  beta <- solve(crossprod(Xhat, X), crossprod(Xhat, iv_data$y))
  score <- Xhat*as.numeric(iv_data$y - X %*% beta)
  bread <- n*solve(crossprod(Xhat, X))
  effects <- marginaleffects::avg_comparisons(
    combined,
    variables = "x",
    newdata = iv_data
  )

  expect_near(
    sandwich::estfun(adjusted),
    score,
    tolerance = 1e-9,
    label = "2SLS scores"
  )
  expect_near(
    sandwich::bread(adjusted),
    bread,
    tolerance = 1e-9,
    label = "2SLS bread"
  )
  expect_true(nrow(effects) == 2L)
  expect_true(all(is.finite(effects$std.error)))

  list(object = combined, effects = effects)
})

test_case("Heteroskedastic probit: scale scores and predictions", {
  base <- Rchoice::hetprob(
    hetero_binary ~ x | z,
    data = extended_data
  )
  adjusted <- Rchoice::hetprob(
    hetero_binary ~ x + z | z,
    data = extended_data,
    link = "probit"
  )
  combined <- suest(base, adjusted)
  effects <- marginaleffects::avg_comparisons(
    combined,
    variables = "x",
    newdata = extended_data
  )

  expect_near(
    stats::vcov(combined)[combined$index[[1L]], combined$index[[1L]]],
    extended_robust_vcov(base),
    tolerance = 1e-8,
    label = "heteroskedastic-probit robust covariance"
  )
  expect_true(any(grepl("het\\.", names(stats::coef(combined)))))
  expect_true(nrow(effects) == 2L)
  expect_true(all(is.finite(effects$std.error)))

  list(object = combined, effects = effects)
})
