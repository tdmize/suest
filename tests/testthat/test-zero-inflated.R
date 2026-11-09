zeroinfl_test_data <- function(n = 1000L) {
  set.seed(7642)
  x <- rnorm(n)
  z <- rnorm(n)
  structural_zero <- rbinom(n, 1, plogis(-1 + 0.45*z))
  count <- rnbinom(n, mu = exp(0.4 + 0.35*x - 0.2*z), size = 2.8)

  data.frame(
    id = seq_len(n),
    x = x,
    z = z,
    y = ifelse(structural_zero == 1, 0, count)
  )
}

# Robust covariance with the exact observed information, from the sandwich
# package's own analytic zeroinfl scores differentiated by Richardson
# extrapolation. (A numerical Hessian of the log likelihood is too noisy for
# probit inflation on some platforms: 1e-6 on macOS.)
zip_exact_vcov <- function(model) {
  k <- length(model$coefficients$count)
  at <- function(p) {
    model$coefficients$count[] <- p[seq_len(k)]
    model$coefficients$zero[] <- p[-seq_len(k)]
    model
  }
  p <- stats::coef(model)
  U <- sandwich::estfun(model)
  H <- numDeriv::jacobian(function(q) colSums(sandwich::estfun(at(q))), p)
  B <- solve(-(H + t(H))/2)
  n <- nrow(U)
  B %*% crossprod(U) %*% B * n/(n - 1)
}

test_that("zero-inflated Poisson uses the exact observed information", {
  skip_if_not_installed("pscl")
  skip_if_not_installed("numDeriv")
  dat <- zeroinfl_test_data()
  base <- pscl::zeroinfl(
    y ~ x | z,
    data = dat,
    dist = "poisson",
    link = "logit",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  adjusted <- pscl::zeroinfl(
    y ~ x + z | x + z,
    data = dat,
    dist = "poisson",
    link = "probit",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  fit <- suest(base, adjusted)

  expect_equal(unname(fit$model_types), rep("zip", 2L))
  components <- .suest_zeroinfl_components(base, "zip")
  expect_equal(
    unname(components$score),
    unname(sandwich::estfun(base)),
    tolerance = 1e-10
  )
  expect_equal(
    unname(vcov(fit)[fit$index[[1L]], fit$index[[1L]], drop = FALSE]),
    unname(zip_exact_vcov(base)),
    tolerance = 1e-7
  )
  expect_equal(
    unname(vcov(fit)[fit$index[[2L]], fit$index[[2L]], drop = FALSE]),
    unname(zip_exact_vcov(adjusted)),
    tolerance = 1e-7
  )
})

test_that("ZINB includes dispersion and matches numerical likelihood scores", {
  skip_if_not_installed("pscl")
  dat <- zeroinfl_test_data()
  model <- pscl::zeroinfl(
    y ~ x + z | z,
    data = dat,
    dist = "negbin",
    link = "logit",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  components <- .suest_zeroinfl_components(model, "zinb")
  parameters <- components$parameters

  expect_true("ln_theta" %in% names(parameters))
  expect_equal(
    components$score[, names(stats::coef(model)), drop = FALSE],
    sandwich::estfun(model),
    tolerance = 1e-10
  )

  # The bread is the inverse observed information of the full likelihood,
  # including log(theta), not pscl's optim() Hessian.
  if (requireNamespace("numDeriv", quietly = TRUE)) {
    X <- model$x$count
    Z <- model$x$zero
    k <- ncol(X)
    m <- ncol(Z)
    total_loglik <- function(p) {
      mu <- exp(as.numeric(X %*% p[seq_len(k)]))
      pz <- model$linkinv(as.numeric(Z %*% p[k + seq_len(m)]))
      theta <- exp(p[k + m + 1L])
      sum(ifelse(
        model$y == 0,
        log(pz + (1 - pz)*stats::dnbinom(0, size = theta, mu = mu)),
        log(1 - pz) + stats::dnbinom(model$y, size = theta, mu = mu, log = TRUE)
      ))
    }
    expect_equal(
      unname(components$bread/nrow(components$score)),
      unname(solve(-numDeriv::hessian(total_loglik, parameters))),
      tolerance = 1e-6
    )
  }

  observation_loglik <- function(par, row) {
    count_names <- names(model$coefficients$count)
    zero_names <- names(model$coefficients$zero)
    count_index <- seq_along(count_names)
    zero_index <- length(count_names) + seq_along(zero_names)
    theta_index <- length(par)
    mu <- exp(sum(model$x$count[row, ]*par[count_index]))
    zero_eta <- sum(model$x$zero[row, ]*par[zero_index])
    probability_zero <- model$linkinv(zero_eta)
    theta <- exp(par[theta_index])
    y <- model$y[row]
    likelihood <- if (y == 0) {
      probability_zero +
        (1 - probability_zero)*stats::dnbinom(0, size = theta, mu = mu)
    } else {
      (1 - probability_zero)*stats::dnbinom(y, size = theta, mu = mu)
    }
    model$weights[row]*log(likelihood)
  }

  step <- 1e-6
  rows <- c(which(model$y == 0)[1L], which(model$y > 0)[1L])
  for (row in rows) {
    numerical <- vapply(seq_along(parameters), function(j) {
      upper <- lower <- parameters
      upper[j] <- upper[j] + step
      lower[j] <- lower[j] - step
      (observation_loglik(upper, row) - observation_loglik(lower, row)) /
        (2*step)
    }, numeric(1))

    expect_equal(
      unname(components$score[row, ]),
      unname(numerical),
      tolerance = 1e-6
    )
  }
})

test_that("ZIP and ZINB work together with marginaleffects", {
  skip_if_not_installed("pscl")
  dat <- zeroinfl_test_data()
  zip <- pscl::zeroinfl(
    y ~ x + z | z,
    data = dat,
    dist = "poisson",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  zinb <- pscl::zeroinfl(
    y ~ x + z | z,
    data = dat,
    dist = "negbin",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  fit <- suest(zip, zinb, model_names = c("ZIP", "ZINB"))
  effects <- marginaleffects::avg_comparisons(
    fit,
    variables = "x",
    newdata = dat
  )

  expect_true("ZINB::ln_theta" %in% names(coef(fit)))
  expect_equal(nrow(effects), 2L)
  expect_true(all(is.finite(effects$estimate)))
  expect_true(all(is.finite(effects$std.error)))
})

# insight cannot find the data of a zeroinfl model fitted inside a function
# and falls back on the model frame, which is fine for these models alone.
without_data_warning <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("could not be extracted reliably", conditionMessage(w), fixed = TRUE))
      invokeRestart("muffleWarning")
  })
}

zeroinfl_type_fits <- function() {
  dat <- zeroinfl_test_data()
  zip <- pscl::zeroinfl(y ~ x + z | z, data = dat, dist = "poisson")
  zinb <- pscl::zeroinfl(y ~ x | x + z, data = dat, dist = "negbin")
  list(
    data = dat,
    zip = zip,
    zinb = zinb,
    fit = suest(zip, zinb, model_names = c("ZIP", "ZINB"))
  )
}

test_that("zero-inflated prediction types match zeroinfl's own predictions", {
  skip_if_not_installed("pscl")
  skip_if_not_installed("marginaleffects")
  s <- zeroinfl_type_fits()
  # zeroinfl alone with the ZIP block of the joint covariance
  index <- s$fit$index[[1L]]
  V <- vcov(s$fit)[index, index]
  dimnames(V) <- list(names(coef(s$zip)), names(coef(s$zip)))

  for (type in c("response", "count", "zero", "prob")) {
    alone <- as.data.frame(without_data_warning(marginaleffects::avg_slopes(
      s$zip, variables = "x", type = type, vcov = V, newdata = s$data
    )))
    joint <- as.data.frame(marginaleffects::avg_slopes(
      s$fit, variables = "x", type = type
    ))
    zip_rows <- as.character(joint$group) == "ZIP" |
      startsWith(as.character(joint$group), "ZIP::")
    joint <- joint[zip_rows, , drop = FALSE]
    if (type == "prob") {
      joint <- joint[match(paste0("ZIP::", alone$group), joint$group), ]
    }

    expect_equal(nrow(joint), nrow(alone))
    expect_equal(joint$estimate, alone$estimate, tolerance = 1e-8)
    expect_equal(joint$std.error, alone$std.error, tolerance = 1e-6)
  }

  for (type in c("response", "count", "zero")) {
    alone <- without_data_warning(
      marginaleffects::avg_predictions(s$zinb, type = type, newdata = s$data)
    )
    joint <- as.data.frame(marginaleffects::avg_predictions(s$fit, type = type))
    expect_equal(
      joint$estimate[joint$group == "ZINB"],
      alone$estimate,
      tolerance = 1e-8
    )
  }
})

test_that("type = 'prob' gives one outcome per count, including dispersion", {
  skip_if_not_installed("pscl")
  skip_if_not_installed("marginaleffects")
  s <- zeroinfl_type_fits()
  probabilities <- as.data.frame(
    marginaleffects::avg_predictions(s$fit, type = "prob")
  )
  counts <- as.character(0:max(s$data$y))

  expect_equal(
    as.character(probabilities$group),
    c(paste0("ZIP::", counts), paste0("ZINB::", counts))
  )
  zinb_rows <- startsWith(as.character(probabilities$group), "ZINB::")
  expect_equal(
    probabilities$estimate[zinb_rows],
    unname(colMeans(predict(s$zinb, type = "prob"))),
    tolerance = 1e-8
  )

  # Pr(y = 0) for the ZINB model depends on log(theta); its standard error
  # comes from the joint covariance, including the dispersion parameter.
  index <- s$fit$index[[2L]]
  parameters <- coef(s$fit)[index]
  average_zero <- function(p) {
    m <- s$zinb
    k <- length(m$coefficients$count)
    m$coefficients$count[] <- p[seq_len(k)]
    m$coefficients$zero[] <- p[k + seq_along(m$coefficients$zero)]
    m$theta <- exp(p[length(p)])
    mean(predict(m, newdata = s$data, type = "prob", at = 0:max(m$y))[, 1L])
  }
  step <- 1e-6
  gradient <- vapply(seq_along(parameters), function(j) {
    upper <- lower <- parameters
    upper[j] <- upper[j] + step
    lower[j] <- lower[j] - step
    (average_zero(upper) - average_zero(lower)) / (2 * step)
  }, numeric(1))
  by_hand <- sqrt(drop(t(gradient) %*% vcov(s$fit)[index, index] %*% gradient))

  expect_equal(
    probabilities$std.error[probabilities$group == "ZINB::0"],
    by_hand,
    tolerance = 1e-5
  )
})

test_that("zero-inflated types work with model-specific samples", {
  skip_if_not_installed("pscl")
  skip_if_not_installed("marginaleffects")
  dat <- zeroinfl_test_data()
  first <- dat[dat$id <= 500, ]
  second <- dat[dat$id > 500, ]
  zip1 <- pscl::zeroinfl(y ~ x + z | z, data = first)
  zip2 <- pscl::zeroinfl(y ~ x + z | z, data = second)
  fit <- suest(zip1, zip2, observation_id = "id")
  newdata <- suest_newdata(fit)

  zero <- as.data.frame(
    marginaleffects::avg_predictions(fit, type = "zero", newdata = newdata)
  )
  expect_equal(
    zero$estimate,
    c(mean(predict(zip1, type = "zero")), mean(predict(zip2, type = "zero"))),
    tolerance = 1e-8
  )

  probabilities <- as.data.frame(
    marginaleffects::avg_predictions(fit, type = "prob", newdata = newdata)
  )
  expect_equal(
    probabilities$estimate[as.character(probabilities$group) == "zip2::0"],
    mean(predict(zip2, type = "prob")[, 1L]),
    tolerance = 1e-8
  )

  difference <- marginaleffects::avg_comparisons(
    fit, variables = "x", type = "count", newdata = newdata,
    hypothesis = ~ pairwise
  )
  expect_equal(nrow(difference), 1L)
  expect_true(is.finite(difference$std.error))
})

test_that("zero-inflated types need every model to be zero-inflated", {
  skip_if_not_installed("pscl")
  skip_if_not_installed("marginaleffects")
  s <- zeroinfl_type_fits()
  dat <- s$data
  zip <- pscl::zeroinfl(y ~ x + z | z, data = dat)
  logit <- glm(I(y > 0) ~ x + z, family = binomial, data = dat)
  mixed <- suest(zip, logit, model_names = c("ZIP", "Logit"))

  expect_error(
    marginaleffects::avg_predictions(mixed, type = "zero"),
    "Not zero-inflated: Logit"
  )
  expect_error(
    marginaleffects::avg_predictions(s$fit, type = "link"),
    "no single link scale"
  )
  expect_error(
    marginaleffects::avg_predictions(s$fit, type = "probs"),
    "'response', 'count', 'zero', or 'prob'"
  )
})
