.suest_pair_key <- function(types) {
  paste(sort(types), collapse = "+")
}

.suest_pair_info <- function(types) {
  key <- .suest_pair_key(types)

  same_type <- length(unique(types)) == 1L
  allowed_mixed <- c(
    "logit+probit",
    "lm+logit",
    "lm+probit",
    "negbin+poisson",
    "multinom+ologit"
  )

  if (!same_type && !key %in% allowed_mixed)
    stop(
      paste0(
        "This model combination is not supported. Allowed cross-model ",
        "comparisons are logit-probit, logit-lm, probit-lm, ",
        "Poisson-negative binomial, and ordered logit-multinomial logit."
      ),
      call. = FALSE
    )

  scale <- if (key == "multinom+ologit" ||
               (same_type && types[1] %in% c("ologit", "oprobit", "multinom"))) {
    "category probabilities"
  } else if (key == "negbin+poisson" ||
             (same_type && types[1] %in% c("poisson", "negbin"))) {
    "expected counts"
  } else if (key %in% c("logit+probit", "lm+logit", "lm+probit") ||
             (same_type && types[1] %in% c("logit", "probit"))) {
    "predicted probabilities"
  } else {
    "fitted values"
  }

  list(
    key = key,
    mixed = !same_type,
    scale = scale
  )
}

.suest_data_source <- function(model) {
  data_call <- model$call$data
  env <- environment(stats::formula(model))
  env_key <- format(env)

  if (is.null(data_call)) {
    paste0("<formula environment>@", env_key)
  } else {
    paste0(
      paste(deparse(data_call, width.cutoff = 500L), collapse = ""),
      "@",
      env_key
    )
  }
}

.suest_extract_parameters <- function(model, type, engine) {
  engine <- .suest_engine_name(engine)
  if (identical(engine, "ordinal::clm")) {
    stats::coef(model)
  } else if (type %in% c("ologit", "oprobit")) {
    c(model$coefficients, model$zeta)
  } else if (type == "multinom") {
    cf <- stats::coef(model)
    b <- as.vector(t(cf))
    names(b) <- as.vector(t(outer(
      rownames(cf),
      colnames(cf),
      paste,
      sep = ":"
    )))
    b
  } else if (type == "negbin") {
    c(stats::coef(model), ln_theta = log(model$theta))
  } else {
    stats::coef(model)
  }
}


.suest_information_inverse <- function(A) {
  A <- (A + t(A)) / 2
  scale <- sqrt(pmax(diag(A), .Machine$double.eps))
  standardized <- A / outer(scale, scale)
  inverse <- solve(standardized) / outer(scale, scale)
  (inverse + t(inverse)) / 2
}

.suest_polr_components <- function(model) {
  mf <- stats::model.frame(model)
  X <- stats::model.matrix(model)[, -1L, drop = FALSE]
  y <- as.integer(stats::model.response(mf))
  weights <- stats::model.weights(mf)
  if (is.null(weights))
    weights <- rep.int(1, nrow(X))

  beta <- model$coefficients
  zeta <- model$zeta
  eta <- as.vector(X %*% beta)
  thresholds <- c(-Inf, zeta, Inf)

  lower <- thresholds[y] - eta
  upper <- thresholds[y + 1L] - eta

  if (identical(model$method, "logistic")) {
    F_lower <- stats::plogis(lower)
    F_upper <- stats::plogis(upper)
    f_lower <- stats::dlogis(lower)
    f_upper <- stats::dlogis(upper)
    fp_lower <- f_lower * (1 - 2 * F_lower)
    fp_upper <- f_upper * (1 - 2 * F_upper)
  } else if (identical(model$method, "probit")) {
    F_lower <- stats::pnorm(lower)
    F_upper <- stats::pnorm(upper)
    f_lower <- stats::dnorm(lower)
    f_upper <- stats::dnorm(upper)
    fp_lower <- ifelse(is.finite(lower), -lower * f_lower, 0)
    fp_upper <- ifelse(is.finite(upper), -upper * f_upper, 0)
  } else {
    stop("Unsupported ordered-response link.", call. = FALSE)
  }

  probability <- F_upper - F_lower
  if (any(!is.finite(probability)) || any(probability <= 0))
    stop(
      "Ordered-model fitted probabilities must be positive and finite.",
      call. = FALSE
    )

  d_eta <- f_lower - f_upper
  score_eta <- d_eta / probability
  U_beta <- X * (weights * score_eta)

  m <- length(zeta)
  U_zeta <- matrix(
    0,
    nrow = nrow(X),
    ncol = m,
    dimnames = list(rownames(X), names(zeta))
  )

  upper_exists <- y <= m
  lower_exists <- y > 1L

  if (any(upper_exists)) {
    U_zeta[
      cbind(which(upper_exists), y[upper_exists])
    ] <- weights[upper_exists] *
      f_upper[upper_exists] /
      probability[upper_exists]
  }

  if (any(lower_exists)) {
    U_zeta[
      cbind(which(lower_exists), y[lower_exists] - 1L)
    ] <- -weights[lower_exists] *
      f_lower[lower_exists] /
      probability[lower_exists]
  }

  U <- cbind(U_beta, U_zeta)
  colnames(U_beta) <- colnames(X)
  colnames(U) <- c(colnames(X), names(zeta))
  rownames(U) <- rownames(X)

  p_eta_eta <- fp_upper - fp_lower
  loglik_eta_eta <- p_eta_eta / probability - score_eta^2
  A_beta_beta <- crossprod(
    X,
    X * (-weights * loglik_eta_eta)
  )

  A_beta_zeta <- matrix(
    0,
    nrow = ncol(X),
    ncol = m,
    dimnames = list(colnames(X), names(zeta))
  )
  A_zeta_zeta <- matrix(
    0,
    nrow = m,
    ncol = m,
    dimnames = list(names(zeta), names(zeta))
  )

  for (threshold in seq_len(m)) {
    is_upper <- y == threshold
    is_lower <- y == threshold + 1L

    derivative <- numeric(nrow(X))

    if (any(is_upper)) {
      derivative[is_upper] <- (
        -fp_upper[is_upper] * probability[is_upper] -
          d_eta[is_upper] * f_upper[is_upper]
      ) / probability[is_upper]^2

      loglik_upper_upper <- (
        fp_upper[is_upper] / probability[is_upper] -
          (f_upper[is_upper] / probability[is_upper])^2
      )

      A_zeta_zeta[threshold, threshold] <-
        A_zeta_zeta[threshold, threshold] -
        sum(weights[is_upper] * loglik_upper_upper)
    }

    if (any(is_lower)) {
      derivative[is_lower] <- (
        fp_lower[is_lower] * probability[is_lower] +
          d_eta[is_lower] * f_lower[is_lower]
      ) / probability[is_lower]^2

      loglik_lower_lower <- (
        -fp_lower[is_lower] / probability[is_lower] -
          (f_lower[is_lower] / probability[is_lower])^2
      )

      A_zeta_zeta[threshold, threshold] <-
        A_zeta_zeta[threshold, threshold] -
        sum(weights[is_lower] * loglik_lower_lower)
    }

    A_beta_zeta[, threshold] <- -colSums(
      X * (weights * derivative)
    )
  }

  if (m >= 2L) {
    for (category in 2L:m) {
      observations <- y == category
      if (any(observations)) {
        cross_second <- (
          f_lower[observations] *
            f_upper[observations] /
            probability[observations]^2
        )
        value <- -sum(weights[observations] * cross_second)
        A_zeta_zeta[category - 1L, category] <- value
        A_zeta_zeta[category, category - 1L] <- value
      }
    }
  }

  A <- rbind(
    cbind(A_beta_beta, A_beta_zeta),
    cbind(t(A_beta_zeta), A_zeta_zeta)
  )
  dimnames(A) <- list(colnames(U), colnames(U))

  n <- nrow(X)
  B <- n * .suest_information_inverse(A)

  list(
    score = U,
    bread = B,
    parameters = .suest_extract_parameters(
      model,
      if (identical(model$method, "logistic")) "ologit" else "oprobit",
      "MASS::polr"
    )
  )
}

.suest_multinom_components <- function(model) {
  mf <- stats::model.frame(model)
  X <- stats::model.matrix(model)
  weights <- stats::model.weights(mf)
  if (is.null(weights))
    weights <- rep.int(1, nrow(X))

  probability <- stats::fitted(model)
  if (is.null(dim(probability)))
    probability <- cbind(1 - probability, probability)

  residual <- stats::residuals(model)
  if (is.null(dim(residual)))
    residual <- cbind(-residual, residual)

  categories <- model$lev
  nonreference <- seq.int(2L, length(categories))
  k <- ncol(X)
  q <- length(nonreference)
  n <- nrow(X)

  parameter_names <- as.vector(t(outer(
    categories[-1L],
    colnames(X),
    paste,
    sep = ":"
  )))

  U <- matrix(
    0,
    nrow = n,
    ncol = q * k,
    dimnames = list(rownames(X), parameter_names)
  )

  for (outcome in seq_len(q)) {
    columns <- (outcome - 1L) * k + seq_len(k)
    U[, columns] <- X *
      (weights * residual[, nonreference[outcome]])
  }

  A <- matrix(
    0,
    nrow = q * k,
    ncol = q * k,
    dimnames = list(parameter_names, parameter_names)
  )

  for (outcome1 in seq_len(q)) {
    p1 <- probability[, nonreference[outcome1]]
    columns1 <- (outcome1 - 1L) * k + seq_len(k)

    for (outcome2 in seq_len(q)) {
      p2 <- probability[, nonreference[outcome2]]
      columns2 <- (outcome2 - 1L) * k + seq_len(k)

      weight <- if (outcome1 == outcome2) {
        p1 * (1 - p1)
      } else {
        -p1 * p2
      }

      A[columns1, columns2] <- crossprod(
        X,
        X * (weights * weight)
      )
    }
  }

  B <- n * .suest_information_inverse(A)

  list(
    score = U,
    bread = B,
    parameters = .suest_extract_parameters(
      model,
      "multinom",
      "nnet::multinom"
    )
  )
}


.suest_negbin_components <- function(model) {
  mf <- stats::model.frame(model)
  X <- stats::model.matrix(model)
  y <- stats::model.response(mf)
  mu <- stats::fitted(model)
  theta <- as.numeric(model$theta)
  n <- nrow(X)

  beta_score_factor <- theta * (y - mu) / (theta + mu)
  U_beta <- X * beta_score_factor

  g <- digamma(y + theta) - digamma(theta) +
    log(theta) + 1 - log(theta + mu) -
    (theta + y) / (theta + mu)
  U_theta <- theta * g

  U <- cbind(U_beta, ln_theta = U_theta)
  colnames(U_beta) <- colnames(X)
  colnames(U) <- c(colnames(X), "ln_theta")
  rownames(U) <- rownames(X)

  weight_beta <- theta * mu * (theta + y) / (theta + mu)^2
  A_bb <- crossprod(X, X * weight_beta)

  h_beta_theta <- theta * mu * (y - mu) / (theta + mu)^2
  A_btheta <- -colSums(X * h_beta_theta)

  g_prime <- trigamma(y + theta) - trigamma(theta) +
    1 / theta - 1 / (theta + mu) +
    (y - mu) / (theta + mu)^2
  A_thetatheta <- -sum(theta * g + theta^2 * g_prime)

  A <- rbind(
    cbind(A_bb, A_btheta),
    c(A_btheta, A_thetatheta)
  )
  A <- (A + t(A)) / 2
  dimnames(A) <- list(colnames(U), colnames(U))

  B <- n * solve(A)
  B <- (B + t(B)) / 2

  list(
    score = U,
    bread = B,
    parameters = .suest_extract_parameters(
      model,
      "negbin",
      "MASS::glm.nb"
    )
  )
}

.suest_model_components <- function(model, type, engine) {
  engine <- .suest_engine_name(engine)
  if (type == "negbin") {
    .suest_negbin_components(model)
  } else if (identical(engine, "ordinal::clm")) {
    list(
      score = sandwich::estfun(model),
      bread = sandwich::bread(model),
      parameters = .suest_extract_parameters(model, type, engine)
    )
  } else if (type %in% c("ologit", "oprobit")) {
    .suest_polr_components(model)
  } else if (type == "multinom") {
    .suest_multinom_components(model)
  } else {
    list(
      score = sandwich::estfun(model),
      bread = sandwich::bread(model),
      parameters = .suest_extract_parameters(model, type, engine)
    )
  }
}


.suest_block_diag <- function(A, B) {
  out <- matrix(0, nrow(A) + nrow(B), ncol(A) + ncol(B))
  out[seq_len(nrow(A)), seq_len(ncol(A))] <- A
  out[nrow(A) + seq_len(nrow(B)), ncol(A) + seq_len(ncol(B))] <- B
  out
}

.suest_align_scores <- function(U, keys, union_keys) {
  out <- matrix(
    0,
    nrow = length(union_keys),
    ncol = ncol(U),
    dimnames = list(union_keys, colnames(U))
  )
  out[match(keys, union_keys), ] <- U
  out
}
