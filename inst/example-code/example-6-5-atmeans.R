# Exact Stata-style atmeans calculation for Example 6.5.
#
# Stata's covariates(atmeans) operates on the model matrix: numeric columns
# and factor-indicator columns are each held at their sample means. This is
# different from marginaleffects newdata = "mean", which uses means for
# numeric variables and modes for categorical variables.

.example65_model_matrix_atmeans <- function(model, age) {
  mf <- stats::model.frame(model)
  tt <- stats::delete.response(stats::terms(model))
  X <- stats::model.matrix(tt, data = mf, contrasts.arg = model$contrasts)
  x <- colMeans(X)

  if ("age" %in% names(x))
    x["age"] <- age
  if ("I(age^2)" %in% names(x))
    x["I(age^2)"] <- age^2

  x
}

.example65_polr_probability_jacobian <- function(model, age) {
  x <- .example65_model_matrix_atmeans(model, age)
  beta <- model$coefficients
  zeta <- model$zeta
  x <- x[names(beta)]

  eta <- sum(x * beta)
  cumulative <- stats::plogis(zeta - eta)
  density <- cumulative * (1 - cumulative)
  categories <- length(zeta) + 1L

  probability <- c(cumulative[1], diff(cumulative),
                   1 - cumulative[length(cumulative)])
  derivative_eta <- c(-density[1],
                      density[-length(density)] - density[-1L],
                      density[length(density)])
  J_beta <- outer(derivative_eta, x)

  J_zeta <- matrix(0, nrow = categories, ncol = length(zeta))
  J_zeta[1L, 1L] <- density[1L]

  if (categories > 2L) {
    for (category in 2L:(categories - 1L)) {
      J_zeta[category, category - 1L] <- -density[category - 1L]
      J_zeta[category, category] <- density[category]
    }
  }

  J_zeta[categories, length(zeta)] <- -density[length(zeta)]
  J <- cbind(J_beta, J_zeta)
  colnames(J) <- c(names(beta), names(zeta))

  list(probability = probability, jacobian = J)
}

.example65_multinom_probability_jacobian <- function(model, age) {
  x <- .example65_model_matrix_atmeans(model, age)
  coefficients <- stats::coef(model)
  x <- x[colnames(coefficients)]

  eta <- c(0, as.numeric(coefficients %*% x))
  eta <- eta - max(eta)
  probability <- exp(eta) / sum(exp(eta))

  categories <- length(probability)
  regressors <- length(x)
  J <- matrix(0, nrow = categories,
              ncol = (categories - 1L) * regressors)

  for (outcome in 2L:categories) {
    derivative <- probability *
      ((seq_len(categories) == outcome) - probability[outcome])
    columns <- (outcome - 2L) * regressors + seq_len(regressors)
    J[, columns] <- outer(derivative, x)
  }

  colnames(J) <- as.vector(t(outer(
    model$lev[-1L], colnames(coefficients), paste, sep = ":"
  )))

  list(probability = probability, jacobian = J)
}

example65_atmeans <- function(object, age_lo = 20, age_hi = 30) {
  if (!inherits(object, "suest_model"))
    stop("object must be a suest_model.", call. = FALSE)
  if (!identical(unname(object$model_types), c("ologit", "multinom")))
    stop("This helper requires an ordered-logit and multinomial-logit pair.",
         call. = FALSE)

  ordered_lo <- .example65_polr_probability_jacobian(object$models[[1]], age_lo)
  ordered_hi <- .example65_polr_probability_jacobian(object$models[[1]], age_hi)
  nominal_lo <- .example65_multinom_probability_jacobian(object$models[[2]], age_lo)
  nominal_hi <- .example65_multinom_probability_jacobian(object$models[[2]], age_hi)

  ordered_effect <- ordered_hi$probability - ordered_lo$probability
  nominal_effect <- nominal_hi$probability - nominal_lo$probability
  ordered_J <- ordered_hi$jacobian - ordered_lo$jacobian
  nominal_J <- nominal_hi$jacobian - nominal_lo$jacobian

  categories <- object$models[[1]]$lev
  K <- length(categories)
  P <- length(stats::coef(object))
  J <- matrix(0, nrow = 2L * K, ncol = P)
  J[seq_len(K), object$index[[1]]] <- ordered_J
  J[K + seq_len(K), object$index[[2]]] <- nominal_J

  estimate <- c(ordered_effect, nominal_effect)
  V <- J %*% stats::vcov(object) %*% t(J)
  V <- (V + t(V)) / 2
  se <- sqrt(diag(V))
  p <- 2 * stats::pnorm(-abs(estimate / se))

  C <- cbind(diag(K), -diag(K))
  difference <- as.numeric(C %*% estimate)
  difference_vcov <- C %*% V %*% t(C)
  difference_vcov <- (difference_vcov + t(difference_vcov)) / 2
  difference_se <- sqrt(diag(difference_vcov))
  difference_p <- 2 * stats::pnorm(-abs(difference / difference_se))

  effects <- data.frame(
    category = rep(categories, 2L),
    model = rep(c("Ordered", "Multinomial"), each = K),
    estimate = estimate,
    std.error = se,
    p.value = p,
    row.names = NULL
  )
  differences <- data.frame(
    category = categories,
    estimate = difference,
    std.error = difference_se,
    p.value = difference_p,
    row.names = NULL
  )

  list(effects = effects, differences = differences)
}
