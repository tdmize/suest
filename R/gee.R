# GEE mean-equation scores with the fitted working correlation held fixed.
# Stata's PA route retains regression parameters, not scale/correlation.
.suest_gee_components <- function(model) {
  X <- stats::model.matrix(model)
  y <- stats::model.response(model$model)
  mu <- model$fitted.values
  derivative <- model$family$mu.eta(model$linear.predictors)
  variance <- model$family$variance(mu)
  phi <- unname(model$geese$gamma)
  if (!is.finite(phi) || phi <= 0 || any(!is.finite(variance) | variance <= 0))
    stop("GEE working variances must be finite and positive.", call. = FALSE)
  rho <- if (model$corstr == "independence") 0 else unname(model$geese$alpha)
  U <- matrix(0, nrow(X), ncol(X), dimnames = dimnames(X))
  information <- matrix(0, ncol(X), ncol(X))
  for (rows in split(seq_len(nrow(X)), model$id)) {
    R <- matrix(rho, length(rows), length(rows))
    diag(R) <- 1
    covariance <- phi * R * outer(sqrt(variance[rows]), sqrt(variance[rows]))
    precision <- chol2inv(chol(covariance))
    D <- derivative[rows] * X[rows, , drop = FALSE]
    U[rows, ] <- D * drop(precision %*% (y[rows] - mu[rows]))
    information <- information + crossprod(D, precision %*% D)
  }
  B <- nrow(X) * solve(information)
  dimnames(B) <- list(colnames(X), colnames(X))
  if (any(!is.finite(U)) || any(!is.finite(B)))
    stop("Unable to construct finite GEE scores and bread.", call. = FALSE)
  list(score = U, bread = B, parameters = stats::coef(model))
}
