offdiag_vcov <- function(object) {
  object$vcov[
    object$index[[1]],
    object$index[[2]],
    drop = FALSE
  ]
}

robust_vcov_direct <- function(model, correction = TRUE) {
  U <- sandwich::estfun(model)
  B <- sandwich::bread(model)
  n <- nrow(U)
  factor <- if (correction) n / (n - 1) else 1
  B %*% crossprod(U) %*% B / n^2 * factor
}
