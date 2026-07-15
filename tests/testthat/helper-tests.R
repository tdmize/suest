offdiag_vcov <- function(object) {
  object$vcov[
    object$index[[1]],
    object$index[[2]],
    drop = FALSE
  ]
}

robust_vcov_direct <- function(model) {
  U <- sandwich::estfun(model)
  B <- sandwich::bread(model)
  n <- nrow(U)
  B %*% crossprod(U) %*% B / n^2 * n / (n - 1)
}
