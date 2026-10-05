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

# Observed-information sandwich from a per-observation log likelihood written
# in one or two linear predictors; derivatives by central differences.
observed_sandwich_direct <- function(ll, X, eta, Z = NULL, zeta = NULL) {
  step <- function(x) 1e-3 * (1 + abs(x))
  d1 <- function(f, x) { h <- step(x); (f(x + h) - f(x - h)) / (2 * h) }
  d2 <- function(f, x) { h <- step(x); (f(x + h) - 2 * f(x) + f(x - h)) / h^2 }
  if (is.null(Z)) {
    f <- function(e) ll(e)
    U <- X * d1(f, eta)
    H <- crossprod(X, X * d2(f, eta))
  } else {
    fe <- function(e) ll(e, zeta)
    fz <- function(z) ll(eta, z)
    he <- step(eta); hz <- step(zeta)
    mixed <- (ll(eta + he, zeta + hz) - ll(eta + he, zeta - hz) -
                ll(eta - he, zeta + hz) + ll(eta - he, zeta - hz)) / (4 * he * hz)
    U <- cbind(X * d1(fe, eta), Z * d1(fz, zeta))
    H <- rbind(cbind(crossprod(X, X * d2(fe, eta)), crossprod(X, Z * mixed)),
               cbind(crossprod(Z, X * mixed), crossprod(Z, Z * d2(fz, zeta))))
  }
  n <- nrow(U)
  Hinv <- solve(-H)
  Hinv %*% crossprod(U) %*% Hinv * n / (n - 1)
}

glm_observed_vcov_direct <- function(model) {
  y <- model$y
  inv <- model$family$linkinv
  ll <- switch(tolower(model$family$family),
    binomial = , quasibinomial = function(e) y * log(inv(e)) + (1 - y) * log(1 - inv(e)),
    gamma = function(e) -y / inv(e) - log(inv(e)),
    gaussian = function(e) -(y - inv(e))^2 / 2)
  observed_sandwich_direct(ll, stats::model.matrix(model), model$linear.predictors)
}

betareg_observed_vcov_direct <- function(model) {
  y <- model$y
  mu <- model$link$mean$linkinv
  phi <- model$link$precision$linkinv
  X <- stats::model.matrix(model, model = "mean")
  Z <- stats::model.matrix(model, model = "precision")
  k <- ncol(X)
  ll <- function(e, z) {
    m <- mu(e); p <- phi(z)
    lgamma(p) - lgamma(m * p) - lgamma((1 - m) * p) +
      (m * p - 1) * log(y) + ((1 - m) * p - 1) * log1p(-y)
  }
  observed_sandwich_direct(ll, X, drop(X %*% coef(model)[seq_len(k)]),
                           Z, drop(Z %*% coef(model)[-seq_len(k)]))
}
