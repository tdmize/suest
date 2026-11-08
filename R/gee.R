# GEE mean-equation scores with the fitted working correlation held fixed.
# Stata's PA route retains regression parameters, not scale/correlation.

# Wave codes as geeglm uses them: the integer codes of factor(waves) over the
# retained sample, or the position within each panel when waves is omitted.
.suest_gee_waves <- function(model) {
  id <- model$id
  if (is.null(model$call$waves))
    return(as.integer(stats::ave(seq_along(id), id, FUN = seq_along)))
  data <- model$data
  if (!is.data.frame(data))
    stop("GEE models with waves must be fit with a data frame in data =.", call. = FALSE)
  values <- tryCatch(eval(model$call$waves, data, environment(stats::formula(model))),
    error = function(e) NULL)
  rows <- match(rownames(model$model), rownames(data))
  if (is.null(values) || length(values) != nrow(data) || anyNA(rows) || anyNA(values[rows]))
    stop("Unable to recover the GEE waves for the retained model sample.", call. = FALSE)
  as.integer(factor(values[rows]))
}

# Per-panel working correlation matrices for ar1 and unstructured fits.
.suest_gee_correlations <- function(model) {
  waves <- .suest_gee_waves(model)
  panels <- split(seq_along(model$id), factor(model$id, levels = unique(model$id)))
  alpha <- model$geese$alpha
  if (model$corstr == "ar1") {
    if (length(alpha) != 1L || !is.finite(alpha) || abs(alpha) >= 1)
      stop("The GEE ar1 correlation must be finite and inside (-1, 1).", call. = FALSE)
    return(lapply(panels, function(rows)
      unname(alpha)^abs(outer(waves[rows], waves[rows], "-"))))
  }
  if (!all(vapply(panels, function(rows) identical(waves[rows], seq_along(rows)),
      logical(1))))
    stop(paste0("GEE unstructured correlation requires each panel's waves to run ",
      "1, 2, ... without gaps or reordering."), call. = FALSE)
  size <- max(lengths(panels))
  pairs <- outer(seq_len(size), seq_len(size), function(j, k) paste0("alpha.", j, ":", k))
  if (length(alpha) != size*(size - 1L)/2L || !all(pairs[upper.tri(pairs)] %in% names(alpha)) ||
      any(!is.finite(alpha)))
    stop("Unable to align the GEE unstructured correlations with the waves.", call. = FALSE)
  R <- diag(size)
  R[upper.tri(R)] <- alpha[pairs[upper.tri(pairs)]]
  R[lower.tri(R)] <- t(R)[lower.tri(R)]
  lapply(panels, function(rows) R[seq_along(rows), seq_along(rows), drop = FALSE])
}

.suest_gee_components <- function(model) {
  X <- stats::model.matrix(model)
  y <- stats::model.response(model$model)
  mu <- model$fitted.values
  derivative <- model$family$mu.eta(model$linear.predictors)
  variance <- model$family$variance(mu)
  phi <- unname(model$geese$gamma)
  if (!is.finite(phi) || phi <= 0 || any(!is.finite(variance) | variance <= 0))
    stop("GEE working variances must be finite and positive.", call. = FALSE)
  panels <- split(seq_len(nrow(X)), factor(model$id, levels = unique(model$id)))
  correlations <- if (model$corstr %in% c("ar1", "unstructured")) {
    .suest_gee_correlations(model)
  } else {
    rho <- if (model$corstr == "independence") 0 else unname(model$geese$alpha)
    lapply(panels, function(rows) {
      R <- matrix(rho, length(rows), length(rows))
      diag(R) <- 1
      R
    })
  }
  U <- matrix(0, nrow(X), ncol(X), dimnames = dimnames(X))
  information <- matrix(0, ncol(X), ncol(X))
  for (i in seq_along(panels)) {
    rows <- panels[[i]]
    covariance <- phi * correlations[[i]] * outer(sqrt(variance[rows]), sqrt(variance[rows]))
    factor <- tryCatch(chol(covariance), error = function(e) NULL)
    if (is.null(factor))
      stop("A GEE working covariance is not positive definite.", call. = FALSE)
    precision <- chol2inv(factor)
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
