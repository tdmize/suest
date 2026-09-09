#' Pool SUEST systems across multiple imputations
#'
#' `suest_mi()` pools complete joint SUEST coefficient vectors and covariance
#' matrices using Rubin's rules. Supply either a list whose elements are
#' compatible [suest()] results fitted in separate imputed datasets, or the
#' `mice::mira` object returned by `with()` when that expression fits a SUEST
#' system. Result lists created by `mitools::with.imputationList()` are also
#' accepted directly.
#'
#' @param fits A list containing at least two compatible `"suest_model"`
#'   objects, one per imputation; a `"mira"` object containing those fits; or
#'   a result list returned by `mitools::with.imputationList()`.
#' @param object,x A `"suest_mi"` object.
#' @param conf.level Confidence level for [summary.suest_mi()].
#' @param ... Additional arguments; currently ignored.
#'
#' @details
#' The pooled covariance is
#' `Ubar + (1 + 1 / m) B`, where `Ubar` is the mean within-imputation joint
#' covariance, `B` is the between-imputation covariance of the complete joint
#' coefficient vector, and `m` is the number of imputations. Because each
#' within-imputation input is a complete SUEST system, cross-model covariance
#' is retained in both components.
#'
#' A `"mira"` object from `mice::with()` is read through its public
#' `analyses` component. The list returned by `mice::getfit()` is accepted as
#' well. The `mice` package is not required when pooling an ordinary list.
#'
#' Results from `mitools::with.imputationList()` use the same list-based route,
#' with their originating call retained in the returned object. Legacy
#' `"imputationResultList"` wrappers are accepted as well. The `mitools`
#' package is not required for an already-created result list.
#'
#' This initial interface supports coefficients, covariance matrices, summaries,
#' and coefficient-level hypotheses. Predictions, slopes, and comparisons must
#' be estimated within each imputation and pooled as estimands; they are not yet
#' implemented for `"suest_mi"` objects.
#'
#' @return A `"suest_mi"` object containing pooled coefficients, total,
#'   within-imputation, and between-imputation covariance matrices, Rubin
#'   large-sample degrees of freedom, and the original SUEST fits.
#'
#' @export
suest_mi <- function(fits) {
  input_source <- "list"
  mice_call <- mice_call1 <- NULL
  mitools_call <- NULL
  if (inherits(fits, "mira")) {
    input_source <- "mice::mira"
    is_wrapper <- all(c("call", "call1", "nmis", "analyses") %in% names(fits))
    if (is_wrapper) {
      mice_call <- fits[["call"]]
      mice_call1 <- fits[["call1"]]
      fits <- fits[["analyses"]]
    } else {
      # mice::getfit() returns the analyses themselves with classes
      # c("mira", "list"), rather than the usual wrapper.
      fits <- unclass(fits)
    }
  } else if (inherits(fits, "imputationResultList")) {
    input_source <- "mitools::imputationList"
    mitools_call <- fits[["call"]]
    fits <- unclass(fits)
    fits[["call"]] <- NULL
  } else if (is.list(fits) && !is.data.frame(fits) &&
      !is.null(attr(fits, "call")) &&
      all(vapply(fits, inherits, logical(1), what = "suest_model"))) {
    input_source <- "mitools::imputationList"
    mitools_call <- attr(fits, "call")
  }
  if (!is.list(fits) || is.data.frame(fits) || length(fits) < 2L)
    stop("'fits' must contain at least two SUEST objects.",
         call. = FALSE)
  if (!all(vapply(fits, inherits, logical(1), what = "suest_model")))
    stop("Every element of 'fits' must be a suest_model object.",
         call. = FALSE)

  reference <- fits[[1L]]
  fields <- c(
    "model_names", "local_names", "model_types", "model_engines",
    "category_levels", "pair_key", "mixed_models", "comparison_scale",
    "weight_type", "nobs_models", "nobs_union", "n_clusters"
  )
  for (i in seq_along(fits)[-1L]) {
    different <- fields[!vapply(fields, function(field)
      identical(fits[[i]][[field]], reference[[field]]), logical(1))]
    if (length(different))
      stop(sprintf(
        "Imputation %d is incompatible with imputation 1: %s differ%s.",
        i, paste(different, collapse = ", "), if (length(different) == 1L) "s" else ""
      ), call. = FALSE)
  }

  coefficient_names <- names(stats::coef(reference))
  if (is.null(coefficient_names) || anyNA(coefficient_names) ||
      any(coefficient_names == "") || anyDuplicated(coefficient_names))
    stop("Every pooled coefficient must have a unique nonempty name.",
         call. = FALSE)

  coefficients <- vapply(seq_along(fits), function(i) {
    value <- stats::coef(fits[[i]])
    if (!identical(names(value), coefficient_names))
      stop(sprintf(
        "Imputation %d has different coefficient names or ordering.", i
      ), call. = FALSE)
    if (any(!is.finite(value)))
      stop(sprintf("Imputation %d has nonfinite coefficients.", i),
           call. = FALSE)
    unname(value)
  }, numeric(length(coefficient_names)))
  rownames(coefficients) <- coefficient_names

  covariances <- lapply(seq_along(fits), function(i) {
    value <- stats::vcov(fits[[i]])
    if (!is.matrix(value) ||
        !identical(rownames(value), coefficient_names) ||
        !identical(colnames(value), coefficient_names))
      stop(sprintf(
        "Imputation %d has a covariance matrix with incompatible names or dimensions.", i
      ), call. = FALSE)
    if (any(!is.finite(value)))
      stop(sprintf("Imputation %d has a nonfinite covariance matrix.", i),
           call. = FALSE)
    (value + t(value))/2
  })

  m <- length(fits)
  pooled <- rowMeans(coefficients)
  within <- Reduce(`+`, covariances)/m
  centered <- sweep(coefficients, 1L, pooled)
  between <- tcrossprod(centered)/(m - 1L)
  total <- within + (1 + 1/m)*between
  total <- (total + t(total))/2
  dimnames(within) <- dimnames(between) <- dimnames(total) <-
    list(coefficient_names, coefficient_names)
  names(pooled) <- coefficient_names
  if (any(diag(total) <= 0) || any(!is.finite(total)))
    stop("Rubin pooling produced a nonpositive or nonfinite variance.",
         call. = FALSE)

  added <- (1 + 1/m)*diag(between)
  within_diagonal <- diag(within)
  relative_increase <- ifelse(within_diagonal > 0,
    added/within_diagonal, Inf)
  degrees_freedom <- ifelse(added == 0, Inf,
    (m - 1L)*(1 + 1/relative_increase)^2)
  fraction_missing <- added/diag(total)
  names(relative_increase) <- names(degrees_freedom) <-
    names(fraction_missing) <- coefficient_names

  out <- list(
    fits = fits,
    coefficients = pooled,
    vcov = total,
    within_vcov = within,
    between_vcov = between,
    degrees_freedom = degrees_freedom,
    relative_increase_variance = relative_increase,
    fraction_missing_information = fraction_missing,
    m = m,
    model_names = reference$model_names,
    model_types = reference$model_types,
    comparison_scale = reference$comparison_scale,
    nobs_models = reference$nobs_models,
    nobs_union = reference$nobs_union,
    input_source = input_source,
    mice_call = mice_call,
    mice_call1 = mice_call1,
    mitools_call = mitools_call,
    call = match.call()
  )
  class(out) <- "suest_mi"
  out
}

#' @rdname suest_mi
#' @export
coef.suest_mi <- function(object, ...) object$coefficients

#' @rdname suest_mi
#' @export
vcov.suest_mi <- function(object, ...) object$vcov

#' @rdname suest_mi
#' @export
nobs.suest_mi <- function(object, ...) object$nobs_union

#' @rdname suest_mi
#' @export
summary.suest_mi <- function(object, conf.level = 0.95, ...) {
  if (length(conf.level) != 1L || !is.finite(conf.level) ||
      conf.level <= 0 || conf.level >= 1)
    stop("'conf.level' must be one number strictly between zero and one.",
         call. = FALSE)
  estimate <- stats::coef(object)
  std.error <- sqrt(diag(stats::vcov(object)))
  df <- object$degrees_freedom
  critical <- stats::qt((1 + conf.level)/2, df = df)
  statistic <- estimate/std.error
  data.frame(
    term = names(estimate), estimate = unname(estimate),
    std.error = unname(std.error), statistic = unname(statistic),
    df = unname(df), p.value = 2*stats::pt(-abs(statistic), df = df),
    conf.low = unname(estimate - critical*std.error),
    conf.high = unname(estimate + critical*std.error),
    row.names = NULL, check.names = FALSE
  )
}

#' @rdname suest_mi
#' @export
print.suest_mi <- function(x, ...) {
  cat("Multiply Imputed Seemingly Unrelated Estimation\n")
  cat("Imputations:", x$m, "\n")
  if (!identical(x$input_source, "list"))
    cat("Input:", x$input_source, "\n")
  cat("Models:", paste(x$model_names, collapse = " + "), "\n")
  cat("Comparison scale:", x$comparison_scale, "\n")
  cat("Union observations per imputation:",
      format(x$nobs_union, big.mark = ","), "\n")
  cat("Parameters:", length(x$coefficients), "\n")
  invisible(x)
}
