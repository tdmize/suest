#' Combine two fitted models with seemingly unrelated estimation
#'
#' `suest()` combines two separately fitted models and constructs a joint
#' model-robust covariance matrix from their observation-level score
#' contributions. The returned object can be passed directly to
#' [marginaleffects::predictions()], [marginaleffects::avg_comparisons()],
#' [marginaleffects::avg_slopes()], and [marginaleffects::hypotheses()].
#'
#' @param model1,model2 Two supported fitted model objects.
#' @param model_names Optional character vector containing two display names.
#'   By default, the object names supplied in the call are used.
#'
#' @return An object of class `"suest_model"` containing the two fitted models,
#'   their joint coefficient vector, a joint model-robust covariance matrix,
#'   and sample-alignment information.
#'
#' @details
#' Exactly two models are supported. Models can use identical, partially
#' overlapping, or completely disjoint samples. When the model calls refer to
#' the same data source, model-frame row names identify overlapping
#' observations. Models fitted from different data objects are treated as
#' disjoint because shared observations cannot be inferred safely without an
#' identifier.
#'
#' Same-family comparisons are supported for all model types listed below.
#' The supported cross-family pairs are logit--probit, logit--linear,
#' probit--linear, Poisson--negative binomial, and ordered
#' logit--multinomial logit.
#'
#' Negative-binomial models include `log(theta)` in the joint parameter vector.
#' Ordered and multinomial models use analytic score and observed-information
#' calculations for stable robust covariance estimation.
#'
#' Offsets, nonunit estimation weights, and aliased parameters are not
#' currently supported.
#'
#' @section Supported models:
#' * [stats::lm()]
#' * binary logit and probit models from [stats::glm()]
#' * Poisson log-link models from [stats::glm()]
#' * negative-binomial log-link models from [MASS::glm.nb()]
#' * ordered logit and probit models from [MASS::polr()]
#' * multinomial logit models from [nnet::multinom()]
#'
#' @examples
#' dat <- mtcars
#' dat$am <- factor(dat$am)
#'
#' model1 <- glm(am ~ wt, family = binomial(), data = dat)
#' model2 <- glm(am ~ wt + hp, family = binomial(), data = dat)
#'
#' fit <- suest(model1, model2, model_names = c("Base", "Adjusted"))
#' fit
#'
#' effects <- marginaleffects::avg_comparisons(
#'   fit,
#'   variables = "wt",
#'   newdata = dat
#' )
#' marginaleffects::hypotheses(
#'   effects,
#'   hypothesis = difference ~ revpairwise
#' )
#'
#' @export
suest <- function(model1, model2, model_names = NULL) {
  if (!requireNamespace("sandwich", quietly = TRUE))
    stop("Package 'sandwich' is required.", call. = FALSE)
  if (!requireNamespace("marginaleffects", quietly = TRUE))
    stop("Package 'marginaleffects' is required.", call. = FALSE)

  models <- list(model1, model2)

  if (is.null(model_names)) {
    model_names <- c(deparse(substitute(model1)), deparse(substitute(model2)))
    if (any(nchar(model_names) > 30) || anyDuplicated(model_names))
      model_names <- c("Model 1", "Model 2")
  }
  if (length(model_names) != 2L || anyNA(model_names) ||
      any(model_names == ""))
    stop("'model_names' must contain two nonempty names.", call. = FALSE)

  model_names <- make.unique(as.character(model_names))
  names(models) <- model_names

  model_types <- vapply(models, .suest_model_type, character(1))

  if (any(model_types == "unsupported"))
    stop(
      paste0(
        "Each model must be a supported lm, binary logit, binary probit, ",
        "Poisson, MASS::glm.nb, MASS::polr, or nnet::multinom model."
      ),
      call. = FALSE
    )

  pair <- .suest_pair_info(model_types)

  if (any(model_types == "multinom") &&
      !requireNamespace("nnet", quietly = TRUE))
    stop(
      "Package 'nnet' is required for multinomial models.",
      call. = FALSE
    )

  categorical <- model_types %in% c("ologit", "oprobit", "multinom")

  if (all(categorical)) {
    levels1 <- models[[1]]$lev
    levels2 <- models[[2]]$lev
    if (!identical(levels1, levels2))
      stop(
        paste0(
          "Categorical-outcome models must use the same outcome categories ",
          "in the same order."
        ),
        call. = FALSE
      )
  }

  model_frames <- lapply(models, stats::model.frame)

  for (i in seq_along(models)) {
    model <- models[[i]]
    type <- model_types[i]
    mf <- model_frames[[i]]

    parameters <- .suest_extract_parameters(model, type)
    if (anyNA(parameters))
      stop(
        sprintf(
          "Model '%s' contains aliased or missing parameters.",
          model_names[i]
        ),
        call. = FALSE
      )

    offset <- stats::model.offset(mf)
    weights <- stats::model.weights(mf)

    if (!is.null(offset) && any(offset != 0))
      stop("Offsets are not supported in this draft.", call. = FALSE)
    if (!is.null(weights) && any(weights != 1))
      stop("Weights are not supported in this draft.", call. = FALSE)

    y <- stats::model.response(mf)

    if (type == "lm" && (!is.numeric(y) || is.matrix(y)))
      stop("Linear-model outcomes must be numeric vectors.",
           call. = FALSE)

    if (type %in% c("logit", "probit")) {
      binary <- if (is.factor(y)) {
        nlevels(y) == 2L
      } else if (is.matrix(y)) {
        FALSE
      } else {
        all(stats::na.omit(unique(y)) %in% c(0, 1))
      }

      if (!binary)
        stop(
          "Binary logit and probit outcomes must be two-level factors or numeric 0/1.",
          call. = FALSE
        )
    }

    if (type %in% c("poisson", "negbin")) {
      count <- is.numeric(y) && !is.matrix(y) && all(is.finite(y)) &&
        all(y >= 0) &&
        all(abs(y - round(y)) < sqrt(.Machine$double.eps))

      if (!count)
        stop(
          "Poisson and negative-binomial outcomes must be nonnegative integer counts.",
          call. = FALSE
        )
    }

    if (type %in% c("ologit", "oprobit") &&
        (!is.ordered(y) || nlevels(y) < 3L))
      stop(
        "Ordered-model outcomes must be ordered factors with at least three levels.",
        call. = FALSE
      )

    if (type == "multinom") {
      cf <- stats::coef(model)
      valid_multinom <- !is.null(model$lev) &&
        length(model$lev) >= 3L &&
        is.matrix(cf) &&
        nrow(cf) == length(model$lev) - 1L

      if (!valid_multinom)
        stop(
          "Multinomial models must have at least three outcome categories.",
          call. = FALSE
        )
    }
  }

  components <- Map(.suest_model_components, models, model_types)
  scores <- lapply(components, `[[`, "score")
  breads <- lapply(components, `[[`, "bread")
  parameters <- lapply(components, `[[`, "parameters")
  local_names <- lapply(parameters, names)

  for (i in seq_along(models)) {
    U <- scores[[i]]
    B <- breads[[i]]
    mf <- model_frames[[i]]

    if (nrow(U) != nrow(mf))
      stop(
        paste0(
          "Each model must provide one score contribution per observation. ",
          "For nnet::multinom(), use the default summ = 0."
        ),
        call. = FALSE
      )

    if (nrow(B) != ncol(U) || ncol(B) != ncol(U))
      stop(
        sprintf(
          "The bread and score dimensions do not match for model '%s'.",
          model_names[i]
        ),
        call. = FALSE
      )

    if (!identical(colnames(U), local_names[[i]]))
      stop(
        sprintf(
          "The score parameter order does not match model '%s'.",
          model_names[i]
        ),
        call. = FALSE
      )

    if (nrow(U) <= 1L)
      stop("Each model must contain at least two observations.",
           call. = FALSE)
  }

  sources <- vapply(models, .suest_data_source, character(1))
  sample_rows <- lapply(model_frames, rownames)

  for (i in seq_along(sample_rows)) {
    if (is.null(sample_rows[[i]]))
      sample_rows[[i]] <- as.character(seq_len(nrow(model_frames[[i]])))
    if (anyDuplicated(sample_rows[[i]]))
      stop(
        sprintf(
          "Model '%s' has duplicated model-frame row names and cannot be aligned automatically.",
          model_names[i]
        ),
        call. = FALSE
      )
  }

  sample_keys <- Map(
    function(source, rows) paste(source, rows, sep = "\r"),
    sources,
    sample_rows
  )

  union_keys <- unique(c(sample_keys[[1]], sample_keys[[2]]))
  overlap_keys <- intersect(sample_keys[[1]], sample_keys[[2]])

  n_model <- vapply(scores, nrow, integer(1))
  n_union <- length(union_keys)
  n_overlap <- length(overlap_keys)

  corrected_scores <- Map(
    function(U, n) U * sqrt(n / (n - 1)),
    scores,
    n_model
  )

  aligned_scores <- Map(
    .suest_align_scores,
    corrected_scores,
    sample_keys,
    MoreArgs = list(union_keys = union_keys)
  )

  joint_breads <- Map(
    function(B, n) B * n_union / n,
    breads,
    n_model
  )

  U <- cbind(aligned_scores[[1]], aligned_scores[[2]])
  B <- .suest_block_diag(joint_breads[[1]], joint_breads[[2]])
  meat <- crossprod(U) / n_union
  V <- (B %*% meat %*% B) / n_union
  V <- (V + t(V)) / 2

  joint_names <- unlist(Map(
    function(model, term) paste0(model, "::", term),
    rep(model_names, lengths(local_names)),
    unlist(local_names, use.names = FALSE)
  ), use.names = FALSE)

  b <- unlist(parameters, use.names = FALSE)
  names(b) <- joint_names
  dimnames(V) <- list(joint_names, joint_names)

  p <- lengths(local_names)
  index <- list(seq_len(p[1]), p[1] + seq_len(p[2]))
  names(index) <- model_names

  theta <- setNames(rep(NA_real_, 2L), model_names)
  for (i in seq_along(models)) {
    if (model_types[i] == "negbin")
      theta[i] <- as.numeric(models[[i]]$theta)
  }
  if (all(is.na(theta)))
    theta <- NULL

  out <- list(
    models = models,
    model_frames = model_frames,
    coefficients = b,
    vcov = V,
    model_names = model_names,
    local_names = local_names,
    index = index,
    model_types = model_types,
    pair_key = pair$key,
    mixed_models = pair$mixed,
    comparison_scale = pair$scale,
    theta = theta,
    sample_sources = sources,
    sample_rows = sample_rows,
    sample_keys = sample_keys,
    nobs_models = setNames(n_model, model_names),
    nobs_union = n_union,
    nobs_overlap = n_overlap,
    call = match.call()
  )
  class(out) <- "suest_model"
  out
}

#' @rdname suest
#' @export
coef.suest_model <- function(object, ...) object$coefficients

#' @rdname suest
#' @export
vcov.suest_model <- function(object, ...) object$vcov

#' @rdname suest
#' @export
nobs.suest_model <- function(object, ...) object$nobs_union

#' @rdname suest
#' @export
print.suest_model <- function(x, ...) {
  cat("Seemingly Unrelated Estimation\n")
  cat("Models:", paste(x$model_names, collapse = " + "), "\n")
  cat(
    "Model types:",
    paste(
      paste0(x$model_names, "=", unname(x$model_types)),
      collapse = ", "
    ),
    "\n"
  )
  cat("Comparison scale:", x$comparison_scale, "\n")

  if (!is.null(x$theta)) {
    keep <- !is.na(x$theta)
    cat(
      "Theta:",
      paste(
        paste0(
          names(x$theta)[keep],
          "=",
          format(x$theta[keep], digits = 5)
        ),
        collapse = ", "
      ),
      "\n"
    )
  }

  cat(
    "Observations:",
    paste(
      paste0(names(x$nobs_models), "=", format(x$nobs_models, big.mark = ",")),
      collapse = ", "
    ),
    "\n"
  )
  cat("Overlapping observations:", format(x$nobs_overlap, big.mark = ","), "\n")
  cat("Union observations:", format(x$nobs_union, big.mark = ","), "\n")
  cat("Parameters:", length(x$coefficients), "\n")
  invisible(x)
}

# Create model-specific newdata using each model's own estimation sample.
# This is useful when models were estimated on different samples and effects
# should be averaged separately within those samples.
