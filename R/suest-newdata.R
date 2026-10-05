#' Stack the estimation samples from a SUEST object
#'
#' Creates a data frame that stacks each model's own estimation sample. Use it
#' when the models were fitted on different samples, so that each model's
#' marginal effects are averaged over its own sample.
#'
#' @param object A `"suest_model"` returned by [suest()].
#'
#' Variables used inside formula transformations, such as `g` in `factor(g)`,
#' are included from the data the model was fit to.
#'
#' @return A data frame with the component model frames stacked vertically and
#'   internal columns `.suest_model` and `.suest_rowid`, used to route rows
#'   to the correct component model. For pweighted and survey models, `.suest_weight`
#'   contains each row's evaluated sampling weight and can be supplied to the
#'   `wts` argument of `marginaleffects` averaging functions.
#'
#' @examples
#' dat <- mtcars
#'
#' model1 <- lm(mpg ~ wt + hp, data = dat, subset = cyl == 4)
#' model2 <- lm(mpg ~ wt + hp, data = dat, subset = cyl != 4)
#' fit <- suest(model1, model2, model_names = c("Four cylinders", "Other"))
#' nd <- suest_newdata(fit)
#' marginaleffects::avg_comparisons(fit, variables = "wt", newdata = nd)
#'
#' @export
suest_newdata <- function(object) {
  if (!inherits(object, "suest_model"))
    stop("'object' must be a suest_model.", call. = FALSE)

  raw_data <- object$raw_data
  if (is.null(raw_data)) raw_data <- vector("list", length(object$model_frames))
  model_frames <- Map(function(x, raw) {
    for (nm in setdiff(names(raw), names(x))) x[[nm]] <- raw[[nm]]
    # Matrix columns such as Surv() responses or poly() terms cannot be
    # stacked; predictions rebuild them from the raw variables.
    keep <- vapply(x, function(v) is.atomic(v) && is.null(dim(v)), logical(1))
    x[keep]
  }, object$model_frames, raw_data)

  reserved <- c(".suest_model", ".suest_rowid")
  if (any(object$weight_type %in% c("pweight", "survey")))
    reserved <- c(reserved, ".suest_weight")
  if (any(vapply(
    model_frames,
    function(x) any(reserved %in% names(x)),
    logical(1)
  )))
    stop(
      "The model data contain a reserved '.suest_' column name.",
      call. = FALSE
    )

  # Stack through character representations so factors with different level
  # sets cannot be corrupted by row binding, then restore factor classes from
  # the component model frames below.
  factor_names <- unique(unlist(lapply(
    model_frames,
    function(x) names(x)[vapply(x, is.factor, logical(1))]
  ), use.names = FALSE))

  factor_specs <- lapply(factor_names, function(nm) {
    present <- lapply(model_frames, function(x) {
      if (nm %in% names(x)) x[[nm]] else NULL
    })
    present <- Filter(Negate(is.null), present)
    factors <- Filter(is.factor, present)

    levs <- unique(unlist(lapply(factors, levels), use.names = FALSE))
    observed <- unique(unlist(lapply(present, as.character), use.names = FALSE))
    levs <- unique(c(levs, observed[!is.na(observed)]))

    same_ordered_levels <- length(factors) > 0L &&
      all(vapply(factors, is.ordered, logical(1))) &&
      all(vapply(
        factors,
        function(x) identical(levels(x), levels(factors[[1L]])),
        logical(1)
      ))

    list(
      levels = levs,
      ordered = length(factors) == length(present) && same_ordered_levels
    )
  })
  names(factor_specs) <- factor_names

  row_offset <- c(0L, cumsum(vapply(
    model_frames,
    nrow,
    integer(1)
  ))[-length(model_frames)])

  frames <- lapply(seq_along(model_frames), function(i) {
    x <- model_frames[[i]]

    for (nm in names(x)) {
      if (is.factor(x[[nm]]))
        x[[nm]] <- as.character(x[[nm]])
    }

    x$.suest_model <- object$model_names[i]
    x$.suest_rowid <- row_offset[i] + seq_len(nrow(x))
    if (any(object$weight_type %in% c("pweight", "survey")))
      x$.suest_weight <- object$model_weights[[i]]
    x
  })

  all_names <- unique(unlist(lapply(frames, names), use.names = FALSE))

  frames <- lapply(frames, function(x) {
    missing <- setdiff(all_names, names(x))
    for (nm in missing)
      x[[nm]] <- NA
    x[all_names]
  })

  out <- do.call(rbind, frames)
  rownames(out) <- NULL

  for (nm in names(factor_specs)) {
    spec <- factor_specs[[nm]]
    values <- as.character(out[[nm]])
    out[[nm]] <- if (isTRUE(spec$ordered)) {
      ordered(values, levels = spec$levels)
    } else {
      factor(values, levels = spec$levels)
    }
  }

  class(out) <- c("suest_newdata", class(out))
  out
}

# marginaleffects extension methods

# Data marginaleffects uses when `newdata` is omitted, and that `datagrid()`
# builds its grid from. When every model has the same estimation sample, one
# row per observation with every model's variables: each model is then
# averaged over its own sample, and grids include every model's predictors.
# Otherwise each model's own sample, stacked as in suest_newdata().
.suest_default_data <- function(object) {
  stacked <- suest_newdata(object)
  keys <- object$sample_keys
  same_sample <- length(keys) > 1L &&
    all(vapply(keys[-1L], identical, logical(1), keys[[1L]]))
  if (!same_sample) return(stacked)

  blocks <- split(stacked, factor(stacked$.suest_model,
    levels = object$model_names))
  out <- blocks[[1L]]
  for (block in blocks[-1L]) {
    for (nm in names(out)) {
      fill <- is.na(out[[nm]]) & !is.na(block[[nm]])
      if (any(fill)) out[[nm]][fill] <- block[[nm]][fill]
    }
  }
  out$.suest_model <- NULL
  out$.suest_rowid <- NULL
  rownames(out) <- NULL
  class(out) <- "data.frame"
  out
}
