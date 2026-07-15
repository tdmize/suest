#' Stack the estimation samples from a SUEST object
#'
#' Creates a data frame containing each component model's own estimation
#' sample. This is useful when the two models were fitted on different samples
#' and marginal effects should be averaged separately within each model's
#' observed covariate distribution.
#'
#' @param object A `"suest_model"` returned by [suest()].
#'
#' @return A data frame with the component model frames stacked vertically and
#'   two internal columns, `.suest_model` and `.suest_rowid`, used to route rows
#'   to the correct component model.
#'
#' @examples
#' dat <- mtcars
#' dat$am <- factor(dat$am)
#'
#' model1 <- glm(am ~ wt + hp, family = binomial(), data = dat, subset = cyl == 4)
#' model2 <- glm(am ~ wt + hp, family = binomial(), data = dat, subset = cyl != 4)
#' fit <- suest(model1, model2, model_names = c("Four cylinders", "Other"))
#' nd <- suest_newdata(fit)
#' marginaleffects::avg_comparisons(fit, variables = "wt", newdata = nd)
#'
#' @export
suest_newdata <- function(object) {
  if (!inherits(object, "suest_model"))
    stop("'object' must be a suest_model.", call. = FALSE)

  reserved <- c(".suest_model", ".suest_rowid")
  if (any(vapply(
    object$model_frames,
    function(x) any(reserved %in% names(x)),
    logical(1)
  )))
    stop(
      "The model data contain a reserved '.suest_' column name.",
      call. = FALSE
    )

  row_offset <- c(0L, cumsum(vapply(
    object$model_frames,
    nrow,
    integer(1)
  ))[-length(object$model_frames)])

  frames <- lapply(seq_along(object$model_frames), function(i) {
    x <- object$model_frames[[i]]

    for (nm in names(x)) {
      if (is.factor(x[[nm]]))
        x[[nm]] <- as.character(x[[nm]])
    }

    x$.suest_model <- object$model_names[i]
    x$.suest_rowid <- row_offset[i] + seq_len(nrow(x))
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
  class(out) <- c("suest_newdata", class(out))
  out
}

# marginaleffects extension methods
