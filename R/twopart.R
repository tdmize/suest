#' Two-part models
#'
#' `suest_twopart()` combines the two parts of a two-part model (Stata's
#' `twopm`): a binary model for whether the outcome is positive, fit to
#' everyone, and a model for the outcome among those with a positive outcome.
#' Predictions, slopes, and comparisons from `marginaleffects` are then on the
#' overall expected outcome, `Pr(y > 0) * E(y | y > 0)`, with standard errors
#' from the joint covariance matrix of the two parts.
#'
#' @param first Binary logit, probit, or complementary log-log model for
#'   whether the outcome is positive.
#' @param second Model for the outcome among the observations with a positive
#'   outcome: a linear model or a GLM, for example a Gamma model or a Poisson
#'   pseudo-likelihood model (`quasipoisson()`) with a log link.
#' @param model_names Names of the two parts. By default, the object names
#'   supplied in the call.
#' @param ... Passed to [suest()], for example `cluster`, `weight_type`, or
#'   `survey_design` with `observation_id`.
#' @return An object of class `c("suest_twopart", "suest_model")`.
#' @details
#' Fit the second part to every observation in the first part's sample that
#' has a positive outcome, for example with `subset = y > 0`. Without
#' `newdata`, results are averaged over the first part's estimation sample.
#' `type = "response"` (the default) gives the overall expected outcome,
#' `type = "prob"` the probability of a positive outcome (the first part), and
#' `type = "conditional"` the expected outcome among those with a positive
#' outcome (the second part).
#' @examples
#' dat <- mtcars
#' dat$y <- ifelse(dat$am == 1, dat$mpg, 0)
#' any <- glm(I(y > 0) ~ wt, family = binomial(), data = dat)
#' amount <- glm(y ~ wt, family = Gamma(link = "log"), data = dat, subset = y > 0)
#' fit <- suest_twopart(any, amount)
#' marginaleffects::avg_slopes(fit, variables = "wt")
#' @export
suest_twopart <- function(first, second, model_names = NULL, ...) {
  if (is.null(model_names))
    model_names <- c(deparse1(substitute(first)), deparse1(substitute(second)))
  fit <- suest(first, second, model_names = model_names, ...)
  if (!fit$model_types[1L] %in% c("logit", "probit", "cloglog"))
    stop("The first part must be a binary logit, probit, or complementary log-log model.",
      call. = FALSE)
  if (!fit$model_types[2L] %in% c("lm", "survey_lm", "glm", "poisson"))
    stop("The second part must be a linear model or a GLM.", call. = FALSE)
  y <- as.numeric(fit$models[[1L]]$y)
  keys <- fit$sample_keys
  if (length(y) != length(keys[[1L]]) || !all(y %in% c(0, 1)) ||
      !setequal(keys[[2L]], keys[[1L]][y == 1]))
    stop(paste0("Fit the second part to the observations in the first part's sample ",
      "with a positive outcome (for example with subset = y > 0). If the parts use ",
      "different data objects, identify the observations with observation_id."), call. = FALSE)
  class(fit) <- c("suest_twopart", class(fit))
  if (anyNA(get_predict.suest_twopart(fit, .suest_twopart_data(fit))$estimate))
    stop(paste0("Some observations in the first part's sample lack a predictor of the ",
      "second part. Fit both parts to the same observations."), call. = FALSE)
  fit
}

# The first part's estimation sample, with every variable of the data object.
.suest_twopart_data <- function(x) {
  first <- x$models[[1L]]
  rows <- rownames(x$model_frames[[1L]])
  data <- if (inherits(first, "svyglm")) first$survey.design$variables else
    tryCatch(eval(first$call$data, environment(stats::formula(first))),
      error = function(e) NULL)
  index <- if (is.data.frame(data)) match(rows, rownames(data))
  if (is.null(index) || anyNA(index)) {
    out <- suest_newdata(x)
    out <- out[out$.suest_model == x$model_names[1L], , drop = FALSE]
    out$.suest_model <- out$.suest_rowid <- NULL
  } else {
    out <- as.data.frame(data)[index, , drop = FALSE]
    if (any(x$weight_type %in% c("pweight", "survey")))
      out$.suest_weight <- x$model_weights[[1L]]
  }
  rownames(out) <- NULL
  class(out) <- "data.frame"
  out
}

get_data.suest_twopart <- function(x, ...) .suest_twopart_data(x)

get_predict.suest_twopart <- function(model, newdata, type = "response", ...) {
  if (is.null(type)) type <- "response"
  if (!type %in% c("response", "prob", "conditional"))
    stop("'type' must be 'response', 'prob', or 'conditional'.", call. = FALSE)
  if (".suest_model" %in% names(newdata))
    stop(paste0("Two-part models use one row per observation: leave newdata out to ",
      "average over the first part's sample, or supply data without .suest_model."),
      call. = FALSE)
  p <- get_predict.suest_model(model, newdata, type = "response")
  n <- nrow(newdata)
  first <- p$estimate[seq_len(n)]
  second <- p$estimate[n + seq_len(n)]
  data.frame(
    rowid = p$rowid[seq_len(n)],
    estimate = switch(type, response = first*second, prob = first, conditional = second)
  )
}

#' @rdname suest_twopart
#' @param x A `"suest_twopart"` object.
#' @export
print.suest_twopart <- function(x, ...) {
  cat("Two-part model: first part ", x$model_names[1L], ", second part ",
    x$model_names[2L], "\n", sep = "")
  NextMethod()
}
