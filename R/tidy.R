#' Tidy and glance methods for regression tables
#'
#' `tidy()` returns one row per parameter, with the model name in a `group`
#' column, so `modelsummary::modelsummary(fit, shape = term ~ group)` puts each
#' model in its own column. Standard errors come from the joint covariance
#' matrix. `glance()` returns the number of observations in the combined
#' sample.
#'
#' @param x A `"suest_model"` or `"suest_mi"` object.
#' @param conf.int Include confidence intervals.
#' @param conf.level Confidence level.
#' @param ... Not used.
#' @return A data frame.
#' @examples
#' dat <- mtcars
#' model1 <- glm(am ~ wt, family = binomial(), data = dat)
#' model2 <- glm(vs ~ wt, family = binomial(), data = dat)
#' fit <- suest(model1, model2)
#' if (requireNamespace("generics", quietly = TRUE)) generics::tidy(fit)
#' @name tidy.suest_model
#' @exportS3Method generics::tidy
tidy.suest_model <- function(x, conf.int = TRUE, conf.level = 0.95, ...) {
  estimate <- stats::coef(x)
  std.error <- sqrt(diag(stats::vcov(x)))
  statistic <- estimate/std.error
  out <- data.frame(
    term = unlist(x$local_names, use.names = FALSE),
    group = rep(x$model_names, lengths(x$local_names)),
    estimate = unname(estimate), std.error = unname(std.error),
    statistic = unname(statistic), p.value = unname(2*stats::pnorm(-abs(statistic))),
    row.names = NULL, check.names = FALSE
  )
  if (isTRUE(conf.int)) {
    critical <- stats::qnorm((1 + conf.level)/2)
    out$conf.low <- out$estimate - critical*out$std.error
    out$conf.high <- out$estimate + critical*out$std.error
  }
  out
}

#' @rdname tidy.suest_model
#' @exportS3Method generics::tidy
tidy.suest_mi <- function(x, conf.int = TRUE, conf.level = 0.95, ...) {
  reference <- x$fits[[1L]]
  out <- summary(x, conf.level = conf.level)
  out$term <- unlist(reference$local_names, use.names = FALSE)
  out <- cbind(out[1L], group = rep(x$model_names, lengths(reference$local_names)), out[-1L])
  if (!isTRUE(conf.int)) out$conf.low <- out$conf.high <- NULL
  out
}

#' @rdname tidy.suest_model
#' @exportS3Method generics::glance
glance.suest_model <- function(x, ...) data.frame(nobs = x$nobs_union)

#' @rdname tidy.suest_model
#' @exportS3Method generics::glance
glance.suest_mi <- function(x, ...) data.frame(nobs = x$nobs_union, nimp = x$m)
