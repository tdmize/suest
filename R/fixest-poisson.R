# Fixed-effects Poisson from fixest::fepois(), Stata's xtpoisson, fe. The
# concentrated likelihood gives the same coefficients, scores, and information
# as the conditional likelihood. Each unit's fixed effect is a function of the
# coefficients, alpha_i(b) = log(sum y_i) - log(sum exp(x_i b)), so predictions
# and their uncertainty follow from the coefficients alone.

.suest_fepois_adapter <- function(model) {
  fe <- model$fixef_vars
  if (length(fe) != 1L || grepl("[\\^\\[]", fe))
    stop("Fixed-effects Poisson models must absorb exactly one plain fixed effect (for example, | id).",
      call. = FALSE)
  if (!identical(model$family$family, "poisson") || !identical(model$family$link, "log"))
    stop("Only Poisson-log fixest::fepois() models are supported.", call. = FALSE)
  if (!is.null(model$weights) || !is.null(model$call$weights))
    stop("Weighted fixest::fepois() models are not supported.", call. = FALSE)
  if (!is.null(model$offset) || !is.null(model$call$offset))
    stop("Offsets are not supported for fixest::fepois() models.", call. = FALSE)
  if (isTRUE(model$is_iv) || !isTRUE(model$convStatus))
    stop("The fixest::fepois() model must converge and must not be an IV model.", call. = FALSE)
  list(engine = "fixest::fepois", type = "panel_poisson_fe")
}

.suest_fepois_data <- function(model) {
  environment <- if (is.environment(model$call_env)) model$call_env else
    environment(stats::formula(model))
  data <- try(eval(model$call$data, envir = environment), silent = TRUE)
  if (inherits(data, "try-error") || !is.data.frame(data))
    stop("Could not recover the original data for the fixest::fepois() model; keep its data object available.",
      call. = FALSE)
  rows <- fixest::obs(model)
  if (length(rows) != stats::nobs(model) || anyNA(rows) || any(rows < 1L | rows > nrow(data)))
    stop("Could not align the fixest::fepois() estimation sample to its data.", call. = FALSE)
  data[rows, , drop = FALSE]
}

.suest_fepois_frame <- function(model) {
  estimation_data <- .suest_fepois_data(model)
  frame <- stats::model.frame(model$fml_all$linear, data = estimation_data,
    na.action = stats::na.pass)
  for (variable in setdiff(names(estimation_data), names(frame)))
    frame[[variable]] <- estimation_data[[variable]]
  rownames(frame) <- rownames(estimation_data)
  frame
}

.suest_fepois_id <- function(model) {
  id <- model$fixef_id[[1L]]
  factor(attr(id, "fixef_names")[id], levels = attr(id, "fixef_names"))
}

.suest_fepois_components <- function(model) {
  parameters <- stats::coef(model)
  score <- as.matrix(sandwich::estfun(model))
  # Inverse information without fixest's small-sample factor, which counts the
  # absorbed fixed effects.
  covariance <- unclass(stats::vcov(model, vcov = "iid",
    ssc = fixest::ssc(adj = FALSE, cluster.adj = FALSE)))
  covariance <- matrix(covariance, nrow(covariance),
    dimnames = list(names(parameters), names(parameters)))
  n <- nrow(score)
  if (is.null(colnames(score))) colnames(score) <- names(parameters)
  if (!identical(colnames(score), names(parameters)) ||
      any(dim(covariance) != length(parameters)) || n != stats::nobs(model))
    stop("The fixest::fepois() scores and covariance do not align.", call. = FALSE)
  rownames(score) <- rownames(.suest_fepois_frame(model))
  list(score = score, bread = n*covariance, parameters = parameters)
}

.suest_fepois_eta <- function(model, newdata) {
  beta <- stats::coef(model)
  X <- stats::model.matrix(model, type = "rhs")
  id <- as.character(.suest_fepois_id(model))
  y <- as.numeric(model$y)
  alpha <- log(tapply(y, id, sum)) - log(tapply(exp(drop(X %*% beta)), id, sum))
  fe <- model$fixef_vars
  if (!fe %in% names(newdata))
    stop(sprintf("Fixed-effects Poisson predictions need the '%s' column in newdata.", fe),
      call. = FALSE)
  X_new <- stats::model.matrix(model, data = newdata, type = "rhs")
  drop(X_new[, names(beta), drop = FALSE] %*% beta) +
    unname(alpha[as.character(newdata[[fe]])])
}
