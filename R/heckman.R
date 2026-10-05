# Heckman selection models fit by maximum likelihood with
# sampleSelection::selection(method = "ml"), Stata's heckman. Ancillary
# parameters use Stata's scale: lnsigma = log(sigma), athrho = atanh(rho).
# Predictions are the outcome equation's linear prediction, Stata's default.

.suest_is_heckman <- function(model)
  inherits(model, "selection") && identical(model$method, "ml") &&
    identical(as.integer(model$tobitType), 2L)

# Give suest's copy a formula so data-source and cluster lookups find the data.
.suest_prepare_heckman <- function(model) {
  if (inherits(model, "selection") && is.list(model) && is.null(model$formula) &&
      !is.null(model$termsS))
    model$formula <- stats::formula(model$termsS)
  model
}

.suest_heckman_adapter <- function(model) {
  if (!.suest_is_heckman(model))
    stop("Only maximum-likelihood type 2 (Heckman) sampleSelection::selection() models are supported.",
      call. = FALSE)
  if (!is.null(model$call$weights))
    stop("Weighted selection models are not supported.", call. = FALSE)
  if (!is.matrix(model$gradientObs) || !is.matrix(model$hessian))
    stop("The selection model does not retain its observation scores and Hessian.", call. = FALSE)
  list(engine = "sampleSelection::selection", type = "heckman")
}

.suest_heckman_data <- function(model) {
  formula <- stats::formula(model$termsS)
  data <- try(eval(model$call$data, envir = environment(formula)), silent = TRUE)
  if (inherits(data, "try-error") || !is.data.frame(data))
    stop("Could not recover the original data for the selection model; keep its data object available.",
      call. = FALSE)
  data
}

# Rows used: complete selection variables, and complete outcome variables
# where the observation is selected, as sampleSelection does.
.suest_heckman_frame <- function(model) {
  data <- .suest_heckman_data(model)
  selection <- stats::model.frame(model$termsS, data, na.action = stats::na.pass)
  selected <- as.logical(stats::model.response(selection))
  outcome_vars <- all.vars(model$termsO)
  complete_selection <- stats::complete.cases(selection)
  complete_outcome <- stats::complete.cases(data[outcome_vars])
  keep <- complete_selection & (!selected | complete_outcome)
  keep[is.na(keep)] <- FALSE
  frame <- data[keep, , drop = FALSE]
  attr(frame, "terms") <- model$termsS
  if (nrow(frame) != nrow(model$gradientObs))
    stop("Could not align the selection model's estimation sample to its data.", call. = FALSE)
  frame
}

.suest_heckman_jacobian <- function(model) {
  raw <- as.numeric(stats::coef(model))
  index <- model$param$index
  k <- length(raw)
  derivative <- rep(1, k)
  derivative[index$errTerms["sigma"]] <- raw[index$errTerms["sigma"]]
  derivative[index$errTerms["rho"]] <- 1 - raw[index$errTerms["rho"]]^2
  derivative
}

.suest_heckman_parameters <- function(model) {
  if (!is.null(model$suest_parameters)) return(model$suest_parameters)
  raw <- stats::setNames(as.numeric(stats::coef(model)), names(stats::coef(model)))
  index <- model$param$index
  out <- raw
  out[index$errTerms["sigma"]] <- log(raw[index$errTerms["sigma"]])
  out[index$errTerms["rho"]] <- atanh(raw[index$errTerms["rho"]])
  names(out)[index$betaS] <- paste0("selection:", names(raw)[index$betaS])
  names(out)[index$betaO] <- paste0("outcome:", names(raw)[index$betaO])
  names(out)[index$errTerms] <- c("lnsigma", "athrho")
  out
}

.suest_heckman_components <- function(model) {
  parameters <- .suest_heckman_parameters(model)
  derivative <- .suest_heckman_jacobian(model)
  score <- sweep(model$gradientObs, 2L, derivative, `*`)
  information <- -model$hessian*outer(derivative, derivative)
  information <- (information + t(information))/2
  frame <- .suest_heckman_frame(model)
  dimnames(score) <- list(rownames(frame), names(parameters))
  B <- nrow(score)*.suest_information_inverse(information)
  dimnames(B) <- list(names(parameters), names(parameters))
  list(score = score, bread = B, parameters = parameters)
}

.suest_heckman_predict <- function(model, newdata) {
  parameters <- .suest_heckman_parameters(model)
  terms <- stats::delete.response(model$termsO)
  frame <- stats::model.frame(terms, newdata, na.action = stats::na.pass)
  X <- stats::model.matrix(terms, frame)
  beta <- parameters[model$param$index$betaO]
  drop(X %*% beta)
}
