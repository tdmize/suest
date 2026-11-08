# Bivariate probit with a different equation for each outcome (Stata's
# biprobit (y1 = ...) (y2 = ...)) from GJRM::gjrm(margins = c("probit",
# "probit"), model = "B") with the default Gaussian copula. Parameters are
# each equation's coefficients ("y1:x", "y2:w") and athrho = atanh(rho),
# which is GJRM's theta.star.

# gjrm() keeps no call to the user's data, so look for the one data frame in
# the formula's environment that holds the estimation sample.
.suest_gjrm_find_data <- function(model, rows) {
  env <- environment(model$formula[[1L]])
  if (!is.environment(env)) return(NULL)
  frames <- list(model$gam1$model, model$gam2$model)
  variables <- unique(unlist(lapply(frames, function(f)
    names(f)[names(f) == make.names(names(f))])))
  found <- character()
  for (name in ls(env)) {
    object <- get(name, envir = env)
    if (!is.data.frame(object) || !all(variables %in% names(object)) ||
        !all(rows %in% rownames(object))) next
    same <- all(vapply(frames, function(f) {
      shared <- intersect(names(f), variables)
      index <- match(rownames(f), rownames(object))
      all(vapply(shared, function(v) identical(as.character(f[[v]]),
        as.character(object[[v]][index])), logical(1)))
    }, logical(1)))
    if (same) found <- c(found, name)
  }
  if (length(found) != 1L) return(NULL)
  list(name = found, data = get(found, envir = env))
}

.suest_wrap_gjrm <- function(model) {
  if (!inherits(model, "gjrm") || inherits(model, "suest_gjrm")) return(model)
  rows <- rownames(model$X1)
  data <- .suest_gjrm_find_data(model, rows)
  if (is.null(data)) {
    source <- model$gam1$model
    for (v in setdiff(names(model$gam2$model), names(source)))
      source[[v]] <- model$gam2$model[[v]][match(rownames(source), rownames(model$gam2$model))]
    call <- call("gjrm")
  } else {
    source <- data$data[match(rows, rownames(data$data)), , drop = FALSE]
    call <- call("gjrm", data = as.name(data$name))
  }
  responses <- vapply(model$formula[1:2], function(f) all.vars(f[[2L]])[1L], character(1))
  # One frame for both equations: the outcome pair, every predictor, and the
  # other columns of the data (for observation_id and cluster).
  labels <- unique(unlist(lapply(model$formula[1:2], function(f)
    attr(stats::terms(f), "term.labels"))))
  joint <- stats::as.formula(paste0("cbind(", responses[1L], ", ", responses[2L], ") ~ ",
    if (length(labels)) paste(labels, collapse = " + ") else "1"),
    env = environment(model$formula[[1L]]))
  frame <- tryCatch(stats::model.frame(joint, data = source, na.action = stats::na.pass),
    error = function(e) NULL)
  if (is.null(frame)) frame <- source
  for (v in setdiff(names(source), names(frame))) frame[[v]] <- source[[v]]
  rownames(frame) <- rows
  structure(list(
    fit = model,
    call = call,
    formula = model$formula[[1L]],
    formulas = model$formula,
    frame = frame,
    responses = responses
  ), class = "suest_gjrm")
}

.suest_gjrm_adapter <- function(wrapper) {
  model <- wrapper$fit
  supported <- identical(model$Model, "B") &&
    identical(as.character(model$margins), c("probit", "probit")) &&
    identical(model$BivD, "N")
  if (!supported)
    stop(paste0("GJRM::gjrm() models are supported as bivariate probits: ",
      "margins = c(\"probit\", \"probit\"), model = \"B\", and the default ",
      "Gaussian copula."), call. = FALSE)
  if (length(model$formula) != 2L || !is.null(model$X3))
    stop("GJRM bivariate probits with covariates in the correlation are not supported.",
      call. = FALSE)
  if (any(c(model$l.sp1, model$l.sp2, model$l.sp3) != 0))
    stop("GJRM models with smooth terms are not supported.", call. = FALSE)
  if (any(model$weights != 1))
    stop("Weighted GJRM models are not supported.", call. = FALSE)
  if (any(model$gam1$offset != 0) || any(model$gam2$offset != 0))
    stop("Offsets are not supported for GJRM models.", call. = FALSE)
  gradient <- model$fit$gradient
  hessian <- model$fit$hessian
  if (!is.numeric(gradient) || max(abs(gradient)) > 1e-3 || !is.matrix(hessian) ||
      min(eigen(hessian, symmetric = TRUE, only.values = TRUE)$values) <= 0)
    stop("The GJRM model must converge with a positive-definite Hessian.", call. = FALSE)
  if (!all(c(model$y1, model$y2) %in% c(0, 1)))
    stop("Bivariate-probit outcomes must be numeric 0/1.", call. = FALSE)
  if (!requireNamespace("mvtnorm", quietly = TRUE))
    stop("Package 'mvtnorm' is required for bivariate probit models.", call. = FALSE)
  list(engine = "GJRM::gjrm", type = "biprobit_gjrm")
}

.suest_gjrm_parameters <- function(wrapper) {
  if (!is.null(wrapper$suest_parameters)) return(wrapper$suest_parameters)
  model <- wrapper$fit
  b <- as.numeric(model$coefficients)
  k1 <- ncol(model$X1)
  k2 <- ncol(model$X2)
  if (length(b) != k1 + k2 + 1L)
    stop("Unable to align the GJRM coefficients with its equations.", call. = FALSE)
  names(b) <- c(paste0(wrapper$responses[1L], ":", colnames(model$X1)),
    paste0(wrapper$responses[2L], ":", colnames(model$X2)), "athrho")
  b
}

# Observation scores of the bivariate probit in (b1, b2, athrho).
.suest_biprobit_observation_scores <- function(b, X1, X2, y1, y2) {
  k1 <- ncol(X1)
  k2 <- ncol(X2)
  eta1 <- as.numeric(X1 %*% b[seq_len(k1)])
  eta2 <- as.numeric(X2 %*% b[k1 + seq_len(k2)])
  rho <- tanh(b[k1 + k2 + 1L])
  q1 <- 2*y1 - 1
  q2 <- 2*y2 - 1
  scale <- sqrt(1 - rho^2)
  probability <- vapply(seq_along(eta1), function(i) {
    r <- q1[i]*q2[i]*rho
    as.numeric(mvtnorm::pmvnorm(upper = c(q1[i]*eta1[i], q2[i]*eta2[i]),
      corr = matrix(c(1, r, r, 1), 2), algorithm = mvtnorm::TVPACK()))
  }, numeric(1))
  first <- q1*stats::dnorm(eta1)*stats::pnorm(q2*(eta2 - rho*eta1)/scale)/probability
  second <- q2*stats::dnorm(eta2)*stats::pnorm(q1*(eta1 - rho*eta2)/scale)/probability
  density <- exp(-(eta1^2 - 2*rho*eta1*eta2 + eta2^2)/(2*scale^2))/(2*pi*scale)
  cbind(X1*first, X2*second, athrho = q1*q2*density*scale^2/probability)
}

.suest_gjrm_components <- function(wrapper) {
  model <- wrapper$fit
  theta <- .suest_gjrm_parameters(wrapper)
  if (abs(tanh(theta[["athrho"]])) >= 1)
    stop("The bivariate-probit correlation must be strictly between -1 and 1.", call. = FALSE)
  score <- function(b) .suest_biprobit_observation_scores(b, model$X1, model$X2,
    model$y1, model$y2)
  U <- score(theta)
  information <- .suest_score_information(theta, function(b) colSums(score(b)))
  B <- nrow(U)*.suest_information_inverse(information)
  dimnames(U) <- list(rownames(model$X1), names(theta))
  dimnames(B) <- list(names(theta), names(theta))
  if (any(!is.finite(U)) || any(!is.finite(B)))
    stop("Unable to construct finite bivariate-probit scores.", call. = FALSE)
  list(score = U, bread = B, parameters = theta)
}

# Joint probability that both outcomes equal 1, as Stata's biprobit default.
.suest_gjrm_predict <- function(wrapper, newdata, type) {
  if (identical(type, "link"))
    stop("Bivariate probit models are supported on the response scale only.", call. = FALSE)
  model <- wrapper$fit
  theta <- .suest_gjrm_parameters(wrapper)
  design <- function(gam, X) {
    terms <- stats::delete.response(stats::terms(gam$pterms))
    frame <- stats::model.frame(terms, newdata, na.action = stats::na.pass,
      xlev = gam$xlevels)
    out <- stats::model.matrix(terms, frame, contrasts.arg = gam$contrasts)
    if (!identical(colnames(out), colnames(X)))
      stop("Unable to construct the bivariate-probit prediction matrix.", call. = FALSE)
    out
  }
  X1 <- design(model$gam1, model$X1)
  X2 <- design(model$gam2, model$X2)
  k1 <- ncol(X1)
  eta1 <- as.numeric(X1 %*% theta[seq_len(k1)])
  eta2 <- as.numeric(X2 %*% theta[k1 + seq_len(ncol(X2))])
  rho <- tanh(theta[["athrho"]])
  sigma <- matrix(c(1, rho, rho, 1), 2)
  vapply(seq_along(eta1), function(i) {
    if (!is.finite(eta1[i]) || !is.finite(eta2[i])) return(NA_real_)
    as.numeric(mvtnorm::pmvnorm(upper = c(eta1[i], eta2[i]), sigma = sigma,
      algorithm = mvtnorm::TVPACK()))
  }, numeric(1))
}
