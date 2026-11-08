# Constant-dispersion negative binomial (Stata's nbreg, dispersion(constant))
# from glmmTMB(family = nbinom1()) without random effects: conditional
# variance mu*(1 + phi), parameters the coefficients and log_phi (Stata's
# lndelta). Scores are analytic; the bread is glmmTMB's covariance (the
# inverse observed information).

.suest_nb1_adapter <- function(model, key) {
  if (!identical(key, "nbinom1"))
    stop(paste0("glmmTMB models without random effects are supported only for the ",
      "constant-dispersion negative binomial, family = nbinom1(); fit other ",
      "single-level models with glm() or MASS::glm.nb()."), call. = FALSE)
  if (length(model$modelInfo$reStruc$ziReStruc) || length(glmmTMB::fixef(model)$zi))
    stop("Zero-inflation and dispersion submodels are not supported for glmmTMB models.",
      call. = FALSE)
  dispersion_terms <- stats::terms(model$modelInfo$allForm$dispformula)
  beta <- glmmTMB::fixef(model)$cond
  log_phi <- model$fit$par[names(model$fit$par) == "betadisp"]
  if (length(attr(dispersion_terms, "term.labels")) || attr(dispersion_terms, "intercept") != 1L ||
      length(attr(dispersion_terms, "offset")) || length(log_phi) != 1L || !is.finite(log_phi))
    stop("glmmTMB nbinom1 requires one estimated, finite constant dispersion (dispformula = ~1).",
      call. = FALSE)
  if (length(model$obj$env$map) ||
      !identical(names(model$fit$par), c(rep("beta", length(beta)), "betadisp")))
    stop("Mapped or constrained glmmTMB nbinom1 parameters are not supported.", call. = FALSE)
  if ("log_phi" %in% names(beta))
    stop("glmmTMB nbinom1 reserves the coefficient name log_phi.", call. = FALSE)
  if (isTRUE(model$modelInfo$REML))
    stop("glmmTMB models must use maximum likelihood.", call. = FALSE)
  if (!is.null(model$modelInfo$priors) && nrow(model$modelInfo$priors))
    stop("glmmTMB models with priors are not supported.", call. = FALSE)
  fit_data <- model$obj$env$data
  if (!is.null(model$call$weights) || length(fit_data$weights) != nrow(model$frame) ||
      any(fit_data$weights != 1))
    stop("Weighted glmmTMB models are not supported.", call. = FALSE)
  if (!is.null(model$call$offset) || length(attr(stats::terms(model), "offset")) ||
      any(fit_data$offset != 0))
    stop("Offsets are not supported for glmmTMB models.", call. = FALSE)
  if (!identical(model$fit$convergence, 0L) || !isTRUE(model$sdr$pdHess))
    stop("The glmmTMB model must converge with a positive-definite Hessian.", call. = FALSE)
  list(engine = "glmmTMB::glmmTMB", type = "nb1")
}

.suest_nb1_parameters <- function(model) {
  if (!is.null(model$suest_parameters)) return(model$suest_parameters)
  c(glmmTMB::fixef(model)$cond, log_phi = unname(glmmTMB::fixef(model)$disp))
}

# Log density and its derivatives in eta and log_phi.
.suest_nb1_scores <- function(y, eta, log_phi) {
  phi <- exp(log_phi)
  k <- exp(eta)/phi
  difference <- digamma(y + k) - digamma(k)
  list(l = stats::dnbinom(y, size = k, mu = exp(eta), log = TRUE),
    eta = k*(difference - log1p(phi)),
    log_phi = -k*difference + k*log1p(phi) - k*phi/(1 + phi) + y/(1 + phi))
}

.suest_nb1_components <- function(model) {
  parameters <- .suest_nb1_parameters(model)
  data <- model$obj$env$data
  X <- as.matrix(data$X)
  y <- as.numeric(data$yobs)
  p <- ncol(X)
  eta <- drop(X %*% parameters[seq_len(p)])
  d <- .suest_nb1_scores(y, eta, parameters[["log_phi"]])
  if (!isTRUE(all.equal(sum(d$l), as.numeric(stats::logLik(model)), tolerance = 1e-8)))
    stop("Could not reproduce the glmmTMB nbinom1 log likelihood for its scores.", call. = FALSE)
  score <- cbind(X*d$eta, log_phi = d$log_phi)
  dimnames(score) <- list(rownames(model$frame), names(parameters))
  covariance <- as.matrix(stats::vcov(model, full = TRUE))
  if (any(dim(covariance) != length(parameters)))
    stop("The glmmTMB nbinom1 covariance does not align with its parameters.", call. = FALSE)
  B <- nrow(X)*covariance
  dimnames(B) <- list(names(parameters), names(parameters))
  if (any(!is.finite(score)) || any(!is.finite(B)))
    stop("Unable to construct finite glmmTMB nbinom1 scores.", call. = FALSE)
  list(score = score, bread = B, parameters = parameters)
}
