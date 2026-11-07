# lme4 models (lmer, glmer). The S4 fit is wrapped in a list so the rest of
# suest can read its call, formula, terms, and frame like other models. Scores
# and observed information come from suest's own per-group likelihood: exact
# for lmer, Laplace or adaptive Gauss-Hermite (the fit's nAGQ) for glmer, so
# they are derivatives of the same likelihood lme4 maximized.

.suest_wrap_lme4 <- function(model) {
  if (!inherits(model, "merMod")) return(model)
  if (!requireNamespace("lme4", quietly = TRUE))
    stop("Package 'lme4' is required for lmer and glmer models.", call. = FALSE)
  frame <- stats::model.frame(model)
  fixed <- stats::terms(model, fixed.only = TRUE)
  formula <- stats::formula(model)
  factors <- vapply(frame, is.factor, logical(1)) &
    names(frame) %in% all.vars(stats::delete.response(fixed))
  structure(list(
    fit = model,
    call = model@call,
    formula = formula,
    terms = fixed,
    frame = frame,
    xlevels = lapply(frame[factors], levels),
    contrasts = attr(lme4::getME(model, "X"), "contrasts"),
    family = stats::family(model),
    engine = if (lme4::isLMM(model)) "lme4::lmer" else "lme4::glmer"
  ), class = "suest_lme4")
}

.suest_lme4_random_info <- function(wrapper) {
  fit <- wrapper$fit
  flist <- lme4::getME(fit, "flist")
  cnms <- lme4::getME(fit, "cnms")
  columns <- if (length(cnms) == 1L) unname(cnms[[1L]]) else character()
  slope <- length(columns) == 2L && identical(columns[1L], "(Intercept)")
  frame <- wrapper$frame
  simple <- length(flist) == 1L && length(cnms) == 1L &&
    (identical(columns, "(Intercept)") || slope) &&
    (!slope || (identical(make.names(columns[2L]), columns[2L]) &&
      columns[2L] %in% names(frame) && is.numeric(frame[[columns[2L]]]) &&
      !is.matrix(frame[[columns[2L]]])))
  if (!simple) {
    structure <- if (!is.null(wrapper$suest_structure)) wrapper$suest_structure else
      .suest_re_structure_lme4(wrapper)
    return(list(name = structure$name, id = structure$id, slope = NULL, structure = structure))
  }
  list(
    name = names(flist)[1L],
    id = droplevels(flist[[1L]]),
    slope = if (slope) columns[2L] else NULL
  )
}

.suest_lme4_adapter <- function(wrapper) {
  fit <- wrapper$fit
  key <- if (identical(wrapper$engine, "lme4::lmer")) {
    "gaussian"
  } else {
    .suest_glmm_key(wrapper$family)
  }
  if (is.na(key) || key %in% c("nbinom2", "gamma"))
    stop(
      paste0(
        "glmer support is currently limited to binomial (logit, probit, ",
        "cloglog) and Poisson-log models; received family '",
        wrapper$family$family, "' with link '", wrapper$family$link, "'."
      ),
      call. = FALSE
    )
  if (identical(wrapper$engine, "lme4::lmer") && lme4::isREML(fit))
    stop("lmer models must be fit by maximum likelihood (REML = FALSE), as Stata's mixed does.",
      call. = FALSE)
  random <- .suest_lme4_random_info(wrapper)
  type <- paste0("glmm_", key,
    if (!is.null(random$structure)) "_re" else if (is.null(random$slope)) "_ri" else "_rs")

  weights <- stats::weights(fit)
  if (!is.null(wrapper$call$weights) || any(weights != 1))
    stop("Weighted lme4 models are not supported.", call. = FALSE)
  offset <- lme4::getME(fit, "offset")
  if (!is.null(wrapper$call$offset) || any(offset != 0))
    stop("Offsets are not supported for lme4 models.", call. = FALSE)
  y <- lme4::getME(fit, "y")
  if (.suest_glmm_binary(type) && !all(y %in% c(0, 1)))
    stop("Grouped-binomial glmer responses are not supported; use Bernoulli 0/1 outcomes.",
      call. = FALSE)
  if (any(names(lme4::fixef(fit)) %in% c("log_sigma", "log_sigma_e", "log_sd_intercept",
      "log_sd_slope", "atanh_rho", random$structure$names)))
    stop("lme4 models reserve the coefficient names log_sigma, log_sigma_e, and the log_sd_ and atanh_rho names of the random effects.",
      call. = FALSE)
  if (any(!is.finite(lme4::getME(fit, "theta"))) || lme4::isSingular(fit))
    stop(
      paste0(
        "The lme4 fit is singular: a random-effect standard deviation is zero ",
        "or a correlation is -1 or 1. suest needs a nonsingular fit."
      ),
      call. = FALSE
    )
  if (!is.null(random$structure)) .suest_re_check_covariances(random$structure)
  list(engine = wrapper$engine, type = type)
}

# Parameters on suest's scale: fixed effects, the log residual SD for lmer,
# then the log random-intercept SD, or the log SDs and the Fisher z of the
# intercept-slope correlation.
.suest_lme4_parameters <- function(wrapper, type) {
  if (!is.null(wrapper$suest_parameters))
    return(wrapper$suest_parameters)
  fit <- wrapper$fit
  beta <- lme4::fixef(fit)
  correlation <- as.data.frame(lme4::VarCorr(fit))
  random <- correlation[correlation$grp != "Residual", , drop = FALSE]
  residual <- if (identical(wrapper$engine, "lme4::lmer"))
    c(log_sigma_e = log(stats::sigma(fit))) else numeric(0)
  if (endsWith(type, "_re")) {
    structure <- .suest_lme4_random_info(wrapper)$structure
    return(c(beta, residual, stats::setNames(structure$values, structure$names)))
  }
  if (.suest_glmm_slope(type)) {
    sd <- random$sdcor[is.na(random$var2)]
    rho <- random$sdcor[!is.na(random$var2)]
    return(c(beta, residual, log_sd_intercept = log(sd[1L]),
      log_sd_slope = log(sd[2L]), atanh_rho = atanh(rho)))
  }
  c(beta, residual, log_sigma = log(random$sdcor[1L]))
}

# Log density of each observation and its first three derivatives in eta.
# dispersion: phi for nbinom2 (variance mu + mu^2/phi), the shape for Gamma.
.suest_glmm_density <- function(family, y, eta, dispersion = NULL) {
  if (family == "logit") {
    p <- stats::plogis(eta)
    v <- p*(1 - p)
    list(l = stats::dbinom(y, 1, p, log = TRUE), d1 = y - p, d2 = -v, d3 = -v*(1 - 2*p))
  } else if (family == "probit") {
    lambda1 <- exp(stats::dnorm(eta, log = TRUE) - stats::pnorm(eta, log.p = TRUE))
    lambda0 <- exp(stats::dnorm(eta, log = TRUE) - stats::pnorm(-eta, log.p = TRUE))
    curve1 <- -lambda1*(eta + lambda1)
    curve0 <- lambda0*(lambda0 - eta)
    list(l = ifelse(y == 1, stats::pnorm(eta, log.p = TRUE), stats::pnorm(-eta, log.p = TRUE)),
      d1 = ifelse(y == 1, lambda1, -lambda0),
      d2 = ifelse(y == 1, curve1, -curve0),
      d3 = ifelse(y == 1, -curve1*(eta + 2*lambda1) - lambda1,
        lambda0 - curve0*(2*lambda0 - eta)))
  } else if (family == "cloglog") {
    t <- exp(eta)
    # With y = 1 and t > 700 the derivatives are below 1e-300.
    tail <- y == 1 & t > 700
    s <- ifelse(tail, 1, t)
    inverse <- 1/expm1(s)
    r <- -1/expm1(-s)
    f1 <- inverse*(1 - s*r)
    f2 <- inverse*(s*r*(2*r - 1) - 2*r)
    list(l = ifelse(y == 1, log(-expm1(-t)), -t),
      d1 = ifelse(tail, 0, ifelse(y == 1, s*inverse, -t)),
      d2 = ifelse(tail, 0, ifelse(y == 1, s*f1, -t)),
      d3 = ifelse(tail, 0, ifelse(y == 1, s*f1 + s^2*f2, -t)))
  } else if (family == "poisson") {
    mu <- exp(eta)
    list(l = stats::dpois(y, mu, log = TRUE), d1 = y - mu, d2 = -mu, d3 = -mu)
  } else if (family == "nbinom2") {
    mu <- exp(eta)
    phi <- dispersion
    list(l = stats::dnbinom(y, size = phi, mu = mu, log = TRUE),
      d1 = phi*(y - mu)/(mu + phi),
      d2 = -phi*mu*(y + phi)/(mu + phi)^2,
      d3 = -phi*mu*(y + phi)*(phi - mu)/(mu + phi)^3)
  } else {
    shape <- dispersion
    ratio <- y*exp(-eta)
    list(l = stats::dgamma(y, shape = shape, rate = shape*exp(-eta), log = TRUE),
      d1 = shape*(ratio - 1), d2 = -shape*ratio, d3 = shape*ratio)
  }
}

# Random-effect covariance from the log SDs (and the Fisher z of the
# intercept-slope correlation).
.suest_re_covariance <- function(theta, slope) {
  if (!slope) return(matrix(exp(2*theta[["log_sigma"]])))
  sd <- exp(c(theta[["log_sd_intercept"]], theta[["log_sd_slope"]]))
  rho <- tanh(theta[["atanh_rho"]])
  matrix(c(sd[1L]^2, rho*sd[1L]*sd[2L], rho*sd[1L]*sd[2L], sd[2L]^2), 2L)
}

# Per-group log likelihood of the suest parameter vector, for all groups at
# once; random-intercept glmer models also return exact group scores.
.suest_lme4_group_loglik <- function(wrapper, type) {
  fit <- wrapper$fit
  X <- as.matrix(lme4::getME(fit, "X"))
  y <- as.numeric(lme4::getME(fit, "y"))
  random <- .suest_lme4_random_info(wrapper)
  p <- ncol(X)
  gaussian <- identical(wrapper$engine, "lme4::lmer")
  family <- .suest_glmm_family(type)
  if (!is.null(random$structure))
    return(.suest_re_general_loglik(random$structure, wrapper$frame, X, y,
      if (gaussian) "gaussian" else family, p + gaussian))
  slope <- if (is.null(random$slope)) NULL else wrapper$frame[[random$slope]]
  nagq <- unname(fit@devcomp$dims["nAGQ"])
  rule <- if (!gaussian && nagq > 1L) {
    hermite <- lme4::GHrule(nagq)
    list(z = hermite[, "z"], w = hermite[, "w"])
  }
  engine <- .suest_re_engine(random$id, X, slope, rule)
  density <- function(eta) .suest_glmm_density(family, y, eta)

  function(parameters, start = NULL, scores = FALSE) {
    eta <- drop(X %*% parameters[seq_len(p)])
    Sigma <- .suest_re_covariance(parameters, !is.null(slope))
    if (gaussian)
      return(list(loglik = engine$gaussian(y, eta, Sigma,
        exp(2*parameters[["log_sigma_e"]]))))
    if (is.null(slope))
      return(engine$intercept(eta, Sigma[1L, 1L], density, start, scores))
    engine$intercept_slope(eta, Sigma, density, start)
  }
}

.suest_lme4_components <- function(wrapper, type) {
  parameters <- .suest_lme4_parameters(wrapper, type)
  evaluate <- .suest_lme4_group_loglik(wrapper, type)
  # lme4's default Laplace fits stop the inner mode search early (tolPwrss),
  # so their log likelihood can differ slightly from the exact Laplace value.
  glmer <- identical(wrapper$engine, "lme4::glmer")
  nagq <- unname(wrapper$fit@devcomp$dims["nAGQ"])
  laplace <- glmer && nagq == 1L
  # For nAGQ > 1 Poisson fits, lme4 1.1 reports the log likelihood less that of
  # the saturated model and lme4 2.1 the full value; accept either.
  y <- as.numeric(lme4::getME(wrapper$fit, "y"))
  saturated <- if (glmer && nagq > 1L && .suest_glmm_family(type) == "poisson")
    sum(stats::dpois(y, y, log = TRUE)) else 0
  base <- evaluate(parameters, scores = TRUE)
  reported <- as.numeric(stats::logLik(wrapper$fit)) + unique(c(0, saturated))
  if (!any(vapply(reported, function(value) isTRUE(all.equal(sum(base$loglik), value,
      tolerance = if (laplace) 1e-4 else 1e-6)), logical(1))))
    stop("Could not reproduce the lme4 log likelihood for its scores.", call. = FALSE)
  .suest_re_components(evaluate, parameters, base,
    .suest_lme4_random_info(wrapper)$id, rownames(wrapper$frame))
}

.suest_lme4_frame <- function(wrapper) {
  frame <- as.data.frame(wrapper$frame)
  original <- try(eval(wrapper$call$data, envir = environment(wrapper$formula)),
    silent = TRUE)
  if (!inherits(original, "try-error") && is.data.frame(original)) {
    index <- match(rownames(frame), rownames(original))
    if (!anyNA(index)) {
      for (variable in setdiff(names(original), names(frame)))
        frame[[variable]] <- original[[variable]][index]
    }
  }
  frame
}

.suest_lme4_eta <- function(wrapper, parameters, newdata) {
  terms <- stats::delete.response(wrapper$terms)
  frame <- stats::model.frame(terms, newdata, na.action = stats::na.pass,
    xlev = wrapper$xlevels)
  X <- stats::model.matrix(terms, frame, contrasts.arg = wrapper$contrasts)
  beta_names <- names(lme4::fixef(wrapper$fit))
  if (!all(beta_names %in% colnames(X)))
    stop("Unable to construct the lme4 prediction matrix.", call. = FALSE)
  as.numeric(X[, beta_names, drop = FALSE] %*% parameters[beta_names])
}
