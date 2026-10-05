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
  if (length(flist) != 1L || length(cnms) != 1L ||
      (!identical(columns, "(Intercept)") && !slope))
    stop(
      paste0(
        "lme4 support is currently limited to one grouping variable with a ",
        "random intercept, or a random intercept and one correlated numeric slope."
      ),
      call. = FALSE
    )
  frame <- wrapper$frame
  if (slope && (!identical(make.names(columns[2L]), columns[2L]) ||
      !columns[2L] %in% names(frame) || !is.numeric(frame[[columns[2L]]]) ||
      is.matrix(frame[[columns[2L]]])))
    stop("The random slope must be one untransformed numeric data column.", call. = FALSE)
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
  type <- paste0("glmm_", key, if (is.null(random$slope)) "_ri" else "_rs")

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
  if (any(names(lme4::fixef(fit)) %in%
      c("log_sigma", "log_sigma_e", "log_sd_intercept", "log_sd_slope", "atanh_rho")))
    stop("lme4 models reserve the coefficient names log_sigma, log_sigma_e, log_sd_intercept, log_sd_slope, and atanh_rho.",
      call. = FALSE)
  theta <- lme4::getME(fit, "theta")
  if (any(!is.finite(theta)) || any(abs(theta[lme4::getME(fit, "lower") == 0]) <=
      sqrt(.Machine$double.eps)))
    stop(
      paste0(
        "The lme4 random-effect standard deviations must be finite and ",
        "strictly away from zero (the fit is singular)."
      ),
      call. = FALSE
    )
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
  if (.suest_glmm_slope(type)) {
    sd <- random$sdcor[is.na(random$var2)]
    rho <- random$sdcor[!is.na(random$var2)]
    return(c(beta, residual, log_sd_intercept = log(sd[1L]),
      log_sd_slope = log(sd[2L]), atanh_rho = atanh(rho)))
  }
  c(beta, residual, log_sigma = log(random$sdcor[1L]))
}

# Log density of each observation and its first two derivatives in eta.
.suest_glmm_density <- function(family, y, eta) {
  if (family == "logit") {
    p <- stats::plogis(eta)
    list(l = stats::dbinom(y, 1, p, log = TRUE), d1 = y - p, d2 = -p*(1 - p))
  } else if (family == "probit") {
    lambda1 <- exp(stats::dnorm(eta, log = TRUE) - stats::pnorm(eta, log.p = TRUE))
    lambda0 <- exp(stats::dnorm(eta, log = TRUE) - stats::pnorm(-eta, log.p = TRUE))
    list(l = ifelse(y == 1, stats::pnorm(eta, log.p = TRUE), stats::pnorm(-eta, log.p = TRUE)),
      d1 = ifelse(y == 1, lambda1, -lambda0),
      d2 = ifelse(y == 1, -lambda1*(eta + lambda1), -lambda0*(lambda0 - eta)))
  } else if (family == "cloglog") {
    t <- exp(eta)
    list(l = ifelse(y == 1, log(-expm1(-t)), -t),
      d1 = ifelse(y == 1, t/expm1(t), -t),
      d2 = ifelse(y == 1, t*(expm1(t) - t*exp(t))/expm1(t)^2, -t))
  } else {
    mu <- exp(eta)
    list(l = stats::dpois(y, mu, log = TRUE), d1 = y - mu, d2 = -mu)
  }
}

# Per-group log-likelihood function of the suest parameter vector.
.suest_lme4_group_loglik <- function(wrapper, type) {
  fit <- wrapper$fit
  X <- lme4::getME(fit, "X")
  y <- as.numeric(lme4::getME(fit, "y"))
  random <- .suest_lme4_random_info(wrapper)
  group <- random$id
  rows <- split(seq_along(y), group)
  slope <- if (is.null(random$slope)) NULL else wrapper$frame[[random$slope]]
  p <- ncol(X)
  gaussian <- identical(wrapper$engine, "lme4::lmer")
  family <- .suest_glmm_family(type)
  nagq <- unname(fit@devcomp$dims["nAGQ"])
  rule <- if (!gaussian && nagq > 1L) lme4::GHrule(nagq) else NULL

  random_covariance <- function(theta) {
    if (is.null(slope)) return(matrix(exp(2*theta[["log_sigma"]])))
    sd <- exp(c(theta[["log_sd_intercept"]], theta[["log_sd_slope"]]))
    rho <- tanh(theta[["atanh_rho"]])
    matrix(c(sd[1L]^2, rho*sd[1L]*sd[2L], rho*sd[1L]*sd[2L], sd[2L]^2), 2L)
  }

  function(parameters, groups = seq_along(rows)) {
    beta <- parameters[seq_len(p)]
    eta_all <- drop(X %*% beta)
    Sigma <- random_covariance(parameters)
    vapply(groups, function(g) {
      i <- rows[[g]]
      Z <- if (is.null(slope)) matrix(1, length(i), 1L) else cbind(1, slope[i])
      eta <- eta_all[i]
      if (gaussian) {
        V <- exp(2*parameters[["log_sigma_e"]])*diag(length(i)) + Z %*% Sigma %*% t(Z)
        R <- chol(V)
        r <- backsolve(R, y[i] - eta, transpose = TRUE)
        return(-0.5*(length(i)*log(2*pi) + 2*sum(log(diag(R))) + sum(r^2)))
      }
      precision <- solve(Sigma)
      q <- ncol(Z)
      h <- function(b) {
        part <- .suest_glmm_density(family, y[i], eta + drop(Z %*% b))
        list(value = sum(part$l) - 0.5*drop(t(b) %*% precision %*% b) -
            0.5*as.numeric(determinant(2*pi*Sigma)$modulus),
          gradient = drop(crossprod(Z, part$d1)) - drop(precision %*% b),
          hessian = crossprod(Z, part$d2*Z) - precision)
      }
      b <- rep(0, q)
      for (iteration in 1:200) {
        current <- h(b)
        step <- solve(current$hessian, current$gradient)
        while (h(b - step)$value < current$value - 1e-12 && max(abs(step)) > 1e-14)
          step <- step/2
        b <- b - step
        if (max(abs(step)) < 1e-12) break
      }
      mode <- h(b)
      curvature <- -mode$hessian
      if (is.null(rule))
        return(mode$value + 0.5*q*log(2*pi) - 0.5*as.numeric(determinant(curvature)$modulus))
      scale <- 1/sqrt(curvature[1L, 1L])
      nodes <- rule[, "z"]
      values <- vapply(nodes, function(z) h(b + scale*z)$value, numeric(1)) -
        stats::dnorm(nodes, log = TRUE) + log(rule[, "w"])
      log(scale) + max(values) + log(sum(exp(values - max(values))))
    }, numeric(1))
  }
}

.suest_lme4_components <- function(wrapper, type) {
  parameters <- .suest_lme4_parameters(wrapper, type)
  loglik <- .suest_lme4_group_loglik(wrapper, type)
  groups <- seq_along(levels(.suest_lme4_random_info(wrapper)$id))
  k <- length(parameters)

  # lme4's default Laplace fits stop the inner mode search early (tolPwrss),
  # so their log likelihood can differ slightly from the exact Laplace value.
  glmer <- identical(wrapper$engine, "lme4::glmer")
  nagq <- unname(wrapper$fit@devcomp$dims["nAGQ"])
  laplace <- glmer && nagq == 1L
  # With nAGQ > 1, lme4 reports the log likelihood less that of the saturated
  # model (zero for Bernoulli outcomes).
  y <- as.numeric(lme4::getME(wrapper$fit, "y"))
  saturated <- if (glmer && nagq > 1L && .suest_glmm_family(type) == "poisson")
    sum(stats::dpois(y, y, log = TRUE)) else 0
  total <- sum(loglik(parameters))
  if (!isTRUE(all.equal(total, as.numeric(stats::logLik(wrapper$fit)) + saturated,
      tolerance = if (laplace) 1e-4 else 1e-6)))
    stop("Could not reproduce the lme4 log likelihood for its scores.", call. = FALSE)

  group_scores <- function(theta) {
    step <- 1e-5*pmax(1, abs(theta))
    out <- vapply(seq_len(k), function(j) {
      up <- down <- theta
      up[j] <- up[j] + step[j]
      down[j] <- down[j] - step[j]
      (loglik(up, groups) - loglik(down, groups))/(2*step[j])
    }, numeric(length(groups)))
    matrix(out, nrow = length(groups), dimnames = list(NULL, names(theta)))
  }
  panel_score <- group_scores(parameters)

  step <- 1e-4*pmax(1, abs(parameters))
  information <- -vapply(seq_len(k), function(j) {
    up <- down <- parameters
    up[j] <- up[j] + step[j]
    down[j] <- down[j] - step[j]
    (colSums(group_scores(up)) - colSums(group_scores(down)))/(2*step[j])
  }, numeric(k))
  information <- (information + t(information))/2
  dimnames(information) <- list(names(parameters), names(parameters))

  frame <- wrapper$frame
  id <- .suest_lme4_random_info(wrapper)$id
  n <- nrow(frame)
  score <- matrix(0, n, k, dimnames = list(rownames(frame), names(parameters)))
  for (g in groups) {
    rows <- which(as.integer(id) == g)
    score[rows[length(rows)], ] <- panel_score[g, ]
  }
  B <- n*.suest_information_inverse(information)
  dimnames(B) <- dimnames(information)
  if (any(!is.finite(score)) || any(!is.finite(B)))
    stop("Unable to construct finite lme4 scores.", call. = FALSE)
  list(score = score, bread = B, parameters = parameters)
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
