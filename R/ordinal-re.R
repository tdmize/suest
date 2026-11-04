# Random-intercept ordered logit and probit from ordinal::clmm(), Stata's
# meologit and meoprobit (and xtologit, xtoprobit). P(Y <= j) =
# F(theta_j - x b - u), u ~ N(0, sigma^2). Scores and observed information come
# from suest's own per-group likelihood (Laplace, or adaptive Gauss-Hermite
# with the fit's nAGQ), which reproduces clmm's log likelihood.

.suest_clmm_adapter <- function(model) {
  if (length(model$gfList) != 1L || length(model$ST) != 1L ||
      !identical(dim(model$ST[[1L]]), c(1L, 1L)))
    stop("ordinal::clmm() support is limited to one grouping variable with a random intercept.",
      call. = FALSE)
  if (!model$link %in% c("logit", "probit"))
    stop("Random-effects ordered models must use a logit or probit link.", call. = FALSE)
  if (!identical(model$threshold, "flexible"))
    stop("Random-effects ordered models must use flexible thresholds.", call. = FALSE)
  if (model$nAGQ < 1L)
    stop("Fit the clmm() model with Laplace (nAGQ = 1) or adaptive quadrature (nAGQ > 1).",
      call. = FALSE)
  if (!is.null(model$call$weights) || !is.null(model$call$offset) ||
      !is.null(model$call$nominal) || !is.null(model$call$scale))
    stop("Weights, offsets, and nominal or scale effects are not supported for clmm() models.",
      call. = FALSE)
  if ("log_sigma" %in% names(model$beta))
    stop("clmm() models reserve the coefficient name log_sigma.", call. = FALSE)
  if (model$ST[[1L]][1L, 1L] <= sqrt(.Machine$double.eps))
    stop("The clmm() random-intercept standard deviation must be away from zero.", call. = FALSE)
  list(engine = "ordinal::clmm", type = if (model$link == "logit") "re_ologit" else "re_oprobit")
}

.suest_clmm_parameters <- function(model) {
  if (!is.null(model$suest_parameters)) return(model$suest_parameters)
  c(model$alpha, model$beta, log_sigma = log(model$ST[[1L]][1L, 1L]))
}

.suest_clmm_id <- function(model) droplevels(factor(model$gfList[[1L]]))

.suest_clmm_design <- function(model, data) {
  terms <- stats::delete.response(stats::terms(model))
  frame <- stats::model.frame(terms, data, na.action = stats::na.pass, xlev = model$xlevels)
  X <- stats::model.matrix(terms, frame)
  X[, colnames(X) != "(Intercept)", drop = FALSE]
}

.suest_clmm_cdf <- function(link) {
  if (link == "logit")
    list(F = stats::plogis, f = stats::dlogis,
      fprime = function(t) stats::dlogis(t)*(1 - 2*stats::plogis(t)),
      fsecond = function(t) {
        f <- stats::dlogis(t)
        f*(1 - 2*stats::plogis(t))^2 - 2*f^2
      })
  else
    list(F = stats::pnorm, f = stats::dnorm,
      fprime = function(t) ifelse(is.finite(t), -t*stats::dnorm(t), 0),
      fsecond = function(t) ifelse(is.finite(t), (t^2 - 1)*stats::dnorm(t), 0))
}

# Per-group log likelihood of the suest parameter vector (thresholds, fixed
# effects, log SD) for all groups at once, with exact group scores.
.suest_clmm_group_loglik <- function(model) {
  X <- .suest_clmm_design(model, model$model)
  y <- as.integer(model$model[[1L]])
  id <- .suest_clmm_id(model)
  K <- length(model$alpha)
  p <- ncol(X)
  cdf <- .suest_clmm_cdf(model$link)
  rule <- if (model$nAGQ > 1L) {
    hermite <- .suest_normal_quadrature(model$nAGQ)
    list(z = hermite$nodes, w = hermite$weights)
  }
  engine <- .suest_re_engine(id, X, rule = rule)
  n <- length(y)
  upper <- which(y <= K)
  lower <- which(y >= 2L)
  # Threshold columns: theta_y moves the upper bound, theta_(y-1) the lower one.
  scatter <- function(on_upper, on_lower) {
    out <- matrix(0, n, K)
    out[cbind(upper, y[upper])] <- on_upper[upper]
    out[cbind(lower, y[lower] - 1L)] <- on_lower[lower]
    out
  }

  function(parameters, start = NULL, scores = FALSE) {
    theta <- c(-Inf, parameters[seq_len(K)], Inf)
    eta <- drop(X %*% parameters[K + seq_len(p)])
    density <- function(eta) {
      # P(Y = y) = F(a) - F(c), a = theta_y - eta, c = theta_(y-1) - eta;
      # upper-tail differences keep precision for high categories.
      a <- theta[y + 1L] - eta
      c <- theta[y] - eta
      P <- ifelse(c > 0, cdf$F(-c) - cdf$F(-a), cdf$F(a) - cdf$F(c))
      fa <- cdf$f(a); fc <- cdf$f(c)
      f1a <- cdf$fprime(a); f1c <- cdf$fprime(c)
      g0 <- fa - fc
      g1 <- f1a - f1c
      d1 <- -g0/P
      list(l = log(P), d1 = d1, d2 = g1/P - d1^2,
        d3 = if (scores) -(cdf$fsecond(a) - cdf$fsecond(c))/P + 3*g0*g1/P^2 - 2*g0^3/P^3,
        extra = function(order) {
          if (order == 0L) return(scatter(fa/P, -fc/P))
          if (order == 1L) return(scatter(-f1a/P + g0*fa/P^2, f1c/P - g0*fc/P^2))
          f2a <- cdf$fsecond(a); f2c <- cdf$fsecond(c)
          g1_over <- (f1a - f1c)/P^2
          scatter(f2a/P - g1_over*fa - 2*g0*f1a/P^2 + 2*g0^2*fa/P^3,
            -f2c/P + g1_over*fc + 2*g0*f1c/P^2 - 2*g0^2*fc/P^3)
        })
    }
    engine$intercept(eta, exp(2*parameters[["log_sigma"]]), density, start, scores)
  }
}

.suest_clmm_components <- function(model) {
  parameters <- .suest_clmm_parameters(model)
  evaluate <- .suest_clmm_group_loglik(model)
  # clmm's Laplace fits stop the inner mode search early (gradTol), so their
  # log likelihood can differ slightly from the exact Laplace value.
  base <- evaluate(parameters, scores = TRUE)
  if (!isTRUE(all.equal(sum(base$loglik), as.numeric(stats::logLik(model)),
      tolerance = if (model$nAGQ == 1L) 1e-4 else 1e-6)))
    stop("Could not reproduce the clmm() log likelihood for its scores.", call. = FALSE)
  .suest_re_components(evaluate, parameters, base, .suest_clmm_id(model),
    rownames(model$model))
}

.suest_clmm_frame <- function(model) {
  frame <- as.data.frame(model$model)
  original <- try(eval(model$call$data, envir = environment(stats::formula(model))),
    silent = TRUE)
  if (!inherits(original, "try-error") && is.data.frame(original)) {
    index <- match(rownames(frame), rownames(original))
    if (!anyNA(index))
      for (variable in setdiff(names(original), names(frame)))
        frame[[variable]] <- original[[variable]][index]
  }
  frame
}

# Population-averaged category probabilities.
.suest_clmm_predict <- function(model, newdata) {
  parameters <- .suest_clmm_parameters(model)
  K <- length(model$alpha)
  X <- .suest_clmm_design(model, newdata)
  beta <- parameters[K + seq_len(ncol(X))]
  eta <- drop(X %*% beta)
  sigma <- exp(parameters[["log_sigma"]])
  cumulative <- vapply(seq_len(K), function(j) {
    t <- parameters[j] - eta
    if (model$link == "probit") stats::pnorm(t/sqrt(1 + sigma^2))
    else .suest_logit_normal_mean(t, rep(sigma, length(t)))
  }, numeric(length(eta)))
  cumulative <- matrix(cumulative, nrow = length(eta))
  probabilities <- cbind(cumulative, 1) - cbind(0, cumulative)
  colnames(probabilities) <- model$y.levels
  probabilities
}
