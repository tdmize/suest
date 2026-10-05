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
      fprime = function(t) stats::dlogis(t)*(1 - 2*stats::plogis(t)))
  else
    list(F = stats::pnorm, f = stats::dnorm,
      fprime = function(t) ifelse(is.finite(t), -t*stats::dnorm(t), 0))
}

.suest_clmm_group_loglik <- function(model) {
  X <- .suest_clmm_design(model, model$model)
  y <- as.integer(model$model[[1L]])
  id <- .suest_clmm_id(model)
  rows <- split(seq_along(y), id)
  K <- length(model$alpha)
  p <- ncol(X)
  cdf <- .suest_clmm_cdf(model$link)
  rule <- if (model$nAGQ > 1L) .suest_normal_quadrature(model$nAGQ) else NULL

  function(parameters, groups = seq_along(rows)) {
    theta <- c(-Inf, parameters[seq_len(K)], Inf)
    beta <- parameters[K + seq_len(p)]
    sigma2 <- exp(2*parameters[["log_sigma"]])
    eta_all <- drop(X %*% beta)
    upper_all <- theta[y + 1L]
    lower_all <- theta[y]
    vapply(groups, function(g) {
      i <- rows[[g]]
      a <- upper_all[i] - eta_all[i]
      c <- lower_all[i] - eta_all[i]
      h <- function(b) {
        # Upper-tail differences keep precision for high categories.
        P <- ifelse(c - b > 0, cdf$F(b - c) - cdf$F(b - a), cdf$F(a - b) - cdf$F(c - b))
        fa <- cdf$f(a - b); fc <- cdf$f(c - b)
        d1 <- -(fa - fc)/P
        d2 <- (cdf$fprime(a - b) - cdf$fprime(c - b))/P - d1^2
        list(value = sum(log(P)) - b^2/(2*sigma2) - 0.5*log(2*pi*sigma2),
          gradient = sum(d1) - b/sigma2, hessian = sum(d2) - 1/sigma2)
      }
      b <- 0
      for (iteration in 1:200) {
        current <- h(b)
        step <- current$gradient/current$hessian
        if (!is.finite(step))
          stop("Could not find the clmm() conditional modes.", call. = FALSE)
        while (abs(step) > 1e-14 && !isTRUE(h(b - step)$value >= current$value - 1e-12))
          step <- step/2
        b <- b - step
        if (abs(step) < 1e-12) break
      }
      mode <- h(b)
      curvature <- -mode$hessian
      if (is.null(rule))
        return(mode$value + 0.5*log(2*pi) - 0.5*log(curvature))
      scale <- 1/sqrt(curvature)
      nodes <- rule$nodes
      values <- vapply(nodes, function(z) h(b + scale*z)$value, numeric(1)) -
        stats::dnorm(nodes, log = TRUE) + log(rule$weights)
      log(scale) + max(values) + log(sum(exp(values - max(values))))
    }, numeric(1))
  }
}

.suest_clmm_components <- function(model) {
  parameters <- .suest_clmm_parameters(model)
  loglik <- .suest_clmm_group_loglik(model)
  id <- .suest_clmm_id(model)
  groups <- seq_along(levels(id))
  k <- length(parameters)
  # clmm's Laplace fits stop the inner mode search early (gradTol), so their
  # log likelihood can differ slightly from the exact Laplace value.
  total <- sum(loglik(parameters))
  if (!isTRUE(all.equal(total, as.numeric(stats::logLik(model)),
      tolerance = if (model$nAGQ == 1L) 1e-4 else 1e-6)))
    stop("Could not reproduce the clmm() log likelihood for its scores.", call. = FALSE)

  group_scores <- function(theta) {
    step <- 1e-5*pmax(1, abs(theta))
    matrix(vapply(seq_len(k), function(j) {
      up <- down <- theta
      up[j] <- up[j] + step[j]
      down[j] <- down[j] - step[j]
      (loglik(up, groups) - loglik(down, groups))/(2*step[j])
    }, numeric(length(groups))), nrow = length(groups))
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

  frame <- model$model
  n <- nrow(frame)
  score <- matrix(0, n, k, dimnames = list(rownames(frame), names(parameters)))
  for (g in groups) {
    rows <- which(as.integer(id) == g)
    score[rows[length(rows)], ] <- panel_score[g, ]
  }
  B <- n*.suest_information_inverse(information)
  dimnames(B) <- dimnames(information)
  list(score = score, bread = B, parameters = parameters)
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
