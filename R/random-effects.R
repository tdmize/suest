# Group log likelihoods of random-intercept (or intercept and slope) models,
# computed for all groups at once. Gaussian models are exact. Otherwise each
# group's integral over its random effects uses Laplace or adaptive
# Gauss-Hermite quadrature around the group's mode. density(eta) returns each
# observation's log density l and its derivatives d1, d2, d3 in eta, and
# optionally extra(order): the derivatives of l, d1 (order 1), or d2 (order 2)
# in parameters other than the fixed effects, one column per parameter.

.suest_re_engine <- function(id, X, slope = NULL, rule = NULL) {
  gid <- as.integer(id)
  G <- nlevels(id)
  gsum <- function(x) {
    out <- unname(rowsum(x, gid, reorder = TRUE))
    if (is.matrix(x)) out else drop(out)
  }
  size <- gsum(rep(1, length(gid)))
  if (!is.null(slope)) {
    slope_sum <- gsum(slope)
    slope_square <- gsum(slope^2)
  }
  # Per-observation derivatives in the extra parameters, then the fixed effects.
  columns <- function(part, order) {
    d <- switch(order + 1L, part$d1, part$d2, part$d3)
    cbind(if (!is.null(part$extra)) part$extra(order), X*d)
  }

  gaussian <- function(y, eta, Sigma, residual) {
    r <- y - eta
    total <- gsum(r)
    square <- gsum(r^2)
    if (is.null(slope)) {
      d <- residual + size*Sigma[1L, 1L]
      return(-0.5*(size*log(2*pi) + (size - 1)*log(residual) + log(d) +
        (square - Sigma[1L, 1L]*total^2/d)/residual))
    }
    # Woodbury identity with M = residual * Sigma^-1 + Z'Z for each group.
    P <- solve(Sigma)
    cross <- gsum(slope*r)
    m11 <- residual*P[1L, 1L] + size
    m12 <- residual*P[1L, 2L] + slope_sum
    m22 <- residual*P[2L, 2L] + slope_square
    det_m <- m11*m22 - m12^2
    quad <- (square - (m22*total^2 - 2*m12*total*cross + m11*cross^2)/det_m)/residual
    -0.5*(size*log(2*pi) + (size - 2)*log(residual) +
      as.numeric(determinant(Sigma)$modulus) + log(det_m) + quad)
  }

  # With scores = TRUE, also the exact group scores in (extra parameters,
  # fixed effects, log SD): the approximation depends on the parameters
  # directly and through each group's mode and curvature.
  intercept <- function(eta, s2, density, start = NULL, scores = FALSE) {
    const <- -0.5*log(2*pi*s2)
    value_of <- function(b, part) gsum(part$l) - 0.5*b^2/s2 + const
    b <- if (is.null(start)) numeric(G) else start
    part <- density(eta + b[gid])
    value <- value_of(b, part)
    for (iteration in 1:200) {
      step <- (gsum(part$d1) - b/s2)/(gsum(part$d2) - 1/s2)
      if (any(!is.finite(step)))
        stop("Could not find the random-effect modes.", call. = FALSE)
      repeat {
        trial <- density(eta + (b - step)[gid])
        trial_value <- value_of(b - step, trial)
        worse <- !(trial_value >= value - 1e-12) & abs(step) > 1e-14
        if (!any(worse)) break
        step[worse] <- step[worse]/2
      }
      b <- b - step
      part <- trial
      value <- trial_value
      if (max(abs(step)) < 1e-12) break
    }
    curvature <- 1/s2 - gsum(part$d2)
    if (scores) {
      d_mode <- cbind(gsum(columns(part, 1L)), 2*b/s2)/curvature
      d_curvature <- -(cbind(gsum(columns(part, 2L)), 2/s2) + gsum(part$d3)*d_mode)
    }
    if (is.null(rule)) {
      loglik <- value + 0.5*log(2*pi) - 0.5*log(curvature)
      score <- if (scores)
        cbind(gsum(columns(part, 0L)), b^2/s2 - 1) - 0.5*d_curvature/curvature
      return(list(loglik = loglik, mode = b, score = score))
    }
    # Quadrature weights are accumulated with a running maximum per group.
    scale <- 1/sqrt(curvature)
    top <- rep(-Inf, G)
    total <- numeric(G)
    if (scores) {
      direct <- 0
      slope_mean <- slope_z <- square <- numeric(G)
    }
    for (k in seq_along(rule$z)) {
      a <- b + scale*rule$z[k]
      node <- density(eta + a[gid])
      v <- value_of(a, node) - stats::dnorm(rule$z[k], log = TRUE) + log(rule$w[k])
      new_top <- pmax(top, v)
      shrink <- exp(top - new_top)
      w <- exp(v - new_top)
      total <- total*shrink + w
      if (scores) {
        direct <- direct*shrink[gid] + columns(node, 0L)*w[gid]
        h_b <- gsum(node$d1) - a/s2
        slope_mean <- slope_mean*shrink + w*h_b
        slope_z <- slope_z*shrink + w*h_b*rule$z[k]
        square <- square*shrink + w*a^2
      }
      top <- new_top
    }
    loglik <- log(scale) + top + log(total)
    if (!scores) return(list(loglik = loglik, mode = b))
    d_scale <- -0.5*scale*d_curvature/curvature
    score <- d_scale/scale + cbind(gsum(direct)/total, square/total/s2 - 1) +
      slope_mean/total*d_mode + slope_z/total*d_scale
    list(loglik = loglik, mode = b, score = score)
  }

  # Random intercept and slope: Laplace around a two-dimensional mode.
  intercept_slope <- function(eta, Sigma, density, start = NULL) {
    P <- solve(Sigma)
    const <- -0.5*as.numeric(determinant(2*pi*Sigma)$modulus)
    b <- if (is.null(start)) matrix(0, G, 2L) else start
    state <- function(b) {
      part <- density(eta + b[gid, 1L] + b[gid, 2L]*slope)
      pb1 <- P[1L, 1L]*b[, 1L] + P[1L, 2L]*b[, 2L]
      pb2 <- P[1L, 2L]*b[, 1L] + P[2L, 2L]*b[, 2L]
      list(value = gsum(part$l) - 0.5*(b[, 1L]*pb1 + b[, 2L]*pb2) + const,
        g1 = gsum(part$d1) - pb1, g2 = gsum(part$d1*slope) - pb2,
        h11 = gsum(part$d2) - P[1L, 1L], h12 = gsum(part$d2*slope) - P[1L, 2L],
        h22 = gsum(part$d2*slope^2) - P[2L, 2L])
    }
    current <- state(b)
    for (iteration in 1:200) {
      det_h <- current$h11*current$h22 - current$h12^2
      step <- cbind(current$h22*current$g1 - current$h12*current$g2,
        current$h11*current$g2 - current$h12*current$g1)/det_h
      if (any(!is.finite(step)))
        stop("Could not find the random-effect modes.", call. = FALSE)
      repeat {
        trial <- state(b - step)
        worse <- !(trial$value >= current$value - 1e-12) &
          pmax(abs(step[, 1L]), abs(step[, 2L])) > 1e-14
        if (!any(worse)) break
        step[worse, ] <- step[worse, ]/2
      }
      b <- b - step
      current <- trial
      if (max(abs(step)) < 1e-12) break
    }
    det_curvature <- current$h11*current$h22 - current$h12^2
    list(loglik = current$value + log(2*pi) - 0.5*log(det_curvature), mode = b)
  }

  list(gaussian = gaussian, intercept = intercept, intercept_slope = intercept_slope)
}

# Central-difference group scores of a group log-likelihood function, starting
# each evaluation from the modes in base (the evaluation at the estimates).
.suest_re_numeric_scores <- function(evaluate, theta, base) {
  step <- 1e-5*pmax(1, abs(theta))
  out <- vapply(seq_along(theta), function(j) {
    up <- down <- theta
    up[j] <- up[j] + step[j]
    down[j] <- down[j] - step[j]
    (evaluate(up, base$mode)$loglik - evaluate(down, base$mode)$loglik)/(2*step[j])
  }, numeric(length(base$loglik)))
  matrix(out, ncol = length(theta))
}

# Observed information from central differences of the summed group scores,
# and each group's score placed on the group's last observation.
.suest_re_components <- function(evaluate, parameters, base, id, rows) {
  k <- length(parameters)
  group_scores <- function(theta) {
    if (!is.null(base$score)) return(evaluate(theta, base$mode, scores = TRUE)$score)
    .suest_re_numeric_scores(evaluate, theta, base)
  }
  panel_score <- if (is.null(base$score)) group_scores(parameters) else base$score
  step <- 1e-4*pmax(1, abs(parameters))
  information <- -vapply(seq_len(k), function(j) {
    up <- down <- parameters
    up[j] <- up[j] + step[j]
    down[j] <- down[j] - step[j]
    (colSums(group_scores(up)) - colSums(group_scores(down)))/(2*step[j])
  }, numeric(k))
  information <- (information + t(information))/2
  dimnames(information) <- list(names(parameters), names(parameters))

  id <- as.integer(id)
  n <- length(id)
  score <- matrix(0, n, k, dimnames = list(rows, names(parameters)))
  last <- which(!duplicated(id, fromLast = TRUE))
  score[last, ] <- panel_score[id[last], ]
  B <- n*.suest_information_inverse(information)
  dimnames(B) <- dimnames(information)
  if (any(!is.finite(score)) || any(!is.finite(B)) || any(diag(B) <= 0))
    stop("Unable to construct finite random-effects scores.", call. = FALSE)
  list(score = score, bread = B, parameters = parameters)
}

# Group scores of a glmmTMB model in its own parameters (fixed effects,
# dispersion, then log SDs and the scaled correlation), by central differences
# of its Laplace group log likelihoods, all groups at once.
.suest_glmmtmb_group_scores <- function(model, type) {
  data <- model$obj$env$data
  X <- as.matrix(data$X)
  y <- as.numeric(data$yobs)
  random <- .suest_glmmtmb_random_info(model)
  slope <- if (is.null(random$slope)) NULL else model$frame[[random$slope]]
  family <- .suest_glmm_family(type)
  par <- model$fit$par
  beta <- which(names(par) == "beta")
  disp <- which(names(par) == "betadisp")
  theta <- which(names(par) == "theta")
  if (length(beta) != ncol(X) || nrow(X) != length(random$id) ||
      length(disp) != !is.na(.suest_glmm_dispersion_name(type)))
    stop("The glmmTMB parameters do not align with its data.", call. = FALSE)
  engine <- .suest_re_engine(random$id, X, slope)

  evaluate <- function(parameters, start = NULL) {
    eta <- drop(X %*% parameters[beta])
    re <- parameters[theta]
    Sigma <- if (is.null(slope)) matrix(exp(2*re)) else {
      sd <- exp(re[1:2])
      rho <- re[3]/sqrt(1 + re[3]^2)
      matrix(c(sd[1]^2, rho*sd[1]*sd[2], rho*sd[1]*sd[2], sd[2]^2), 2L)
    }
    if (family == "gaussian")
      return(list(loglik = engine$gaussian(y, eta, Sigma, exp(2*parameters[disp]))))
    dispersion <- if (length(disp)) exp(parameters[disp])
    density <- function(eta) .suest_glmm_density(family, y, eta, dispersion)
    if (is.null(slope)) engine$intercept(eta, Sigma[1L, 1L], density, start)
    else engine$intercept_slope(eta, Sigma, density, start)
  }
  base <- evaluate(par)
  if (!isTRUE(all.equal(sum(base$loglik), as.numeric(stats::logLik(model)), tolerance = 1e-6)))
    stop("Could not reproduce the glmmTMB log likelihood for its scores.", call. = FALSE)
  score <- .suest_re_numeric_scores(evaluate, unname(par), base)
  rownames(score) <- levels(random$id)
  score
}
