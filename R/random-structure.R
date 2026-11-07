# General random-effects structures (types ending in _re): one top-level
# grouping factor and at most one factor nested in it. Each random-effects term
# is an unstructured (or, for glmmTMB diag(), diagonal) block over its columns:
# the intercept and untransformed numeric slopes. Terms on the same factor are
# independent blocks. Each top-level group's log likelihood is the Laplace
# approximation around the joint mode of its random effects (exact for
# Gaussian models), computed for all groups at once with batched small-matrix
# solves. Parameters are each term's log SDs and the Fisher z (atanh) of its
# correlations, as in Stata's lns and atr.

# terms: list of list(group, columns, covariance, diagonal); factors: named list
# of grouping factors (one value per observation); frame: the estimation frame.
.suest_re_structure <- function(terms, factors, frame) {
  groups <- unique(vapply(terms, `[[`, "", "group"))
  if (length(groups) > 2L)
    stop("Random effects are limited to one grouping variable and at most one grouping variable nested in it.",
      call. = FALSE)
  ids <- lapply(factors[groups], function(f) droplevels(factor(f)))
  if (any(vapply(ids, anyNA, logical(1))))
    stop("The random-effects grouping variables do not align with the estimation sample.", call. = FALSE)
  if (length(groups) == 2L) {
    nested <- function(inner, outer) all(tapply(as.integer(outer), inner, function(v) length(unique(v))) == 1L)
    top <- if (nested(ids[[1L]], ids[[2L]])) 2L else if (nested(ids[[2L]], ids[[1L]])) 1L else
      stop("Random-effects grouping variables must be nested (crossed random effects are not supported).",
        call. = FALSE)
    groups <- groups[c(top, 3L - top)]
    ids <- ids[groups]
  }
  for (term in terms) {
    bad <- setdiff(term$columns, "(Intercept)")
    if (any(!bad %in% names(frame)) || any(make.names(bad) != bad) ||
        any(!vapply(frame[bad], function(x) is.numeric(x) && !is.matrix(x) && all(is.finite(x)), logical(1))))
      stop("Random slopes must be untransformed numeric data columns.", call. = FALSE)
  }
  # Center and spread of each column, for the boundary check.
  terms <- lapply(terms, function(t) {
    t$center <- vapply(t$columns, function(col) if (col == "(Intercept)") 0 else
      mean(frame[[col]]), numeric(1))
    t$spread <- vapply(t$columns, function(col) if (col == "(Intercept)") 1 else
      stats::sd(frame[[col]]), numeric(1))
    t
  })
  level <- vapply(terms, function(t) match(t$group, groups), integer(1))
  terms <- terms[order(level)]
  level <- sort(level)
  label <- function(x) if (x == "(Intercept)") "intercept" else x
  names_of <- function(t) {
    g <- gsub(":", ".", t$group, fixed = TRUE)
    cols <- vapply(t$columns, label, "")
    sds <- paste0("log_sd_", g, "_", cols)
    q <- length(cols)
    pairs <- if (q > 1L && !isTRUE(t$diagonal)) utils::combn(q, 2L) else matrix(0L, 2L, 0L)
    c(sds, if (ncol(pairs)) paste0("atanh_rho_", g, "_", cols[pairs[1L, ]], "_", cols[pairs[2L, ]]))
  }
  values_of <- function(t) {
    S <- as.matrix(t$covariance)
    sd <- sqrt(diag(S))
    q <- length(sd)
    pairs <- if (q > 1L && !isTRUE(t$diagonal)) utils::combn(q, 2L) else matrix(0L, 2L, 0L)
    c(log(sd), atanh(S[t(pairs)]/(sd[pairs[1L, ]]*sd[pairs[2L, ]])))
  }
  counts <- vapply(terms, function(t) length(names_of(t)), integer(1))
  top_id <- ids[[1L]]
  parent <- if (length(groups) == 2L)
    as.integer(top_id)[match(levels(ids[[2L]]), as.character(ids[[2L]]))] else NULL
  list(groups = groups, ids = ids, terms = terms, level = level, parent = parent,
    names = unlist(lapply(terms, names_of)), values = unlist(lapply(terms, values_of)),
    counts = counts, id = top_id, name = groups[1L])
}

.suest_re_structure_lme4 <- function(wrapper) {
  fit <- wrapper$fit
  cnms <- lme4::getME(fit, "cnms")
  covariances <- unclass(lme4::VarCorr(fit))
  terms <- lapply(seq_along(cnms), function(i)
    list(group = names(cnms)[i], columns = unname(cnms[[i]]), covariance = covariances[[i]],
      diagonal = FALSE))
  .suest_re_structure(terms, lme4::getME(fit, "flist"), wrapper$frame)
}

.suest_re_structure_glmmtmb <- function(model) {
  random <- model$modelInfo$reTrms$cond
  blocks <- model$modelInfo$reStruc$condReStruc
  codes <- vapply(blocks, function(b) names(b$blockCode), "")
  if (any(!codes %in% c("us", "diag")))
    stop("glmmTMB random effects must be unstructured (us) or diagonal (diag).", call. = FALSE)
  covariances <- glmmTMB::VarCorr(model)$cond
  terms <- lapply(seq_along(random$cnms), function(i)
    list(group = names(random$cnms)[i], columns = unname(random$cnms[[i]]),
      covariance = covariances[[i]], diagonal = codes[i] == "diag"))
  .suest_re_structure(terms, random$flist, model$frame)
}

# Block-diagonal covariance of each level's random effects from the
# parameters (log SDs, then atanh correlations, term by term).
.suest_re_level_covariances <- function(structure, theta) {
  position <- 0L
  blocks <- lapply(seq_along(structure$terms), function(i) {
    t <- structure$terms[[i]]
    q <- length(t$columns)
    values <- theta[position + seq_len(structure$counts[i])]
    position <<- position + structure$counts[i]
    sd <- exp(values[seq_len(q)])
    R <- diag(q)
    if (q > 1L && !isTRUE(t$diagonal)) {
      pairs <- utils::combn(q, 2L)
      R[t(pairs)] <- tanh(values[-seq_len(q)])
      R[t(pairs[2:1, , drop = FALSE])] <- tanh(values[-seq_len(q)])
    }
    outer(sd, sd)*R
  })
  lapply(seq_along(structure$groups), function(l) {
    parts <- blocks[structure$level == l]
    q <- sum(vapply(parts, nrow, integer(1)))
    out <- matrix(0, q, q)
    at <- 0L
    for (B in parts) {
      index <- at + seq_len(nrow(B))
      out[index, index] <- B
      at <- at + nrow(B)
    }
    out
  })
}

# Each level's random-effect design (one column per random effect) from data.
.suest_re_design <- function(structure, data) {
  lapply(seq_along(structure$groups), function(l) {
    columns <- unlist(lapply(structure$terms[structure$level == l], `[[`, "columns"))
    vapply(columns, function(column)
      if (column == "(Intercept)") rep(1, nrow(data)) else as.numeric(data[[column]]),
      numeric(nrow(data)))
  })
}

# Variance of each observation's random part, for population-averaged predictions.
.suest_re_prediction_variance <- function(structure, theta, newdata) {
  G <- .suest_re_level_covariances(structure, theta)
  Z <- .suest_re_design(structure, newdata)
  Reduce(`+`, Map(function(z, S) rowSums((matrix(z, nrow(newdata)) %*% S)*matrix(z, nrow(newdata))), Z, G))
}

# Batched Cholesky factors of N symmetric q x q matrices (an N x q x q array);
# NULL when one is not positive definite.
.suest_bchol <- function(A) {
  N <- dim(A)[1L]; q <- dim(A)[2L]
  L <- array(0, dim(A))
  for (j in seq_len(q)) {
    prev <- seq_len(j - 1L)
    s <- A[, j, j] - if (j > 1L) rowSums(matrix(L[, j, prev], N)^2) else 0
    if (any(!(s > 0))) return(NULL)
    L[, j, j] <- sqrt(s)
    for (i in seq_len(q)[seq_len(q) > j])
      L[, i, j] <- (A[, i, j] - if (j > 1L) rowSums(matrix(L[, i, prev], N)*matrix(L[, j, prev], N)) else 0)/L[, j, j]
  }
  L
}

# Solves (L L') x = b for each of the N systems; b is N x q x m.
.suest_bsolve <- function(L, b) {
  N <- dim(L)[1L]; q <- dim(L)[2L]; m <- dim(b)[3L]
  y <- x <- array(0, c(N, q, m))
  for (i in seq_len(q)) {
    acc <- matrix(b[, i, ], N)
    for (k in seq_len(i - 1L)) acc <- acc - L[, i, k]*matrix(y[, k, ], N)
    y[, i, ] <- acc/L[, i, i]
  }
  for (i in rev(seq_len(q))) {
    acc <- matrix(y[, i, ], N)
    for (k in seq_len(q)[seq_len(q) > i]) acc <- acc - L[, k, i]*matrix(x[, k, ], N)
    x[, i, ] <- acc/L[, i, i]
  }
  x
}

.suest_blogdet <- function(L) 2*rowSums(log(matrix(vapply(seq_len(dim(L)[2L]),
  function(j) L[, j, j], numeric(dim(L)[1L])), dim(L)[1L])))

# Per-group sums of w * a_j * b_k for all column pairs: an N x qa x qb array.
.suest_bouter <- function(A, B, w, id, N) {
  out <- array(0, c(N, ncol(A), ncol(B)))
  for (j in seq_len(ncol(A))) for (k in seq_len(ncol(B)))
    out[, j, k] <- rowsum(w*A[, j]*B[, k], id, reorder = TRUE)[, 1L]
  out
}

# Batched matrix products: N x a x b times N x b x c.
.suest_bprod <- function(A, B) {
  out <- array(0, c(dim(A)[1L], dim(A)[2L], dim(B)[3L]))
  for (i in seq_len(dim(A)[2L])) for (k in seq_len(dim(B)[3L]))
    out[, i, k] <- rowSums(matrix(A[, i, ], dim(A)[1L])*matrix(B[, , k], dim(A)[1L]))
  out
}

# Laplace engine for a structure and fixed-effect design. The returned function
# takes the fixed-part linear predictor, the level covariances, the density
# (log density l and derivatives d1, d2 in eta), and optional starting modes;
# it returns each top-level group's log likelihood and the modes.
.suest_re_laplace <- function(structure, data) {
  Z <- .suest_re_design(structure, data)
  ids <- lapply(structure$ids, as.integer)
  N <- vapply(structure$ids, nlevels, integer(1))
  nested <- length(Z) == 2L
  parent <- structure$parent
  top_sum <- function(x) rowsum(x, ids[[1L]], reorder = TRUE)[, 1L]
  up <- function(x) rowsum(x, parent, reorder = TRUE)[, 1L]
  random_part <- function(modes)
    rowSums(Z[[1L]]*modes$u[ids[[1L]], , drop = FALSE]) +
      if (nested) rowSums(Z[[2L]]*modes$v[ids[[2L]], , drop = FALSE]) else 0

  function(eta, G, density, start = NULL) {
    P <- lapply(G, solve)
    logdet_G <- vapply(G, function(S) as.numeric(determinant(S)$modulus), numeric(1))
    q <- vapply(G, nrow, integer(1))
    modes <- if (is.null(start)) list(u = matrix(0, N[1L], q[1L]),
      v = if (nested) matrix(0, N[2L], q[2L])) else start
    objective <- function(modes, part) {
      out <- top_sum(part$l) - 0.5*rowSums((modes$u %*% P[[1L]])*modes$u) - 0.5*logdet_G[1L]
      if (nested) out <- out + up(-0.5*rowSums((modes$v %*% P[[2L]])*modes$v) - 0.5*logdet_G[2L])
      out
    }
    broadcast <- function(M, n) aperm(array(M, c(dim(M), n)), c(3L, 1L, 2L))
    step_of <- function(modes, part) {
      w <- -part$d2
      A <- broadcast(P[[1L]], N[1L]) + .suest_bouter(Z[[1L]], Z[[1L]], w, ids[[1L]], N[1L])
      gu <- rowsum(part$d1*Z[[1L]], ids[[1L]], reorder = TRUE) - modes$u %*% P[[1L]]
      if (!nested) {
        LS <- .suest_bchol(A)
        if (is.null(LS)) stop("The random-effects curvature is not positive definite.", call. = FALSE)
        return(list(du = matrix(.suest_bsolve(LS, array(gu, c(dim(gu), 1L))), N[1L]),
          logdet = .suest_blogdet(LS)))
      }
      B <- broadcast(P[[2L]], N[2L]) + .suest_bouter(Z[[2L]], Z[[2L]], w, ids[[2L]], N[2L])
      C <- .suest_bouter(Z[[1L]], Z[[2L]], w, ids[[2L]], N[2L])
      gv <- rowsum(part$d1*Z[[2L]], ids[[2L]], reorder = TRUE) - modes$v %*% P[[2L]]
      LB <- .suest_bchol(B)
      if (is.null(LB)) stop("The random-effects curvature is not positive definite.", call. = FALSE)
      X <- .suest_bsolve(LB, aperm(C, c(1L, 3L, 2L)))
      yv <- matrix(.suest_bsolve(LB, array(gv, c(dim(gv), 1L))), N[2L])
      CX <- .suest_bprod(C, X)
      S <- A - array(apply(CX, 2:3, up), dim(A))
      r <- gu - apply(matrix(.suest_bprod(C, array(yv, c(dim(yv), 1L))), N[2L]), 2L, up)
      LS <- .suest_bchol(S)
      if (is.null(LS)) stop("The random-effects curvature is not positive definite.", call. = FALSE)
      du <- matrix(.suest_bsolve(LS, array(r, c(N[1L], q[1L], 1L))), N[1L])
      dv <- yv - matrix(.suest_bprod(X, array(du[parent, , drop = FALSE], c(N[2L], q[1L], 1L))), N[2L])
      list(du = du, dv = dv, logdet = .suest_blogdet(LS) + up(.suest_blogdet(LB)))
    }
    part <- density(eta + random_part(modes))
    value <- objective(modes, part)
    for (iteration in 1:200) {
      step <- step_of(modes, part)
      size <- max(abs(c(step$du, step$dv)))
      if (!is.finite(size)) stop("Could not find the random-effect modes.", call. = FALSE)
      scale <- rep(1, N[1L])
      repeat {
        trial_modes <- list(u = modes$u + scale*step$du,
          v = if (nested) modes$v + scale[parent]*step$dv)
        trial <- density(eta + random_part(trial_modes))
        trial_value <- objective(trial_modes, trial)
        worse <- !(trial_value >= value - 1e-12) & scale > 1e-10
        if (!any(worse)) break
        scale[worse] <- scale[worse]/2
      }
      modes <- trial_modes
      part <- trial
      value <- trial_value
      if (size < 1e-11) break
    }
    list(loglik = value - 0.5*step_of(modes, part)$logdet, mode = modes)
  }
}

# Group log likelihood of a general structure in suest's parameters: the fixed
# effects, the dispersion (when offset > ncol(X)), then the random effects.
.suest_re_general_loglik <- function(structure, frame, X, y, family, offset) {
  laplace <- .suest_re_laplace(structure, frame)
  p <- ncol(X)
  re <- offset + seq_along(structure$names)
  function(parameters, start = NULL, scores = FALSE) {
    eta <- drop(X %*% parameters[seq_len(p)])
    G <- .suest_re_level_covariances(structure, parameters[re])
    density <- if (family == "gaussian") {
      s <- exp(parameters[[p + 1L]])
      function(e) list(l = stats::dnorm(y, e, s, log = TRUE), d1 = (y - e)/s^2,
        d2 = rep(-1/s^2, length(y)))
    } else {
      dispersion <- if (offset > p) exp(parameters[[p + 1L]])
      function(e) .suest_glmm_density(family, y, e, dispersion)
    }
    laplace(eta, G, density, start)
  }
}

.suest_glmmtmb_re_components <- function(model, type) {
  parameters <- .suest_extract_parameters(model, type, "glmmTMB::glmmTMB")
  data <- model$obj$env$data
  X <- as.matrix(data$X)
  structure <- .suest_glmmtmb_random_info(model)$structure
  offset <- ncol(X) + !is.na(.suest_glmm_dispersion_name(type))
  evaluate <- .suest_re_general_loglik(structure, model$frame, X, as.numeric(data$yobs),
    .suest_glmm_family(type), offset)
  base <- evaluate(parameters)
  if (!isTRUE(all.equal(sum(base$loglik), as.numeric(stats::logLik(model)), tolerance = 1e-6)))
    stop("Could not reproduce the glmmTMB log likelihood for its scores.", call. = FALSE)
  .suest_re_components(evaluate, parameters, base, structure$id,
    rownames(.suest_model_frame(model, "glmmTMB::glmmTMB")))
}

# Random-effect covariances must be away from zero variances and singular
# correlations (in a standardized basis).
.suest_re_check_covariances <- function(structure) {
  for (t in structure$terms) {
    S <- as.matrix(t$covariance)
    sd <- sqrt(diag(S))
    bad <- any(!is.finite(S)) || any(!is.finite(t$spread)) || any(sd <= sqrt(.Machine$double.eps))
    if (!bad) {
      # The covariance of the term's part of the linear predictor in a
      # centered, standardized basis, so predictor units do not matter.
      q <- length(sd)
      T <- diag(t$spread, q)
      i <- match("(Intercept)", t$columns)
      if (!is.na(i)) T[i, ] <- ifelse(seq_len(q) == i, 1, t$center)
      E <- eigen(T %*% S %*% t(T), symmetric = TRUE, only.values = TRUE)$values
      R <- eigen(S/outer(sd, sd), symmetric = TRUE, only.values = TRUE)$values
      bad <- min(E) <= sqrt(.Machine$double.eps)*max(1, E) || min(R) <= 1e-4
    }
    if (bad)
      stop(paste0("The random-effect covariance is nearly singular (a variance near zero or ",
        "correlations near -1 or 1); suest needs a fit away from that boundary. ",
        "Simplify the random-effects structure."), call. = FALSE)
  }
}
