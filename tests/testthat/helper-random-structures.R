# Data and an independent likelihood for random-effects structures with
# several slopes and a nested level: 40 schools, 3 classes in each, 6 students
# in each class.
random_structure_data <- function(schools = 40L, classes = 3L, size = 6L, seed = 3104) {
  set.seed(seed)
  n <- schools*classes*size
  d <- data.frame(school = factor(rep(seq_len(schools), each = classes*size)),
    class = factor(rep(seq_len(classes), each = size, times = schools)),
    x1 = rnorm(n), x2 = rnorm(n))
  R <- matrix(c(1, .3, -.2, .3, 1, .1, -.2, .1, 1), 3)
  sd <- c(.8, .5, .45)
  U <- matrix(rnorm(schools*3), schools) %*% chol(outer(sd, sd)*R)
  s <- as.integer(d$school)
  k <- (s - 1L)*classes + as.integer(d$class)
  v <- rnorm(schools*classes, sd = .6)[k]
  re <- U[s, 1] + U[s, 2]*d$x1 + U[s, 3]*d$x2
  d$yc <- 1 + .5*d$x1 - .3*d$x2 + re + v + rnorm(n)
  d$yb <- rbinom(n, 1L, plogis(-.2 + .6*d$x1 - .4*d$x2 + re + v))
  d$count <- rpois(n, exp(.2 + .4*d$x1 - .2*d$x2 + .6*re + .6*v))
  d$yo <- ordered(cut(.6*d$x1 - .3*d$x2 + re + v + rlogis(n), c(-Inf, -1, .5, 2, Inf),
    labels = c("a", "b", "c", "d")))
  d$nb <- rnbinom(n, size = 2, mu = exp(.3 + .3*d$x1 + .6*re + .6*v))
  d$positive <- rgamma(n, shape = 3, scale = exp(.2 + .3*d$x1 + .5*re + .5*v)/3)
  d$id <- seq_len(n)
  d
}

# glmmTMB's covariance parameters are log SDs and, for unstructured terms, the
# entries of a scaled Cholesky factor; suest uses log SDs and Fisher z
# correlations. This maps glmmTMB's to suest's, named as suest names them.
rs_glmmtmb_covariance <- function(model, theta) {
  info <- model$modelInfo$reTrms$cond
  codes <- names(glmmTMB:::.valid_covstruct)[match(vapply(model$modelInfo$reStruc$condReStruc,
    `[[`, 0, "blockCode"), glmmTMB:::.valid_covstruct)]
  at <- 0L
  unlist(unname(Map(function(cols, code, group) {
    q <- length(cols); k <- if (code == "diag") q else q*(q + 1L)/2L
    part <- theta[at + seq_len(k)]; at <<- at + k
    label <- ifelse(cols == "(Intercept)", "intercept", cols); g <- gsub(":", ".", group)
    out <- stats::setNames(part[seq_len(q)], paste0("log_sd_", g, "_", label))
    if (code == "us" && q > 1L) {
      pairs <- utils::combn(q, 2L)
      out <- c(out, stats::setNames(atanh(glmmTMB:::get_cor(part[-seq_len(q)])),
        paste0("atanh_rho_", g, "_", label[pairs[1, ]], "_", label[pairs[2, ]])))
    }
    out
  }, info$cnms, codes, names(info$cnms))))
}

# The Jacobian of suest's parameters with respect to glmmTMB's (fixed effects,
# dispersion, covariance), with suest's names on the rows.
rs_glmmtmb_jacobian <- function(model) {
  par <- model$fit$par
  theta <- par[names(par) == "theta"]
  k <- sum(names(par) != "theta")
  J <- diag(length(par))
  J[-seq_len(k), -seq_len(k)] <- numDeriv::jacobian(function(t)
    rs_glmmtmb_covariance(model, t), theta)
  rownames(J) <- c(names(par)[seq_len(k)], names(rs_glmmtmb_covariance(model, theta)))
  J
}

# One random-effects term: the name used in suest's parameter names, the
# grouping factor, and the columns ("intercept" or a data column).
rs_term <- function(name, factor, columns, diagonal = FALSE)
  list(name = name, factor = droplevels(factor(factor)), columns = columns, diagonal = diagonal)

rs_covariance <- function(theta, term) {
  q <- length(term$columns)
  sd <- exp(theta[paste0("log_sd_", term$name, "_", term$columns)])
  R <- diag(q)
  if (q > 1L && !term$diagonal) for (i in 1:(q - 1L)) for (j in (i + 1L):q)
    R[i, j] <- R[j, i] <- tanh(theta[[paste0("atanh_rho_", term$name, "_",
      term$columns[i], "_", term$columns[j])]])
  outer(sd, sd)*R
}

# Conditional log density of each observation and its first two derivatives
# with respect to the linear predictor.
rs_density <- function(family, y, e, theta, thresholds) {
  switch(family,
    logit = {
      p <- plogis(e)
      list(l = y*e - pmax(e, 0) - log1p(exp(-abs(e))), d1 = y - p, d2 = -p*(1 - p))
    },
    poisson = {
      m <- exp(e)
      list(l = y*e - m - lgamma(y + 1), d1 = y - m, d2 = -m)
    },
    ologit = {
      cut <- c(-Inf, theta[thresholds], Inf)
      a <- cut[y] - e; b <- cut[y + 1L] - e
      Fa <- plogis(a); Fb <- plogis(b); fa <- dlogis(a); fb <- dlogis(b)
      ga <- ifelse(is.finite(a), fa*(1 - 2*Fa), 0); gb <- ifelse(is.finite(b), fb*(1 - 2*Fb), 0)
      P <- Fb - Fa
      list(l = log(P), d1 = -(fb - fa)/P, d2 = ((gb - ga)*P - (fb - fa)^2)/P^2)
    })
}

# Each school's log likelihood: exact for Gaussian models, otherwise the
# Laplace approximation around the joint mode of all of the school's random
# effects (both levels), found by Newton's method on dense matrices.
rs_school_loglik <- function(theta, spec, schools = levels(spec$school)) {
  X <- spec$X
  eta <- drop(X %*% theta[colnames(X)])
  vapply(schools, function(s) {
    rows <- spec$school == s
    blocks <- lapply(spec$terms, function(term) {
      f <- droplevels(term$factor[rows])
      Zt <- sapply(term$columns, function(col)
        if (col == "intercept") rep(1, sum(rows)) else spec$data[[col]][rows])
      Zt <- matrix(Zt, sum(rows))
      Z <- do.call(cbind, lapply(levels(f), function(l) Zt*(f == l)))
      list(Z = Z, G = kronecker(diag(nlevels(f)), rs_covariance(theta, term)))
    })
    Z <- do.call(cbind, lapply(blocks, `[[`, "Z"))
    sizes <- vapply(blocks, function(b) ncol(b$G), integer(1)); G <- diag(0, sum(sizes))
    at <- cumsum(c(0L, sizes))
    for (b in seq_along(blocks)) G[at[b] + seq_len(sizes[b]), at[b] + seq_len(sizes[b])] <- blocks[[b]]$G
    y <- spec$y[rows]; e0 <- eta[rows]
    if (spec$family == "gaussian") {
      V <- Z %*% G %*% t(Z) + diag(exp(2*theta[["log_sigma_e"]]), sum(rows))
      r <- y - e0
      return(-(length(r)*log(2*pi) + as.numeric(determinant(V, logarithm = TRUE)$modulus) +
        sum(r*solve(V, r)))/2)
    }
    P <- solve(G)
    objective <- function(u) {
      k <- rs_density(spec$family, y, e0 + drop(Z %*% u), theta, spec$thresholds)
      list(value = sum(k$l) - sum(u*(P %*% u))/2,
        gradient = drop(crossprod(Z, k$d1)) - drop(P %*% u),
        hessian = crossprod(Z, Z*k$d2) - P)
    }
    u <- rep(0, ncol(Z)); current <- objective(u)
    for (iteration in 1:100) {
      step <- -solve(current$hessian, current$gradient)
      t <- 1
      repeat {
        candidate <- objective(u + t*step)
        if (candidate$value >= current$value - 1e-12 || t < 1e-10) break
        t <- t/2
      }
      u <- u + t*step; current <- candidate
      if (max(abs(t*step)) < 1e-12) break
    }
    stopifnot(max(abs(current$gradient)) < 1e-8)
    current$value - as.numeric(determinant(G, logarithm = TRUE)$modulus)/2 -
      as.numeric(determinant(-current$hessian, logarithm = TRUE)$modulus)/2
  }, numeric(1))
}

rs_numeric_scores <- function(theta, spec, schools, h = 1e-5) {
  out <- vapply(seq_along(theta), function(j) {
    up <- down <- theta
    up[j] <- up[j] + h; down[j] <- down[j] - h
    (rs_school_loglik(up, spec, schools) - rs_school_loglik(down, spec, schools))/(2*h)
  }, numeric(length(schools)))
  matrix(out, length(schools), dimnames = list(schools, names(theta)))
}

rs_spec <- function(data, formula, y, family, terms, thresholds = NULL) {
  X <- stats::model.matrix(formula, data)
  list(data = data, X = X, y = y, family = family, terms = terms,
    thresholds = thresholds, school = droplevels(data$school))
}
