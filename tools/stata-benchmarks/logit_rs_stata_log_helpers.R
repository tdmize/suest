# Parse printed full-precision Stata matrices, with strict shape checks.
logitrs_read_matrix <- function(block, name, panels = FALSE) {
  headers <- grep("^(symmetric )?[er]\\([[:alnum:]_]+\\)\\[", block)
  start <- headers[grepl(paste0(name, "["), block[headers], fixed = TRUE)]
  stopifnot(length(start) == 1L)
  next_header <- headers[headers > start]
  end <- if (length(next_header)) next_header[1L]-1L else length(block)
  rows <- list(); equation <- ""
  for (line in block[(start+1L):end]) {
    token <- strsplit(trimws(line), "[[:space:]]+")[[1L]]
    if (length(token) == 1L && grepl(":$", token)) equation <- token
    if (length(token) < 2L) next
    values <- suppressWarnings(as.numeric(token[-1L]))
    if (anyNA(values)) next
    if (panels) {
      key <- if (grepl(":", token[1L])) token[1L] else paste0(equation, token[1L])
      rows[[key]] <- c(rows[[key]], values)
    } else rows[[length(rows)+1L]] <- values
  }
  shape <- as.integer(strsplit(sub(".*\\[([0-9,]+)\\].*", "\\1", block[start]), ",")[[1L]])
  if (name %in% c("e(b)", "r(b)", "e(gradient)")) {
    stopifnot(shape[1L] == 1L, length(rows) == 1L, length(rows[[1L]]) == shape[2L])
    return(rows[[1L]])
  }
  p <- length(rows)
  stopifnot(all(shape == p), identical(unname(lengths(rows)), seq_len(p)))
  V <- matrix(0, p, p)
  for (i in seq_len(p)) V[i, seq_len(i)] <- rows[[i]]
  V[upper.tri(V)] <- t(V)[upper.tri(V)]
  V
}

logitrs_names <- c("(Intercept)", "x", "z", "log_sd_intercept", "log_sd_slope", "atanh_rho")
logitrs_transform <- function(b, index) {
  # index: intercept, x, z, var(intercept), var(slope), covariance.
  A <- b[index[4L]]; B <- b[index[5L]]; C <- b[index[6L]]; rho <- C/sqrt(A*B)
  stopifnot(A > 0, B > 0, abs(rho) < 1)
  J <- matrix(0, 6, length(b)); J[cbind(1:3, index[1:3])] <- 1
  J[4, index[4]] <- 1/(2*A); J[5, index[5]] <- 1/(2*B)
  J[6, index[4:6]] <- c(-rho/(2*A), -rho/(2*B), 1/sqrt(A*B))/(1-rho^2)
  list(coefficients = setNames(unname(c(b[index[1:3]], .5*log(A), .5*log(B), atanh(rho))), logitrs_names),
    jacobian = J)
}

logitrs_component <- function(block) {
  margins <- grep("^MARGINAL (MEAN|SLOPE X):", block); stopifnot(length(margins) == 2L)
  parameters <- block[seq_len(margins[1L]-1L)]
  b <- logitrs_read_matrix(parameters, "e(b)"); V <- logitrs_read_matrix(parameters, "e(V)")
  stopifnot(length(b) == 6L, nrow(V) == 6L, all(is.finite(b)), all(is.finite(V)))
  raw_names <- c("x", "z", "_cons", "var_slope", "var_intercept", "cov_intercept_slope")
  names(b) <- raw_names; dimnames(V) <- list(raw_names, raw_names)
  trans <- logitrs_transform(b, c(3, 1, 2, 5, 4, 6))
  W <- trans$jacobian %*% V %*% t(trans$jacobian); dimnames(W) <- list(logitrs_names, logitrs_names)
  stopifnot(min(eigen(W, symmetric = TRUE)$values) > 0)
  status <- grep("^CONVERGED=1 N=", block, value = TRUE); stopifnot(length(status) == 1L)
  list(coefficients = trans$coefficients, vcov = W,
    logLik = as.numeric(sub(".*LL= *", "", status)),
    nobs = as.integer(sub(".* N=([0-9]+).*", "\\1", status)),
    n_groups = as.integer(logitrs_read_matrix(parameters, "e(N_g)")[1, 1]),
    stata_raw = list(coefficients = b, vcov = V), jacobian = trans$jacobian,
    native_margins = setNames(lapply(1:2, function(i) {
      end <- if (i == 1L) margins[2L]-1L else length(block)
      m <- block[margins[i]:end]
      value <- c(estimate = logitrs_read_matrix(m, "r(b)"), variance = logitrs_read_matrix(m, "r(V)")[1, 1])
      stopifnot(all(is.finite(value)), value["variance"] >= 0)
      value
    }), c("mean", "slope_x")))
}

logitrs_joint <- function(block, outcomes, nobs, n_clusters) {
  b <- logitrs_read_matrix(block, "e(b)", TRUE); V <- logitrs_read_matrix(block, "e(V)", TRUE)
  B <- logitrs_read_matrix(block, "e(V_modelbased)", TRUE)
  stopifnot(length(b) == 16L, nrow(V) == 16L, nrow(B) == 16L,
    all(is.finite(b)), all(is.finite(V)), all(is.finite(B)),
    all(b[c(3, 4, 8, 9)] == 1), all(V[c(3, 4, 8, 9), ] == 0), all(B[c(3, 4, 8, 9), ] == 0))
  raw_names <- c(paste0(outcomes[1], ":", c("x", "z", "I1[id]", "c.x#S1[id]", "_cons")),
    paste0(outcomes[2], ":", c("x", "z", "I2[id]", "c.x#S2[id]", "_cons")),
    "/:var(I1[id])", "/:var(S1[id])", "/:var(I2[id])", "/:var(S2[id])",
    "/:cov(I1[id],S1[id])", "/:cov(I2[id],S2[id])")
  names(b) <- raw_names; dimnames(V) <- dimnames(B) <- list(raw_names, raw_names)
  a <- logitrs_transform(b, c(5, 1, 2, 11, 12, 15))
  z <- logitrs_transform(b, c(10, 6, 7, 13, 14, 16))
  converted <- setNames(c(a$coefficients, z$coefficients),
    unlist(lapply(c("Y1", "Y2"), function(y) paste0(y, "::", logitrs_names))))
  J <- rbind(a$jacobian, z$jacobian); dimnames(J) <- list(names(converted), raw_names)
  W <- J %*% V %*% t(J); Bc <- J %*% B %*% t(J)
  stopifnot(min(eigen(W, symmetric = TRUE)$values) > 0, min(eigen(Bc, symmetric = TRUE)$values) > 0)
  status <- grep(paste0("^CONVERGED=1 N=", nobs, " N_clust=", n_clusters, " LL="), block, value = TRUE)
  stopifnot(length(status) == 1L)
  list(coefficients = converted, vcov = W, modelbased = Bc,
    logLik = as.numeric(sub(".*LL= *", "", status)), nobs = as.integer(nobs), n_clusters = as.integer(n_clusters),
    jacobian = J, stata_raw = list(coefficients = b, vcov = V, modelbased = B,
      gradient = logitrs_read_matrix(block, "e(gradient)", TRUE)))
}
