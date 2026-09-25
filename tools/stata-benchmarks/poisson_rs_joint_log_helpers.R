# Shared parser for the returned Poisson random-slope joint logs.
read_poisson_rs_joint <- function(lines, case, outcomes, n_clusters, revision,
  source_log, source_md5) {
  start <- grep(paste0("^=== ", toupper(case), " JOINT LAPLACE FROM COMPONENT FITS"), lines)
  stopifnot(length(start) == 1L)
  end <- grep("^=== SAMPLE CASE:", lines)
  end <- end[end > start]
  block <- lines[start:if (length(end)) end[1L]-1L else length(lines)]
  # Concatenate wrapped column panels using equation plus printed row label.
  # In this joint model the two wrapped covariance labels end in distinct names.
  read_matrix <- function(name, optional = FALSE) {
    headers <- grep("^(symmetric )?e\\([[:alnum:]_]+\\)\\[", block)
    start <- headers[grepl(paste0(name, "["), block[headers], fixed = TRUE)]
    if (!length(start) && optional) return(NULL)
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
      key <- if (grepl(":", token[1L])) token[1L] else paste0(equation, token[1L])
      rows[[key]] <- c(rows[[key]], values)
    }
    shape <- as.integer(strsplit(sub(".*\\[([0-9,]+)\\].*", "\\1", block[start]), ",")[[1L]])
    if (shape[1L] == 1L) {
      stopifnot(length(rows) == 1L, length(rows[[1L]]) == shape[2L])
      return(rows[[1L]])
    }
    p <- length(rows)
    stopifnot(all(shape == p), identical(unname(lengths(rows)), seq_len(p)))
    V <- matrix(0, p, p)
    for (i in seq_len(p)) V[i, seq_len(i)] <- rows[[i]]
    V[upper.tri(V)] <- t(V)[upper.tri(V)]
    V
  }
  b <- read_matrix("e(b)"); V <- read_matrix("e(V)")
  stopifnot(length(b) == 16L, nrow(V) == 16L, all(is.finite(b)), all(is.finite(V)),
    all(b[c(3, 4, 8, 9)] == 1), all(V[c(3, 4, 8, 9), ] == 0))
  native_names <- c(paste0(outcomes[1], ":", c("x", "z", "I1[id]", "c.x#S1[id]", "_cons")),
    paste0(outcomes[2], ":", c("x", "z", "I2[id]", "c.x#S2[id]", "_cons")),
    "/:var(I1[id])", "/:var(S1[id])", "/:var(I2[id])", "/:var(S2[id])",
    "/:cov(I1[id],S1[id])", "/:cov(I2[id],S2[id])")
  names(b) <- native_names; dimnames(V) <- list(native_names, native_names)
  J <- matrix(0, 12, 16); converted <- numeric(12)
  for (i in 1:2) {
    idx <- if (i == 1L) c(5, 1, 2, 11, 12, 15) else c(10, 6, 7, 13, 14, 16)
    rows <- (i-1L)*6L + 1:6
    A <- b[idx[4]]; B <- b[idx[5]]; C <- b[idx[6]]; rho <- C/sqrt(A*B)
    stopifnot(A > 0, B > 0, abs(rho) < 1)
    converted[rows] <- c(b[idx[1:3]], .5*log(A), .5*log(B), atanh(rho))
    J[cbind(rows[1:3], idx[1:3])] <- 1
    J[rows[4], idx[4]] <- 1/(2*A); J[rows[5], idx[5]] <- 1/(2*B)
    J[rows[6], idx[4:6]] <- c(-rho/(2*A), -rho/(2*B), 1/sqrt(A*B))/(1-rho^2)
  }
  names(converted) <- unlist(lapply(c("Y1", "Y2"), function(y)
    paste0(y, "::", c("(Intercept)", "x", "z", "log_sd_intercept", "log_sd_slope", "atanh_rho"))))
  dimnames(J) <- list(names(converted), native_names)
  W <- J %*% V %*% t(J)
  stopifnot(min(eigen(W, symmetric = TRUE)$values) > 0)
  status <- grep(paste0("^CONVERGED=1 N=800 N_clust=", n_clusters, " LL="), block, value = TRUE)
  stopifnot(length(status) == 1L)
  V_modelbased <- read_matrix("e(V_modelbased)", optional = TRUE)
  fixture <- list(coefficients = converted, vcov = W,
    logLik = as.numeric(sub(".*LL= *", "", status)), nobs = 800L, n_clusters = as.integer(n_clusters),
    jacobian = J, stata_raw = list(coefficients = b, vcov = V,
      gradient = read_matrix("e(gradient)"), modelbased = V_modelbased),
    modelbased = if (!is.null(V_modelbased)) J %*% V_modelbased %*% t(J) else NULL,
    metadata = list(Stata = "19.5", revision = revision, source_log = basename(source_log),
      source_md5 = source_md5, fitting_integration = "Laplace",
      adapt_tolerance = if (revision == 3L) 1e-8 else 1e-12,
      covariance_validation = "pending; a completed Stata run does not establish covariance agreement"))
  fixture
}
