# Run here with the returned revision-2 log; this records components only.
args <- commandArgs(trailingOnly = TRUE)
log_file <- if (length(args)) args[1] else "suest_r_glmmtmb_poisson_rs_crosslang_benchmark_v2.log"
lines <- readLines(log_file, warn = FALSE)
stopifnot("BENCHMARK_FILE_REVISION=2" %in% lines,
  "POISSON_RS_COMPONENT_BENCHMARK_COMPLETE=1" %in% lines)
starts <- grep("^=== .* ===$", lines)
sections <- setNames(lapply(seq_along(starts), function(i) {
  end <- if (i < length(starts)) starts[i+1L]-1L else length(lines)
  lines[starts[i]:end]
}), sub("^=== (.*) ===$", "\\1", lines[starts]))
stopifnot(!anyDuplicated(names(sections)))

# These six-column matrices fit in one panel at linesize 255. Preserve row
# order: Stata wraps both variance and covariance labels onto _cons[id]).
read_matrix <- function(block, name) {
  headers <- grep("^(symmetric )?[er]\\([bV]\\)\\[", block)
  start <- headers[grepl(paste0(name, "["), block[headers], fixed = TRUE)]
  stopifnot(length(start) == 1L)
  next_header <- headers[headers > start]
  end <- if (length(next_header)) next_header[1L]-1L else length(block)
  rows <- list()
  for (line in block[(start+1L):end]) {
    token <- strsplit(trimws(line), "[[:space:]]+")[[1L]]
    if (length(token) < 2L) next
    values <- suppressWarnings(as.numeric(token[-1L]))
    if (!anyNA(values)) rows[[length(rows)+1L]] <- values
  }
  shape <- as.integer(strsplit(sub(".*\\[([0-9,]+)\\].*", "\\1", block[start]), ",")[[1L]])
  if (name %in% c("e(b)", "r(b)")) {
    stopifnot(shape[1L] == 1L, length(rows) == 1L, length(rows[[1L]]) == shape[2L])
    return(rows[[1L]])
  }
  p <- length(rows)
  stopifnot(all(shape == p), identical(lengths(rows), seq_len(p)))
  V <- matrix(0, p, p)
  for (i in seq_len(p)) V[i, seq_len(i)] <- rows[[i]]
  V[upper.tri(V)] <- t(V)[upper.tri(V)]
  V
}

outcomes <- c("y1", "y2", "y1_partial", "y2_partial", "y1_left", "y2_right")
components <- setNames(lapply(outcomes, function(y) {
  attempts <- sections[[paste0("FIT ATTEMPTS: ", y)]]
  stopifnot("FINAL_LAPLACE_ACCEPTED=1 START_ROUTE=adaptive_start" %in% attempts)
  block <- sections[[paste0("LAPLACE COMPONENT: ", y)]]
  margin_starts <- grep("^MARGINAL (MEAN|SLOPE X):", block)
  parameter_block <- if (length(margin_starts)) block[seq_len(margin_starts[1]-1L)] else block
  b <- read_matrix(parameter_block, "e(b)"); V <- read_matrix(parameter_block, "e(V)")
  stopifnot(length(b) == 6L, nrow(V) == 6L, all(is.finite(b)), all(is.finite(V)))
  # Returned stripes: x,z,_cons,var(x[id]),var(_cons[id]),cov(x[id],_cons[id]).
  raw_names <- c("x", "z", "_cons", "var_slope", "var_intercept", "cov_intercept_slope")
  names(b) <- raw_names; dimnames(V) <- list(raw_names, raw_names)
  raw <- list(coefficients = b, vcov = V)
  A <- b[5]; B <- b[4]; C <- b[6]; rho <- C/sqrt(A*B)
  stopifnot(A > 0, B > 0, abs(rho) < 1)
  J <- matrix(0, 6, 6); J[1, 3] <- J[2, 1] <- J[3, 2] <- 1
  J[4, 5] <- 1/(2*A); J[5, 4] <- 1/(2*B)
  J[6, c(5, 4, 6)] <- c(-rho/(2*A), -rho/(2*B), 1/sqrt(A*B))/(1-rho^2)
  converted <- unname(c(b[c(3, 1, 2)], .5*log(A), .5*log(B), atanh(rho)))
  names(converted) <- c("(Intercept)", "x", "z", "log_sd_intercept", "log_sd_slope", "atanh_rho")
  V <- J %*% V %*% t(J); dimnames(V) <- list(names(converted), names(converted))
  stopifnot(min(eigen(V, symmetric = TRUE)$values) > 0)
  nline <- grep("^N=", block, value = TRUE); stopifnot(length(nline) == 1L)
  result <- list(coefficients = converted, vcov = V,
    logLik = as.numeric(sub(".*LL= *", "", nline)),
    nobs = as.integer(sub("^N=([0-9]+).*", "\\1", nline)), stata_raw = raw)
  if (length(margin_starts)) {
    stopifnot(length(margin_starts) == 2L)
    result$native_margins <- setNames(lapply(seq_along(margin_starts), function(i) {
      end <- if (i == 1L) margin_starts[2L]-1L else length(block)
      m <- block[margin_starts[i]:end]
      stopifnot(any(grepl("^MARGINAL_.*_ACCEPTED=1$", m)))
      c(estimate = read_matrix(m, "r(b)"), variance = read_matrix(m, "r(V)")[1, 1])
    }), c("mean", "slope_x"))
  }
  result
}), outcomes)
fixture <- list(data = read.csv("suest_r_glmmtmb_poisson_rs_crosslang_benchmark.csv"),
  components = components, metadata = list(Stata = "19.5", revision = 2L,
    source_log = basename(log_file), source_md5 = unname(tools::md5sum(log_file)),
    fitting_integration = "Laplace", initialization = "Stata adaptive seven-point fit",
    joint_validation = "pending: balanced cold attempt interrupted"))
saveRDS(fixture, "../../tests/testthat/fixtures/glmmtmb-poisson-rs-stata-components.rds")
cat("POISSON_RS_STATA_COMPONENT_FIXTURE_COMPLETE=1\n")
