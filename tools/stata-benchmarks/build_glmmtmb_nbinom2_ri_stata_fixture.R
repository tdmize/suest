# Run from tools/stata-benchmarks, passing the returned revision-2 Stata log.
args <- commandArgs(trailingOnly = TRUE)
log_file <- if (length(args)) args[1] else
  "suest_r_glmmtmb_nbinom2_ri_crosslang_benchmark_v2.log"
lines <- readLines(log_file, warn = FALSE)
stopifnot(all(c("BENCHMARK_FILE_REVISION=2",
  "NBINOM2_PRIMARY_LAPLACE_BENCHMARK_COMPLETE=1",
  "GLMMTMB_NBINOM2_RI_CROSSLANG_BENCHMARK_COMPLETE=1") %in% lines))
starts <- grep("^=== .* ===$", lines)
sections <- setNames(lapply(seq_along(starts), function(i) {
  end <- if (i < length(starts)) starts[i+1L]-1L else length(lines)
  lines[starts[i]:end]
}), sub("^=== (.*) ===$", "\\1", lines[starts]))
stopifnot(!anyDuplicated(names(sections)))

# Stata wraps wide matrices by columns. Reassemble each labelled row before
# checking its lower-triangle length; standalone equation labels disambiguate
# repeated coefficient names in the adaptive suest2 matrices.
read_matrix <- function(block, name) {
  headers <- grep("^(symmetric )?[er]\\([bV]\\)\\[", block)
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
    key <- if (grepl(":", token[1L])) token[1L] else paste0(equation, token[1L])
    rows[[key]] <- c(rows[[key]], values)
  }
  shape <- as.integer(strsplit(sub(".*\\[([0-9,]+)\\].*", "\\1",
    block[start]), ",", fixed = TRUE)[[1L]])
  if (name %in% c("e(b)", "r(b)")) {
    stopifnot(length(rows) == 1L, length(rows[[1L]]) == shape[2L])
    return(unname(rows[[1L]]))
  }
  p <- length(rows)
  stopifnot(all(shape == p), identical(unname(lengths(rows)), seq_len(p)))
  V <- matrix(0, p, p)
  for (i in seq_len(p)) V[i, seq_len(i)] <- rows[[i]]
  V[upper.tri(V)] <- t(V)[upper.tri(V)]
  V
}

convert <- function(block, type) {
  b <- read_matrix(block, "e(b)"); V <- read_matrix(block, "e(V)")
  stopifnot(length(b) == nrow(V), all(is.finite(b)), all(is.finite(V)))
  raw <- list(coefficients = b, vcov = V)
  if (type == "component") {
    stopifnot(length(b) == 5L)
    dispersion <- 4L; variance <- 5L; order <- c(3, 1, 2, 4, 5)
  } else if (type == "gsem") {
    stopifnot(length(b) == 12L, all(b[c(3, 7)] == 1), all(V[c(3, 7), ] == 0))
    dispersion <- c(9L, 10L); variance <- c(11L, 12L)
    order <- c(4, 1, 2, 9, 11, 8, 5, 6, 10, 12)
  } else {
    stopifnot(type == "adaptive", length(b) == 10L)
    dispersion <- c(4L, 9L); variance <- c(5L, 10L)
    order <- c(3, 1, 2, 4, 5, 8, 6, 7, 9, 10)
  }
  stopifnot(all(b[variance] > 0))
  J <- diag(length(b)); diag(J)[dispersion] <- -1
  diag(J)[variance] <- 1/(2*b[variance])
  V <- (J %*% V %*% J)[order, order]
  b[dispersion] <- -b[dispersion]; b[variance] <- .5*log(b[variance])
  b <- b[order]
  parameter_names <- c("(Intercept)", "x", "z", "log_phi", "log_sigma")
  if (type != "component") parameter_names <- unlist(lapply(c("Y1", "Y2"),
    function(s) paste0(s, "::", parameter_names)))
  names(b) <- parameter_names; dimnames(V) <- list(parameter_names, parameter_names)
  stopifnot(all(eigen(V, symmetric = TRUE, only.values = TRUE)$values > 0))
  nline <- grep("^N=", block, value = TRUE)
  stopifnot(length(nline) == 1L)
  ll <- if (grepl("LL=", nline)) as.numeric(sub(".*LL= *", "", nline)) else NA_real_
  list(coefficients = b, vcov = V, logLik = ll,
    nobs = as.integer(sub("^N=([0-9]+).*", "\\1", nline)),
    n_clusters = if (grepl("N_clust=", nline))
      as.integer(sub(".*N_clust=([0-9]+).*", "\\1", nline)) else NULL,
    stata_raw = raw)
}

outcomes <- c("y1", "y2", "y1_partial", "y2_partial", "y1_left", "y2_right")
components <- setNames(lapply(outcomes, function(y) {
  block <- sections[[paste0("LAPLACE COMPONENT: ", y)]]
  result <- convert(block, "component")
  if (y %in% c("y1", "y2")) {
    margin_starts <- grep("^MARGINAL (MEAN|SLOPE X):", block)
    stopifnot(length(margin_starts) == 2L)
    margins <- lapply(seq_along(margin_starts), function(i) {
      end <- if (i == 1L) margin_starts[2L]-1L else length(block)
      margin <- block[margin_starts[i]:end]
      c(estimate = read_matrix(margin, "r(b)"), variance = read_matrix(margin, "r(V)")[1, 1])
    })
    result$native_margins <- list(
      estimates = setNames(vapply(margins, `[[`, numeric(1), "estimate"), c("mean", "slope_x")),
      variances = setNames(vapply(margins, `[[`, numeric(1), "variance"), c("mean", "slope_x")))
  }
  result
}), outcomes)
gsem_names <- c(balanced = "BALANCED JOINT GSEM, LAPLACE",
  partial_higher = "PARTIAL OVERLAP, HIGHER-CLUSTER GSEM LAPLACE",
  disjoint = "DISJOINT-GROUP GSEM, LAPLACE")
gsem <- lapply(gsem_names, function(s) convert(sections[[s]], "gsem"))
adaptive_names <- c(balanced = "balanced", partial_higher = "partial", disjoint = "disjoint")
adaptive <- lapply(adaptive_names, function(s) {
  block <- sections[[paste0("OPTIONAL ADAPTIVE 12-POINT SUEST2: ", s)]]
  stopifnot("SUEST2_ADAPTIVE_RC=0" %in% block)
  convert(block, "adaptive")
})
stopifnot("SUEST2_LAPLACE_RC=322" %in% sections[["OPTIONAL SUEST2 LAPLACE CHECK"]])
dat <- read.csv("suest_r_glmmtmb_nbinom2_ri_crosslang_benchmark_v2.csv")
fixture <- list(data = dat, components = components, gsem = gsem, adaptive = adaptive,
  metadata = list(Stata = "19.5", revision = 2L, source_log = basename(log_file),
    source_md5 = unname(tools::md5sum(log_file)), laplace_suest2_return_code = 322L,
    disjoint_groups = "odd versus even id",
    primary_integration = "Laplace", secondary_integration = "adaptive 12-point",
    gsem_correction = c(balanced = 80/79, partial_higher = 16/15, disjoint = 80/79),
    suest_correction = c(balanced = 80/79, partial_higher = 16/15, disjoint = 40/39),
    nobs_overlap = c(balanced = 480L, partial_higher = 444L, disjoint = 0L)))
saveRDS(fixture, "../../tests/testthat/fixtures/glmmtmb-nbinom2-ri-stata.rds")
cat("NBINOM2_REVISION2_STATA_FIXTURE_COMPLETE=1\n")
