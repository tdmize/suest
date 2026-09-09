# Build the permanent survey cross-language fixture from the returned Stata log.
# Run from the package root:
# Rscript tools/stata-benchmarks/build_survey_stata_fixture.R path/to/log

args <- commandArgs(trailingOnly = TRUE)
log_file <- if (length(args)) args[1L] else
  "tools/stata-benchmarks/suest_r_survey_linear_benchmark_v2.log"
lines <- sub("\r$", "", readLines(log_file, warn = FALSE))

read_matrix <- function(block, marker) {
  start <- grep(paste0("^", marker), block)[1L]
  if (is.na(start)) return(NULL)
  block <- block[start:length(block)]
  b_start <- grep("^e\\(b\\)\\[1,[0-9]+\\]$", block)[1L]
  v_start <- grep("^symmetric e\\(V\\)\\[[0-9]+,[0-9]+\\]$", block)[1L]
  n <- as.integer(sub("^e\\(b\\)\\[1,([0-9]+)\\]$", "\\1", block[b_start]))
  b_line <- grep("^y1 +", block[seq.int(b_start + 1L, v_start - 1L)],
    value = TRUE)[1L]
  b <- as.numeric(strsplit(trimws(sub("^y1", "", b_line)), " +")[[1L]])
  rows <- list()
  for (line in block[seq.int(v_start + 1L, length(block))]) {
    tokens <- strsplit(trimws(line), " +")[[1L]]
    if (length(tokens) < 2L) next
    values <- suppressWarnings(as.numeric(tokens[-1L]))
    if (anyNA(values) || length(values) != length(rows) + 1L) next
    rows[[length(rows) + 1L]] <- values
    if (length(rows) == n) break
  }
  stopifnot(length(b) == n, length(rows) == n)
  V <- matrix(0, n, n)
  for (i in seq_len(n)) V[i, seq_len(i)] <- rows[[i]]
  V <- V + t(V) - diag(diag(V))
  list(b = b, V = V)
}

case_starts <- grep("^SURVEY_CASE=[1-6]$", lines)
reference <- lapply(seq_along(case_starts), function(i) {
  end <- if (i < length(case_starts)) case_starts[i + 1L] - 1L else
    grep("^SURVEY_LINEAR_BENCHMARK_V2_COMPLETE=1$", lines)[1L] - 1L
  block <- lines[case_starts[i]:end]
  native_a <- read_matrix(block, "NATIVE_MODEL=A")
  native_b <- read_matrix(block, "NATIVE_MODEL=B")
  stacked <- read_matrix(block, "STACKED_JOINT_MODEL")
  rc <- as.integer(sub("^SUEST2_RC=", "", grep("^SUEST2_RC=", block,
    value = TRUE)[1L]))
  joint <- if (rc == 0L) read_matrix(block, "JOINT_MODEL") else NULL
  order_a <- c(2L, 1L)
  order_b <- c(3L, 1L, 2L)
  order_joint <- c(2L, 1L, 5L, 3L, 4L)
  list(
    b = stacked$b,
    V = stacked$V,
    native_a = native_a$V[order_a, order_a, drop = FALSE],
    native_b = native_b$V[order_b, order_b, drop = FALSE],
    suest2_rc = rc,
    suest2_b = if (rc == 0L) joint$b[order_joint] else NULL,
    suest2_V = if (rc == 0L)
      joint$V[order_joint, order_joint, drop = FALSE] else NULL
  )
})

fixture <- list(
  data = read.csv("tools/stata-benchmarks/suest_r_survey_linear_benchmark.csv"),
  reference = reference,
  source = list(
    stata = "19.5",
    suest2 = "1.0.0 30aug2026",
    benchmark_revision = 2L,
    returned_log = basename(log_file),
    coefficient_order = c("A::(Intercept)", "A::x", "B::(Intercept)",
      "B::x", "B::z")
  )
)
saveRDS(fixture, "tests/testthat/fixtures/survey-stata.rds", version = 3)
