# Run here with the returned log and optional output RDS path.
source("logit_rs_stata_log_helpers.R")
args <- commandArgs(trailingOnly = TRUE)
log_file <- args[1L]; lines <- readLines(log_file, warn = FALSE)
stopifnot("BENCHMARK_FILE_REVISION=1" %in% lines, "LOGIT_RS_BALANCED_COMPLETE=1" %in% lines,
  all(paste0("LOGIT_RS_COMPONENT_", c("y1", "y2"), "_COMPLETE=1") %in% lines))
gaps <- as.numeric(sub("^INITIAL_LOGLIK_GAP= *", "", grep("^INITIAL_LOGLIK_GAP=", lines, value = TRUE)))
changes <- as.numeric(sub("^MAX_START_PARAMETER_CHANGE= *", "", grep("^MAX_START_PARAMETER_CHANGE=", lines, value = TRUE)))
stopifnot(length(gaps) == 3L, length(changes) == 3L, all(is.finite(gaps)),
  max(abs(gaps)) < 1e-4, max(abs(changes)) <= 1e-10,
  sum(lines == "COMPONENT_FIT_RC=0") == 2L, "BALANCED_FIT_RC=0" %in% lines,
  !any(grepl("^[A-Z_]+_FAILED:", lines)))
components <- setNames(lapply(c("y1", "y2"), function(y) {
  start <- which(lines == paste0("=== LAPLACE COMPONENT: ", y, " ==="))
  end <- which(lines == paste0("LOGIT_RS_COMPONENT_", y, "_COMPLETE=1"))
  stopifnot(length(start) == 1L, length(end) == 1L, end > start)
  result <- logitrs_component(lines[start:end])
  stopifnot(result$nobs == 1200L, result$n_groups == 100L)
  result
}), c("y1", "y2"))
start <- which(lines == "=== BALANCED JOINT GSEM LAPLACE (MAXIMUM 15 ITERATIONS) ===")
stopifnot(length(start) == 1L)
joint <- logitrs_joint(lines[start:length(lines)], c("y1", "y2"), 1200, 100)
fixture <- list(data = read.csv("suest_r_glmmtmb_logit_rs_crosslang_benchmark.csv"),
  components = components, joint = joint, metadata = list(Stata = "19.5", revision = 1L,
    source_log = basename(log_file), source_md5 = unname(tools::md5sum(log_file)),
    integration = "Laplace", zero_step_loglik_gaps = gaps, zero_step_parameter_changes = changes, joint_adapt_tolerance = 1e-12,
    validation = "raw returned results; assess numerical agreement separately"))
saveRDS(fixture, if (length(args) > 1L) args[2L] else "../../tests/testthat/fixtures/glmmtmb-logit-rs-stata-balanced-v1.rds")
cat("LOGIT_RS_BALANCED_FIXTURE_COMPLETE=1\n")
