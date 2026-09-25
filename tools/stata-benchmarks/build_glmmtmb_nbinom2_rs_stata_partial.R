# Run here with the returned log and optional output RDS path.
source("nbinom2_rs_stata_log_helpers.R")
args <- commandArgs(trailingOnly = TRUE)
log_file <- args[1L]; lines <- readLines(log_file, warn = FALSE)
stopifnot("BENCHMARK_FILE_REVISION=2" %in% lines, "NB2_RS_PARTIAL_HIGHER_COMPLETE=1" %in% lines,
  all(paste0("NB2_RS_COMPONENT_", c("y1_partial", "y2_partial"), "_COMPLETE=1") %in% lines))
components <- setNames(lapply(c("y1_partial", "y2_partial"), function(y) {
  start <- which(lines == paste0("=== LAPLACE COMPONENT: ", y, " ==="))
  end <- which(lines == paste0("NB2_RS_COMPONENT_", y, "_COMPLETE=1"))
  stopifnot(length(start) == 1L, length(end) == 1L, end > start)
  nb2rs_component(lines[start:end])
}), c("y1_partial", "y2_partial"))
start <- which(lines == "=== PARTIAL JOINT GSEM LAPLACE (MAXIMUM 15 ITERATIONS) ===")
stopifnot(length(start) == 1L)
joint <- nb2rs_joint(lines[start:length(lines)], c("y1_partial", "y2_partial"), 1000, 20)
fixture <- list(data = read.csv("suest_r_glmmtmb_nbinom2_rs_crosslang_benchmark.csv"),
  components = components, joint = joint, metadata = list(Stata = "19.5", revision = 2L,
    source_log = basename(log_file), source_md5 = unname(tools::md5sum(log_file)),
    integration = "Laplace", joint_adapt_tolerance = 1e-12,
    validation = "raw returned results; assess numerical agreement separately"))
saveRDS(fixture, if (length(args) > 1L) args[2L] else "../../tests/testthat/fixtures/glmmtmb-nbinom2-rs-stata-partial-v2.rds")
cat("NB2_RS_PARTIAL_FIXTURE_COMPLETE=1\n")
