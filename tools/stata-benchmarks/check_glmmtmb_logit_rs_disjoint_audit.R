# Run here; optional output RDS. Independently audit LOGIT Laplace curvature/scores.
source("logit_rs_joint_audit_helpers.R")
f <- readRDS("../../tests/testthat/fixtures/glmmtmb-logit-rs-stata-disjoint-v1.rds")
r <- readRDS("glmmtmb_logit_rs_r_reference.rds")
result <- setNames(lapply(c("Stata", "R"), function(point) {
  cat("=== EVALUATION POINT:", point, "===\n")
  audit_logit_rs_joint(f$joint, r$systems$disjoint, f$data,
    c("y1_left", "y2_right"), "id", f$components, evaluate_at = point)
}), c("Stata", "R"))
for (i in 1:2) {
  y <- c("y1_left", "y2_right")[i]; B <- result$R$components[[i]]$bread_native_delta
  cat(y, "independent_vs_R_native_max_covariance=", max(abs(B-r$components[[y]]$vcov)),
    "relative_SE=", max(abs(sqrt(diag(B)/diag(r$components[[y]]$vcov))-1)), "\n")
}
stopifnot(all(result$R$metadata$component_clusters == 50L),
  result$R$metadata$n_clusters == 100L,
  abs(result$R$metadata$normalization-(50/49)/(100/99)) < 1e-14,
  max(abs(result$R$native_coordinate_vcov[1:6, 7:12])) == 0)
# Preserve the original bounds and distinguish R validation from raw Stata diagnostics.
metrics <- result$R$comparison[c("sandwich_independent_vs_R", "native_coordinate_sandwich_vs_R"), ]
raw <- result$Stata$comparison["stata_bread_independent_meat_vs_stata", ]
result$validation <- list(absolute_bound = 1e-6, relative_se_bound = 1e-4,
  independent_R_pass = all(metrics[, "max_abs_covariance"] < 1e-6) &&
    all(metrics[, "max_relative_diagonal_sd"] < 1e-4),
  raw_reconstruction_absolute_pass = unname(raw["max_abs_covariance"] < 1e-6),
  raw_reconstruction_relative_se_pass = unname(raw["max_relative_diagonal_sd"] < 1e-4))
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) saveRDS(result, args[1L])
stopifnot(result$validation$independent_R_pass)
saveRDS(result, "../../tests/testthat/fixtures/glmmtmb-logit-rs-disjoint-audit.rds")
cat("LOGIT_RS_DISJOINT_RAW_RECONSTRUCTION_ABSOLUTE_PASS=",
  as.integer(result$validation$raw_reconstruction_absolute_pass), "\n", sep = "")
cat("LOGIT_RS_DISJOINT_RAW_RECONSTRUCTION_RELATIVE_SE_PASS=",
  as.integer(result$validation$raw_reconstruction_relative_se_pass), "\n", sep = "")
cat("LOGIT_RS_DISJOINT_INDEPENDENT_AUDIT_COMPLETE=1\n")
