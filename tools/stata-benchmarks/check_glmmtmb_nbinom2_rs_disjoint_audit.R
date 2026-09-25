# Run here; optional output RDS. Independently audit NB2 Laplace curvature/scores.
source("nbinom2_rs_joint_audit_helpers.R")
f <- readRDS("../../tests/testthat/fixtures/glmmtmb-nbinom2-rs-stata-disjoint-v3.rds")
r <- readRDS("glmmtmb_nbinom2_rs_r_reference.rds")
result <- setNames(lapply(c("Stata", "R"), function(point) {
  cat("=== EVALUATION POINT:", point, "===\n")
  audit_nbinom2_rs_joint(f$joint, r$systems$disjoint, f$data,
    c("y1_left", "y2_right"), "id", f$components, evaluate_at = point)
}), c("Stata", "R"))
for (i in 1:2) {
  y <- c("y1_left", "y2_right")[i]; B <- result$R$components[[i]]$bread_native_delta
  cat(y, "independent_vs_R_native_max_covariance=", max(abs(B-r$components[[y]]$vcov)),
    "relative_SE=", max(abs(sqrt(diag(B)/diag(r$components[[y]]$vcov))-1)), "\n")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) saveRDS(result, args[1L])
stopifnot(all(result$R$metadata$component_clusters == 50L),
  abs(result$R$metadata$normalization-99/98) < 1e-14,
  max(abs(result$R$independent_meat[1:7, 8:14])) == 0,
  max(abs(result$R$independent_vcov[1:7, 8:14])) == 0)
# Match both the parameter point and native information coordinates. The
# unadjusted Fisher-coordinate audit remains recorded separately, including its discrepancy.
metrics <- rbind(result$R$comparison["native_coordinate_sandwich_vs_R", ],
  result$Stata$comparison["stata_bread_independent_meat_vs_stata", ])
stopifnot(all(metrics[, "max_abs_covariance"] < 1e-6),
  all(metrics[, "max_relative_diagonal_sd"] < 1e-4))
saveRDS(result, "../../tests/testthat/fixtures/glmmtmb-nbinom2-rs-disjoint-audit.rds")
cat("NB2_RS_DISJOINT_INDEPENDENT_AUDIT_COMPLETE=1\n")
