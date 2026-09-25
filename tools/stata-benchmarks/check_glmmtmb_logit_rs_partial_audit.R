# Run here; optional output RDS. Independently audit LOGIT Laplace curvature/scores.
source("logit_rs_joint_audit_helpers.R")
f <- readRDS("../../tests/testthat/fixtures/glmmtmb-logit-rs-stata-partial-v1.rds")
r <- readRDS("glmmtmb_logit_rs_r_reference.rds")
result <- setNames(lapply(c("Stata", "R"), function(point) {
  cat("=== EVALUATION POINT:", point, "===\n")
  audit_logit_rs_joint(f$joint, r$systems$partial_higher, f$data,
    c("y1_partial", "y2_partial"), "higher", f$components, evaluate_at = point)
}), c("Stata", "R"))
for (i in 1:2) {
  y <- c("y1_partial", "y2_partial")[i]; B <- result$R$components[[i]]$bread_native_delta
  cat(y, "independent_vs_R_native_max_covariance=", max(abs(B-r$components[[y]]$vcov)),
    "relative_SE=", max(abs(sqrt(diag(B)/diag(r$components[[y]]$vcov))-1)), "\n")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) saveRDS(result, args[1L])
# Retain the earlier strict covariance/SE tolerances, now at matching points.
metrics <- rbind(result$R$comparison["sandwich_independent_vs_R", ],
  result$R$comparison["native_coordinate_sandwich_vs_R", ],
  result$Stata$comparison["stata_bread_independent_meat_vs_stata", ])
stopifnot(all(metrics[, "max_abs_covariance"] < 1e-6),
  all(metrics[, "max_relative_diagonal_sd"] < 1e-4))
saveRDS(result, "../../tests/testthat/fixtures/glmmtmb-logit-rs-partial-audit.rds")
cat("LOGIT_RS_PARTIAL_INDEPENDENT_AUDIT_COMPLETE=1\n")
