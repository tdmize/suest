# Run here after building the returned balanced-v4 and remaining-v5 fixtures.
# Optional argument: output RDS. Raw GSEM covariance is diagnostic evidence.
source("poisson_rs_joint_audit_helpers.R")
f <- readRDS("../../tests/testthat/fixtures/glmmtmb-poisson-rs-stata-components.rds")
s <- c(list(balanced = readRDS("../../tests/testthat/fixtures/glmmtmb-poisson-rs-stata-balanced-v4.rds")),
  readRDS("../../tests/testthat/fixtures/glmmtmb-poisson-rs-stata-remaining-v5.rds"))
r <- readRDS("glmmtmb_poisson_rs_r_reference.rds")$systems
outcomes <- list(balanced = c("y1", "y2"), partial_higher = c("y1_partial", "y2_partial"),
  disjoint = c("y1_left", "y2_right"))
results <- setNames(lapply(names(s), function(case) {
  cat("\n=== AUDIT:", case, "===\n")
  audit_poisson_rs_joint(s[[case]], r[[case]], f$data, outcomes[[case]],
    if (case == "partial_higher") "higher" else "id", f$components)
}), names(s))
cat("DISJOINT_MAX_ABS_CROSS_COVARIANCE\n")
print(c(R = max(abs(r$disjoint$vcov[1:6, 7:12])),
  raw_stata = max(abs(s$disjoint$vcov[1:6, 7:12])),
  normalized_stata = max(abs(s$disjoint$vcov[1:6, 7:12]))*99/98,
  independent = max(abs(results$disjoint$independent_vcov[1:6, 7:12])),
  component_bread_stata_meat = max(abs(results$disjoint$component_bread_stata_meat_vcov[1:6, 7:12]))))
stopifnot(all(r$disjoint$vcov[1:6, 7:12] == 0),
  all(results$disjoint$independent_vcov[1:6, 7:12] == 0))
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) saveRDS(results, args[1])
cat("POISSON_RS_ALL_JOINT_AUDITS_COMPLETE=1\n")
