# Round 10 (November 7, 2026): linear models with two random slopes and with a
# level nested in the grouping variable (lmer; Stata mixed), with suest2 and
# each model alone, and the survey complementary log-log pair. The nonlinear
# random-slope models of this round are left out: R uses the Laplace
# approximation for them and Stata's comparison fits did not match it.
source("helpers.R")
source("rs_reference_v01.R")
source("compare_v10.R")
write_gaps(out, c("max_b", "max_se_rel", "max_V_rel"))
