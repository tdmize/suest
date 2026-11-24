# Round 18 (November 8, 2026): Gompertz and generalized gamma survival models
# (flexsurv; Stata streg) with ancillary covariates and delayed entry, and
# bivariate probit with different predictors in each equation (GJRM; Stata
# biprobit), with suest2 and each model alone, and average predictions.
source("helpers.R")
source("make_v18.R")
source("compare_v18.R")
write_gaps(exact, c("max_b", "max_se_rel", "max_V_rel"))
write_gaps(other[1:3, ], c("abs_diff", "se_rel"), check = "check", append = TRUE)
