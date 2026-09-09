version 16
clear all
set more off

capture log close _all
log using suest_r_pweight_lnvar_scale_diagnostic.log, replace text name(pwlnvar)

di as result "DIAGNOSTIC_FILE_REVISION=1"
di as txt "SUEST2 PWEIGHT LINEAR ANCILLARY SCALE DIAGNOSTIC"
di as txt "Stata version: `c(stata_version)'"
which suest2

import delimited using suest_r_cluster_crosslang_benchmark.csv, clear varnames(1)
generate double pw10 = 10*pw

tempname B1 V1 B10 V10
tempvar resid1 wrss1 resid10 wrss10

di as txt "=== ORIGINAL SCALE ==="
quietly regress y x z [pweight=pw]
scalar original_lnvar1 = ln(e(rmse)^2)
scalar k1 = e(rank)
predict double `resid1' if e(sample), residuals
generate double `wrss1' = pw*`resid1'^2 if e(sample)
quietly summarize pw if e(sample), meanonly
scalar sumw1 = r(sum)
quietly summarize `wrss1' if e(sample), meanonly
scalar formula_lnvar1 = ln(r(sum)/(sumw1-k1))
estimates store scale1_lm

quietly logit b x z [pweight=pw]
estimates store scale1_logit
suest2 scale1_lm scale1_logit, cluster(cluster4)
matrix `B1' = e(b)
matrix `V1' = e(V)
scalar suest_lnvar1 = el(`B1',1,4)
capture noisily suest2_cleanup, force

di as txt "=== TEN TIMES SCALE ==="
quietly regress y x z [pweight=pw10]
scalar original_lnvar10 = ln(e(rmse)^2)
scalar k10 = e(rank)
predict double `resid10' if e(sample), residuals
generate double `wrss10' = pw10*`resid10'^2 if e(sample)
quietly summarize pw10 if e(sample), meanonly
scalar sumw10 = r(sum)
quietly summarize `wrss10' if e(sample), meanonly
scalar formula_lnvar10 = ln(r(sum)/(sumw10-k10))
estimates store scale10_lm

quietly logit b x z [pweight=pw10]
estimates store scale10_logit
suest2 scale10_lm scale10_logit, cluster(cluster4)
matrix `B10' = e(b)
matrix `V10' = e(V)
scalar suest_lnvar10 = el(`B10',1,4)
capture noisily suest2_cleanup, force

di as txt "=== SCALE COMPARISONS ==="
di as result "SUM_WEIGHTS_1=" %21.15g sumw1
di as result "SUM_WEIGHTS_10=" %21.15g sumw10
di as result "ORIGINAL_PWEIGHT_LNVAR_1=" %21.15g original_lnvar1
di as result "ORIGINAL_PWEIGHT_LNVAR_10=" %21.15g original_lnvar10
di as result "ORIGINAL_PWEIGHT_LNVAR_DIFF=" %21.15g (original_lnvar10-original_lnvar1)
di as result "IWEIGHT_FORMULA_LNVAR_1=" %21.15g formula_lnvar1
di as result "IWEIGHT_FORMULA_LNVAR_10=" %21.15g formula_lnvar10
di as result "SUEST2_LNVAR_1=" %21.15g suest_lnvar1
di as result "SUEST2_LNVAR_10=" %21.15g suest_lnvar10
di as result "SUEST2_LNVAR_DIFF=" %21.15g (suest_lnvar10-suest_lnvar1)
di as result "SUEST2_MINUS_FORMULA_1=" %21.15g (suest_lnvar1-formula_lnvar1)
di as result "SUEST2_MINUS_FORMULA_10=" %21.15g (suest_lnvar10-formula_lnvar10)
di as result "LM_X_COEF_DIFF=" %21.15g (el(`B10',1,1)-el(`B1',1,1))
di as result "LM_X_VARIANCE_DIFF=" %21.15g (el(`V10',1,1)-el(`V1',1,1))
di as result "LNVAR_VARIANCE_DIFF=" %21.15g (el(`V10',4,4)-el(`V1',4,4))
di as result "LOGIT_X_VARIANCE_DIFF=" %21.15g (el(`V10',5,5)-el(`V1',5,5))

di as result "PWEIGHT_LNVAR_SCALE_DIAGNOSTIC_COMPLETE=1"
log close pwlnvar
