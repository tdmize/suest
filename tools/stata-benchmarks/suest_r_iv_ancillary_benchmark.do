version 16
clear all
set more off

capture log close _all
log using suest_r_iv_ancillary_benchmark.log, replace text name(ivancillary)

di as txt "SUEST R IV AND ANCILLARY BENCHMARK"
di as txt "Stata version: `c(stata_version)'"
which suest2
import delimited using suest_r_iv_ancillary_benchmark.csv, clear varnames(1)

di as txt "=== CASE 1: IDENTICAL-SAMPLE 2SLS ==="
quietly ivregress 2sls y x (endogenous = z1)
estimates store iv_base
quietly ivregress 2sls y x (endogenous = z1 z2)
estimates store iv_adjusted
suest2 iv_base iv_adjusted, vce(robust)
di as txt "N=" e(N)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily margins, dydx(x) predict(model(iv_base) xb)
capture noisily margins, dydx(x) predict(model(iv_adjusted) xb)
capture noisily suest2_cleanup, force

di as txt "=== CASE 2: PARTIALLY OVERLAPPING 2SLS ==="
quietly ivregress 2sls y x (endogenous = z1 z2) if id <= 900
estimates store iv_left
quietly ivregress 2sls y x (endogenous = z1 z2) if id > 300
estimates store iv_right
suest2 iv_left iv_right, vce(robust)
di as txt "N=" e(N)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as txt "=== CASE 3: REGRESS ANCILLARY PARAMETER ==="
quietly regress y x
estimates store lm_base
quietly regress y x z1
estimates store lm_adjusted
suest2 lm_base lm_adjusted
di as txt "N=" e(N)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as txt "=== CASE 4: THREE-MODEL SYNTAX ==="
quietly regress y x
estimates store lm_one
quietly regress y x z1
estimates store lm_two
quietly regress y x z1 z2
estimates store lm_three
capture noisily suest2 lm_one lm_two lm_three
local three_model_rc = _rc
di as txt "THREE_MODEL_RC=`three_model_rc'"
if `three_model_rc' == 0 {
    di as txt "N=" e(N)
    matrix list e(b), format(%21.15g)
    matrix list e(V), format(%21.15g)
}
capture noisily suest2_cleanup, force

di as txt "=== CASE 5: HETEROSKEDASTIC PROBIT ==="
quietly hetprobit y_het x, het(z1)
estimates store het_base
quietly hetprobit y_het x z2, het(z1 z2)
estimates store het_adjusted
suest2 het_base het_adjusted, vce(robust)
di as txt "N=" e(N)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily margins, dydx(x) predict(model(het_base) pr)
capture noisily margins, dydx(x) predict(model(het_adjusted) pr)
capture noisily suest2_cleanup, force

di as result "IV_ANCILLARY_BENCHMARK_COMPLETE=1"
log close ivancillary
