version 16
clear all
set more off

capture log close _all
log using suest_r_pweight_extended_benchmark.log, replace text name(pwextended)

di as txt "SUEST R EXTENDED PWEIGHT BENCHMARK"
di as txt "Stata version: `c(stata_version)'"
which suest2
import delimited using suest_r_pweight_extended_benchmark.csv, clear varnames(1)

di as txt "=== CASE 1: PWEIGHT NEGATIVE BINOMIAL ==="
quietly nbreg y_nb x1 x2 i.z [pweight=pw]
estimates store pw_nb_base
quietly nbreg y_nb x1 x2 i.z mediator [pweight=pw]
estimates store pw_nb_adjusted
suest2 pw_nb_base pw_nb_adjusted
di as txt "N=" e(N)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as txt "=== CASE 2: THREE-MODEL PWEIGHT CATEGORICAL SYSTEM ==="
quietly ologit y_ord x1 x2 i.z [pweight=pw]
estimates store pw_ologit
quietly oprobit y_ord x1 x2 i.z [pweight=pw]
estimates store pw_oprobit
quietly mlogit y_ord x1 x2 i.z [pweight=pw], baseoutcome(1)
estimates store pw_mlogit
suest2 pw_ologit pw_oprobit pw_mlogit
di as txt "N=" e(N)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as result "PWEIGHT_EXTENDED_BENCHMARK_COMPLETE=1"
log close pwextended
