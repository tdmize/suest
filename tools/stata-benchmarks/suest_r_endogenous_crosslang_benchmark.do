version 16
clear all
set more off

capture log close _all
log using suest_r_endogenous_crosslang_benchmark.log, replace text name(endogenous)

di as result "BENCHMARK_FILE_REVISION=1"
di as txt "SUEST R BIPROBIT/IVPROBIT CROSS-LANGUAGE BENCHMARK"
di as txt "Stata version: `c(stata_version)'"
which suest2

import delimited using suest_r_endogenous_crosslang_benchmark.csv, clear varnames(1)

di as txt "=== CASE 1: BIPROBIT, SAME REGRESSORS WITHIN EACH MODEL ==="
quietly biprobit (binary1 = x) (binary2 = x)
estimates store biprobit_base
quietly biprobit (binary1 = x z) (binary2 = x z)
estimates store biprobit_adjusted
suest2 biprobit_base biprobit_adjusted
di as txt "N=" e(N)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as txt "=== CASE 2: IVPROBIT ML, IDENTICAL SAMPLE ==="
quietly ivprobit binary_iv x (endogenous = z1), ml
estimates store ivprobit_base
quietly ivprobit binary_iv x (endogenous = z1 z2), ml
estimates store ivprobit_adjusted
suest2 ivprobit_base ivprobit_adjusted
di as txt "N=" e(N)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as result "ENDOGENOUS_CROSSLANG_BENCHMARK_COMPLETE=1"
log close endogenous
