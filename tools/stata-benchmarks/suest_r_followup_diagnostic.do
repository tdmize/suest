version 16
clear all
set more off
set linesize 255
capture log close _all
log using suest_r_followup_diagnostic.log, replace text name(followup)
di as result "FOLLOWUP_REVISION=1"
which suest2

* Isolate the native GEE score failure without changing any package code.
import delimited using suest_r_gee_crosslang_benchmark.csv, clear varnames(1)
xtset id time
foreach corr in independent exchangeable {
    di as txt "=== NATIVE GEE SCORE: `corr', FULL DATA WITH IF SAMPLE ==="
    quietly xtlogit binary x z if id >= 21 & time != 4, pa corr(`corr')
    estimates store diagnostic_logit
    generate byte diagnostic_sample = e(sample)
    capture noisily predict double diagnostic_score if diagnostic_sample, score
    di as result "NATIVE_SCORE_FULL_RC=" _rc
    capture drop diagnostic_score
    preserve
    keep if diagnostic_sample
    sort id time
    capture noisily predict double diagnostic_score, score
    di as result "NATIVE_SCORE_KEPT_RC=" _rc
    restore
    estimates restore diagnostic_logit
    quietly xtreg y x z if id <= 100, pa corr(`corr')
    estimates store diagnostic_gaussian
    capture noisily suest2 diagnostic_gaussian diagnostic_logit, cluster(higher)
    di as result "JOINT_GEE_RC=" _rc
    capture noisily suest2_cleanup, force
    drop diagnostic_sample
}

* Capture the native ML bread used by suest2 for the discrepant second model.
import delimited using suest_r_panel_fe_crosslang_benchmark.csv, clear varnames(1)
xtset id time
foreach y in y2 y2_partial {
    di as txt "=== NATIVE ML BREAD: `y' ==="
    quietly xtreg `y' x z i.time, mle
    matrix list e(b), format(%21.15g)
    matrix list e(V), format(%21.15g)
}

* This unbalanced RE case avoids the singular between design in plm.
di as txt "=== UNBALANCED RE WITHOUT TIME INDICATORS ==="
quietly xtreg y1_partial x z, re
estimates store re_simple1
quietly xtreg y2_partial x z, re
estimates store re_simple2
suest2 re_simple1 re_simple2, cluster(higher)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as result "FOLLOWUP_DIAGNOSTIC_COMPLETE=1"
log close followup
