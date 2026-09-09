version 16.0
clear all
set more off
set linesize 255
set rmsg on
set seed 3762026

capture log close _all
log using "suest_pweight_poisson_benchmark.log", text replace name(poissonpw)

display as text "=============================================================================="
display as text "SUEST PWEIGHT POISSON BENCHMARK"
display as text "Started: " c(current_date) " " c(current_time)
display as text "Stata version: " c(stata_version)
display as text "=============================================================================="

use "suest_pweight_benchmark_data.dta", clear

* Generate a reproducible count outcome and save the expanded benchmark data.
generate double mu_count = exp(0.35 + 0.28*x1 - 0.18*x2 + 0.32*z + ///
    0.22*mediator + 0.08*x1*z)
generate long y_count = rpoisson(mu_count)
drop mu_count
compress
save "suest_pweight_poisson_benchmark_data.dta", replace

capture program drop dump_estimation
program define dump_estimation
    syntax, Label(string)

    display as text _newline "------------------------------------------------------------------------------"
    display as text "`label'"
    display as text "------------------------------------------------------------------------------"
    estimates replay
    display as text "e(cmd)       = `e(cmd)'"
    display as text "e(cmdline)   = `e(cmdline)'"
    display as text "e(N)         = " %21.15g e(N)
    capture display as text "e(N_clust)   = " %21.15g e(N_clust)
    display as text "e(vce)       = `e(vce)'"
    display as text "e(vcetype)   = `e(vcetype)'"
    display as text "e(wtype)     = `e(wtype)'"
    display as text "e(wexp)      = `e(wexp)'"

    tempname b V
    matrix `b' = e(b)
    matrix `V' = e(V)

    display as text _newline "COEFFICIENT VECTOR e(b)"
    matrix list `b', format(%21.15g)

    display as text _newline "VARIANCE-COVARIANCE MATRIX e(V)"
    matrix list `V', format(%21.15g)
end

capture program drop run_suest
program define run_suest
    syntax namelist(min=2 max=2), Label(string)

    capture noisily suest `namelist'
    local rc = _rc
    display as result "`label' RETURN CODE = `rc'"

    if `rc' == 0 {
        dump_estimation, label("`label'")
        capture noisily suest, coeflegend
    }
end

display as text _newline "=============================================================================="
display as text "0. DATA CHECKS"
display as text "=============================================================================="

describe y_count
summarize y_count x1 x2 z mediator pw
tabulate y_count if y_count <= 10
count if y_count == 0
display as text "Zero-count observations = " r(N)

display as text _newline "=============================================================================="
display as text "1. IDENTICAL-SAMPLE POISSON MODELS"
display as text "=============================================================================="

estimates clear

poisson y_count i.z x1 x2 [pweight=pw]
estimates store poisson_pw_base
dump_estimation, label("PWEIGHT REFERENCE: Poisson base")

poisson y_count i.z x1 x2 mediator [pweight=pw]
estimates store poisson_pw_adjusted
dump_estimation, label("PWEIGHT REFERENCE: Poisson adjusted")

poisson y_count i.z x1 x2 [iweight=pw]
estimates store poisson_iw_base

poisson y_count i.z x1 x2 mediator [iweight=pw]
estimates store poisson_iw_adjusted

run_suest poisson_iw_base poisson_iw_adjusted, ///
    label("SUEST PWEIGHT IDENTICAL-SAMPLE POISSON")

display as text _newline "=============================================================================="
display as text "2. PARTIALLY OVERLAPPING POISSON SAMPLES"
display as text "=============================================================================="

estimates clear

poisson y_count i.z x1 x2 [iweight=pw] if sample_a
estimates store poisson_partial_a

poisson y_count i.z x1 x2 mediator [iweight=pw] if sample_b
estimates store poisson_partial_b

run_suest poisson_partial_a poisson_partial_b, ///
    label("SUEST PWEIGHT PARTIAL-OVERLAP POISSON")

display as text _newline "=============================================================================="
display as text "3. DISJOINT POISSON SAMPLES"
display as text "=============================================================================="

estimates clear

poisson y_count i.z x1 x2 [iweight=pw] if sample_left
estimates store poisson_disjoint_left

poisson y_count i.z x1 x2 mediator [iweight=pw] if sample_right
estimates store poisson_disjoint_right

run_suest poisson_disjoint_left poisson_disjoint_right, ///
    label("SUEST PWEIGHT DISJOINT POISSON")

display as text _newline "=============================================================================="
display as text "4. POISSON WEIGHT-RULE CHECKS"
display as text "=============================================================================="

display as text _newline "4A. IDENTICAL VALUES UNDER DIFFERENT VARIABLE NAMES"
estimates clear
poisson y_count x1 x2 [iweight=pw]
estimates store poisson_same_a
poisson y_count x1 x2 mediator [iweight=pw_copy]
estimates store poisson_same_b
run_suest poisson_same_a poisson_same_b, ///
    label("SUEST PWEIGHT POISSON IDENTICAL VALUES / DIFFERENT NAMES")

display as text _newline "4B. SAME VALUES ON OVERLAP, DIFFERENT OFF OVERLAP"
estimates clear
poisson y_count x1 x2 [iweight=pw] if sample_a
estimates store poisson_overlap_a
poisson y_count x1 x2 mediator [iweight=pw_partial_b] if sample_b
estimates store poisson_overlap_b
run_suest poisson_overlap_a poisson_overlap_b, ///
    label("SUEST PWEIGHT POISSON SAME ON OVERLAP / DIFFERENT OFF OVERLAP")

display as text _newline "4C. DISJOINT SAMPLES WITH DIFFERENT WEIGHTS"
estimates clear
poisson y_count x1 x2 [iweight=pw_left] if sample_left
estimates store poisson_disjoint_diff_a
poisson y_count x1 x2 mediator [iweight=pw_right] if sample_right
estimates store poisson_disjoint_diff_b
run_suest poisson_disjoint_diff_a poisson_disjoint_diff_b, ///
    label("SUEST PWEIGHT POISSON DISJOINT / DIFFERENT WEIGHTS")

display as text _newline "4D. SAME SAMPLE, WEIGHTS DIFFER BY A CONSTANT FACTOR"
estimates clear
poisson y_count x1 x2 [iweight=pw]
estimates store poisson_scale_a
poisson y_count x1 x2 mediator [iweight=pw10]
estimates store poisson_scale_b
run_suest poisson_scale_a poisson_scale_b, ///
    label("SUEST PWEIGHT POISSON DIFFERENT WEIGHT SCALE")

display as text _newline "=============================================================================="
display as text "5. COMPLETED"
display as text "=============================================================================="
display as text "Upload both files:"
display as result "  suest_pweight_poisson_benchmark.log"
display as result "  suest_pweight_poisson_benchmark_data.dta"
display as text "Completed: " c(current_date) " " c(current_time)

log close poissonpw
