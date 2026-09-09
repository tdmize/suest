version 16.0
clear all
set more off
set linesize 255
set rmsg on

capture log close _all
log using "suest_pweight_binary_benchmark.log", text replace name(binarypw)

display as text "=============================================================================="
display as text "SUEST PWEIGHT BINARY-MODEL BENCHMARK"
display as text "Started: " c(current_date) " " c(current_time)
display as text "Stata version: " c(stata_version)
display as text "=============================================================================="

use "suest_pweight_benchmark_data.dta", clear

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
display as text "1. IDENTICAL-SAMPLE LOGIT MODELS"
display as text "=============================================================================="

estimates clear

logit y_binary i.z x1 x2 [pweight=pw]
estimates store logit_pw_base
dump_estimation, label("PWEIGHT REFERENCE: logit base")

logit y_binary i.z x1 x2 mediator [pweight=pw]
estimates store logit_pw_adjusted
dump_estimation, label("PWEIGHT REFERENCE: logit adjusted")

logit y_binary i.z x1 x2 [iweight=pw]
estimates store logit_iw_base

logit y_binary i.z x1 x2 mediator [iweight=pw]
estimates store logit_iw_adjusted

run_suest logit_iw_base logit_iw_adjusted, ///
    label("SUEST PWEIGHT IDENTICAL-SAMPLE LOGIT")

display as text _newline "=============================================================================="
display as text "2. IDENTICAL-SAMPLE PROBIT MODELS"
display as text "=============================================================================="

estimates clear

probit y_binary i.z x1 x2 [pweight=pw]
estimates store probit_pw_base
dump_estimation, label("PWEIGHT REFERENCE: probit base")

probit y_binary i.z x1 x2 mediator [pweight=pw]
estimates store probit_pw_adjusted
dump_estimation, label("PWEIGHT REFERENCE: probit adjusted")

probit y_binary i.z x1 x2 [iweight=pw]
estimates store probit_iw_base

probit y_binary i.z x1 x2 mediator [iweight=pw]
estimates store probit_iw_adjusted

run_suest probit_iw_base probit_iw_adjusted, ///
    label("SUEST PWEIGHT IDENTICAL-SAMPLE PROBIT")

display as text _newline "=============================================================================="
display as text "3. IDENTICAL-SAMPLE LOGIT-PROBIT PAIR"
display as text "=============================================================================="

estimates clear

logit y_binary i.z x1 x2 mediator [iweight=pw]
estimates store mixed_logit

probit y_binary i.z x1 x2 mediator [iweight=pw]
estimates store mixed_probit

run_suest mixed_logit mixed_probit, ///
    label("SUEST PWEIGHT IDENTICAL-SAMPLE LOGIT-PROBIT")

display as text _newline "=============================================================================="
display as text "4. PARTIALLY OVERLAPPING LOGIT SAMPLES"
display as text "=============================================================================="

estimates clear

logit y_binary i.z x1 x2 [iweight=pw] if sample_a
estimates store logit_partial_a

logit y_binary i.z x1 x2 mediator [iweight=pw] if sample_b
estimates store logit_partial_b

run_suest logit_partial_a logit_partial_b, ///
    label("SUEST PWEIGHT PARTIAL-OVERLAP LOGIT")

display as text _newline "=============================================================================="
display as text "5. DISJOINT PROBIT SAMPLES"
display as text "=============================================================================="

estimates clear

probit y_binary i.z x1 x2 [iweight=pw] if sample_left
estimates store probit_disjoint_left

probit y_binary i.z x1 x2 mediator [iweight=pw] if sample_right
estimates store probit_disjoint_right

run_suest probit_disjoint_left probit_disjoint_right, ///
    label("SUEST PWEIGHT DISJOINT PROBIT")

display as text _newline "=============================================================================="
display as text "6. IDENTICAL-SAMPLE LINEAR-LOGIT PAIR"
display as text "=============================================================================="

estimates clear

regress y_linear x1 x2 i.z [iweight=pw]
estimates store mixed_linear_logit_lm

logit y_binary i.z x1 x2 [iweight=pw]
estimates store mixed_linear_logit_bin

run_suest mixed_linear_logit_lm mixed_linear_logit_bin, ///
    label("SUEST PWEIGHT IDENTICAL-SAMPLE LINEAR-LOGIT")

display as text _newline "=============================================================================="
display as text "7. IDENTICAL-SAMPLE LINEAR-PROBIT PAIR"
display as text "=============================================================================="

estimates clear

regress y_linear x1 x2 i.z [iweight=pw]
estimates store mixed_linear_probit_lm

probit y_binary i.z x1 x2 [iweight=pw]
estimates store mixed_linear_probit_bin

run_suest mixed_linear_probit_lm mixed_linear_probit_bin, ///
    label("SUEST PWEIGHT IDENTICAL-SAMPLE LINEAR-PROBIT")

display as text _newline "=============================================================================="
display as text "8. WEIGHT-RULE CHECKS FOR BINARY MODELS"
display as text "=============================================================================="

display as text _newline "8A. SAME VALUES UNDER DIFFERENT VARIABLE NAMES"
estimates clear
logit y_binary x1 x2 [iweight=pw]
estimates store binary_same_a
logit y_binary x1 x2 mediator [iweight=pw_copy]
estimates store binary_same_b
run_suest binary_same_a binary_same_b, ///
    label("SUEST PWEIGHT BINARY IDENTICAL VALUES / DIFFERENT NAMES")

display as text _newline "8B. SAME ON OVERLAP, DIFFERENT OFF OVERLAP"
estimates clear
logit y_binary x1 x2 [iweight=pw] if sample_a
estimates store binary_overlap_a
logit y_binary x1 x2 mediator [iweight=pw_partial_b] if sample_b
estimates store binary_overlap_b
run_suest binary_overlap_a binary_overlap_b, ///
    label("SUEST PWEIGHT BINARY SAME ON OVERLAP / DIFFERENT OFF OVERLAP")

display as text _newline "8C. DISJOINT SAMPLES WITH DIFFERENT WEIGHTS"
estimates clear
probit y_binary x1 x2 [iweight=pw_left] if sample_left
estimates store binary_disjoint_a
probit y_binary x1 x2 mediator [iweight=pw_right] if sample_right
estimates store binary_disjoint_b
run_suest binary_disjoint_a binary_disjoint_b, ///
    label("SUEST PWEIGHT BINARY DISJOINT / DIFFERENT WEIGHTS")

display as text _newline "=============================================================================="
display as text "9. COMPLETED"
display as text "=============================================================================="
display as text "Upload:"
display as result "  suest_pweight_binary_benchmark.log"
display as text "Completed: " c(current_date) " " c(current_time)

log close binarypw
