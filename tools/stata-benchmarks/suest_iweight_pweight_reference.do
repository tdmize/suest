version 16.0
clear all
set more off
set linesize 255
set rmsg on

capture log close _all
log using "suest_iweight_pweight_reference.log", text replace name(iwbench)

display as text "=============================================================================="
display as text "SUEST IWEIGHT WORKAROUND / PWEIGHT REFERENCE BENCHMARK"
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
display as text "1. IDENTICAL-SAMPLE LINEAR MODELS"
display as text "=============================================================================="

estimates clear

regress y_linear c.x1##c.x1 x2 i.z [pweight=pw]
estimates store lin_pw_base
dump_estimation, label("PWEIGHT REFERENCE: linear base")

regress y_linear c.x1##c.x1 x2 i.z mediator [pweight=pw]
estimates store lin_pw_adjusted
dump_estimation, label("PWEIGHT REFERENCE: linear adjusted")

regress y_linear c.x1##c.x1 x2 i.z [iweight=pw]
estimates store lin_iw_base
dump_estimation, label("IWEIGHT INPUT TO SUEST: linear base")

regress y_linear c.x1##c.x1 x2 i.z mediator [iweight=pw]
estimates store lin_iw_adjusted
dump_estimation, label("IWEIGHT INPUT TO SUEST: linear adjusted")

run_suest lin_iw_base lin_iw_adjusted, ///
    label("SUEST IWEIGHT IDENTICAL-SAMPLE LINEAR")

display as text _newline "=============================================================================="
display as text "2. IDENTICAL-SAMPLE BINARY LOGIT MODELS"
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
dump_estimation, label("IWEIGHT INPUT TO SUEST: logit base")

logit y_binary i.z x1 x2 mediator [iweight=pw]
estimates store logit_iw_adjusted
dump_estimation, label("IWEIGHT INPUT TO SUEST: logit adjusted")

run_suest logit_iw_base logit_iw_adjusted, ///
    label("SUEST IWEIGHT IDENTICAL-SAMPLE LOGIT")

display as text _newline "=============================================================================="
display as text "3. PARTIAL OVERLAP WITH A COMMON WEIGHT"
display as text "=============================================================================="

estimates clear

regress y_linear c.x1##c.x1 x2 i.z [iweight=pw] if sample_a
estimates store partial_common_a
dump_estimation, label("IWEIGHT partial sample A, common pw")

regress y_linear c.x1##c.x1 x2 i.z mediator [iweight=pw] if sample_b
estimates store partial_common_b
dump_estimation, label("IWEIGHT partial sample B, common pw")

run_suest partial_common_a partial_common_b, ///
    label("SUEST IWEIGHT PARTIAL OVERLAP / COMMON WEIGHT")

display as text _newline "=============================================================================="
display as text "4. PARTIAL OVERLAP: SAME VALUES ON OVERLAP, DIFFERENT OFF OVERLAP"
display as text "=============================================================================="

count if sample_a & sample_b & pw != pw_partial_b
display as text "Different weights on overlap = " r(N)
count if !sample_a & sample_b & pw != pw_partial_b
display as text "Different weights on model-B-only observations = " r(N)

estimates clear

regress y_linear x1 x2 [iweight=pw] if sample_a
estimates store partial_diff_a
dump_estimation, label("IWEIGHT partial sample A, pw")

regress y_linear x1 x2 mediator [iweight=pw_partial_b] if sample_b
estimates store partial_diff_b
dump_estimation, label("IWEIGHT partial sample B, pw_partial_b")

run_suest partial_diff_a partial_diff_b, ///
    label("SUEST IWEIGHT SAME ON OVERLAP / DIFFERENT OFF OVERLAP")

display as text _newline "=============================================================================="
display as text "5. DISJOINT SAMPLES WITH A COMMON WEIGHT"
display as text "=============================================================================="

estimates clear

regress y_linear c.x1##c.x1 x2 i.z [iweight=pw] if sample_left
estimates store disjoint_common_left
dump_estimation, label("IWEIGHT disjoint left, common pw")

regress y_linear c.x1##c.x1 x2 i.z mediator [iweight=pw] if sample_right
estimates store disjoint_common_right
dump_estimation, label("IWEIGHT disjoint right, common pw")

run_suest disjoint_common_left disjoint_common_right, ///
    label("SUEST IWEIGHT DISJOINT / COMMON WEIGHT")

display as text _newline "=============================================================================="
display as text "6. DISJOINT SAMPLES WITH ENTIRELY DIFFERENT WEIGHTS"
display as text "=============================================================================="

estimates clear

regress y_linear x1 x2 [iweight=pw_left] if sample_left
estimates store disjoint_diff_left
dump_estimation, label("IWEIGHT disjoint left, pw_left")

regress y_linear x1 x2 mediator [iweight=pw_right] if sample_right
estimates store disjoint_diff_right
dump_estimation, label("IWEIGHT disjoint right, pw_right")

run_suest disjoint_diff_left disjoint_diff_right, ///
    label("SUEST IWEIGHT DISJOINT / DIFFERENT WEIGHTS")

display as text _newline "=============================================================================="
display as text "7. SAME SAMPLE: IDENTICAL VALUES UNDER DIFFERENT VARIABLE NAMES"
display as text "=============================================================================="

estimates clear

regress y_linear x1 x2 [iweight=pw]
estimates store same_values_a
regress y_linear x1 x2 mediator [iweight=pw_copy]
estimates store same_values_b

run_suest same_values_a same_values_b, ///
    label("SUEST IWEIGHT IDENTICAL VALUES / DIFFERENT NAMES")

display as text _newline "=============================================================================="
display as text "8. SAME SAMPLE: WEIGHTS DIFFER BY A CONSTANT FACTOR"
display as text "=============================================================================="

estimates clear

regress y_linear x1 x2 [iweight=pw]
estimates store scale_a
dump_estimation, label("IWEIGHT model using pw")

regress y_linear x1 x2 mediator [iweight=pw10]
estimates store scale_b
dump_estimation, label("IWEIGHT model using 10*pw")

run_suest scale_a scale_b, ///
    label("SUEST IWEIGHT DIFFERENT WEIGHT SCALE")

display as text _newline "=============================================================================="
display as text "9. COMPLETED"
display as text "=============================================================================="
display as text "Upload:"
display as result "  suest_iweight_pweight_reference.log"
display as text "Completed: " c(current_date) " " c(current_time)

log close iwbench
