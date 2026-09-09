version 16.0
clear all
set more off
set linesize 255
set rmsg on
set seed 3762026

capture log close _all
log using "suest_pweight_benchmark.log", text replace name(pwbench)

display as text "=============================================================================="
display as text "SUEST PWEIGHT BENCHMARK SUITE"
display as text "Started: " c(current_date) " " c(current_time)
display as text "Stata version: " c(stata_version)
display as text "Operating system: " c(os)
display as text "=============================================================================="

display as text _newline "Installed commands"
capture noisily which mecompare
display as text "which mecompare return code = " _rc
capture noisily which melincom
display as text "which melincom return code = " _rc

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

capture program drop dump_postestimation
program define dump_postestimation
    syntax, Label(string)

    display as text _newline "------------------------------------------------------------------------------"
    display as text "`label'"
    display as text "------------------------------------------------------------------------------"

    capture noisily ereturn list
    capture noisily return list

    capture matrix __post_b = e(b)
    if !_rc {
        display as text _newline "POSTESTIMATION e(b)"
        matrix list __post_b, format(%21.15g)
    }

    capture matrix __post_V = e(V)
    if !_rc {
        display as text _newline "POSTESTIMATION e(V)"
        matrix list __post_V, format(%21.15g)
    }

    capture matrix drop __post_b __post_V
end

display as text _newline "=============================================================================="
display as text "0. GENERATE AND SAVE REPRODUCIBLE BENCHMARK DATA"
display as text "=============================================================================="

set obs 2400
generate long id = _n
generate double x1 = rnormal()
generate double x2 = rnormal()
generate byte z = runiform() > 0.48
generate double mediator = 0.55*x1 - 0.30*x2 + 0.35*z + rnormal()

generate double pw = exp(0.30*x2 + 0.20*z + 0.15*rnormal())
replace pw = max(pw, 0.05)
generate double pw_copy = pw
generate double pw10 = 10*pw

generate byte sample_a = id <= 1800
generate byte sample_b = id >= 601
generate byte sample_left = id <= 1200
generate byte sample_right = id > 1200

* Same weight on the partial-overlap region (IDs 601-1800), different only
* on observations unique to model B (IDs 1801-2400).
generate double pw_partial_b = pw
replace pw_partial_b = pw*(1.25 + 0.10*abs(x1)) if id > 1800

* Entirely different weights for two disjoint samples.
generate double pw_left = exp(0.25*x1 - 0.15*x2 + 0.10*rnormal())
generate double pw_right = exp(-0.20*x1 + 0.30*x2 + 0.10*rnormal())

generate double y_linear = 1.10 + 0.75*x1 - 0.18*x1^2 - 0.40*x2 + ///
    0.55*z + 0.65*mediator + (0.70 + 0.25*abs(x1))*rnormal()

generate double xb = -0.35 + 0.85*x1 - 0.30*x2 + 0.65*z + ///
    0.45*mediator + 0.25*x1*z
generate byte y_binary = runiform() < invlogit(xb)
drop xb

order id y_linear y_binary x1 x2 z mediator pw pw_copy pw10 ///
    pw_partial_b pw_left pw_right sample_a sample_b sample_left sample_right
compress
save "suest_pweight_benchmark_data.dta", replace

describe
summarize y_linear y_binary x1 x2 z mediator pw pw_copy pw10 ///
    pw_partial_b pw_left pw_right
tabulate sample_a sample_b
tabulate sample_left sample_right

display as text _newline "Partial-overlap weight agreement check"
count if sample_a & sample_b
display as text "overlap N = " r(N)
count if sample_a & sample_b & pw != pw_partial_b
display as text "different weights on overlap = " r(N)
count if !sample_a & sample_b & pw != pw_partial_b
display as text "different weights on model-B-only observations = " r(N)

display as text _newline "=============================================================================="
display as text "1. IDENTICAL-SAMPLE WEIGHTED LINEAR MODELS"
display as text "=============================================================================="

estimates clear

regress y_linear c.x1##c.x1 x2 i.z [pweight=pw]
estimates store lin_base
dump_estimation, label("Linear base model used by suest")

regress y_linear c.x1##c.x1 x2 i.z mediator [pweight=pw]
estimates store lin_adjusted
dump_estimation, label("Linear adjusted model used by suest")

capture noisily suest lin_base lin_adjusted
local rc = _rc
display as result "SUEST IDENTICAL-SAMPLE LINEAR RETURN CODE = `rc'"
if `rc' == 0 {
    dump_estimation, label("SUEST: identical-sample weighted linear models")
}

regress y_linear c.x1##c.x1 x2 i.z [pweight=pw], vce(robust)
estimates store lin_base_robust
dump_estimation, label("Linear base model with explicit vce(robust)")

regress y_linear c.x1##c.x1 x2 i.z mediator [pweight=pw], vce(robust)
estimates store lin_adjusted_robust
dump_estimation, label("Linear adjusted model with explicit vce(robust)")

capture noisily mecompare x1, models(lin_base_robust lin_adjusted_robust) amount(1)
local rc = _rc
display as result "MECOMPARE IDENTICAL-SAMPLE LINEAR RETURN CODE = `rc'"
if `rc' == 0 {
    dump_postestimation, label("MECOMPARE: identical-sample weighted linear models, x1 + 1")
}

display as text _newline "=============================================================================="
display as text "2. IDENTICAL-SAMPLE WEIGHTED BINARY LOGIT MODELS"
display as text "=============================================================================="

estimates clear

logit y_binary i.z x1 x2 [pweight=pw]
estimates store logit_base
dump_estimation, label("Logit base model used by suest")

logit y_binary i.z x1 x2 mediator [pweight=pw]
estimates store logit_adjusted
dump_estimation, label("Logit adjusted model used by suest")

capture noisily suest logit_base logit_adjusted
local rc = _rc
display as result "SUEST IDENTICAL-SAMPLE LOGIT RETURN CODE = `rc'"
if `rc' == 0 {
    dump_estimation, label("SUEST: identical-sample weighted logit models")
}

logit y_binary i.z x1 x2 [pweight=pw], vce(robust)
estimates store logit_base_robust
dump_estimation, label("Logit base model with explicit vce(robust)")

logit y_binary i.z x1 x2 mediator [pweight=pw], vce(robust)
estimates store logit_adjusted_robust
dump_estimation, label("Logit adjusted model with explicit vce(robust)")

capture noisily mecompare i.z, models(logit_base_robust logit_adjusted_robust)
local rc = _rc
display as result "MECOMPARE IDENTICAL-SAMPLE LOGIT RETURN CODE = `rc'"
if `rc' == 0 {
    dump_postestimation, label("MECOMPARE: identical-sample weighted logit models, z 1 - 0")
}

display as text _newline "=============================================================================="
display as text "3. PARTIALLY OVERLAPPING SAMPLES WITH ONE COMMON PWEIGHT"
display as text "=============================================================================="

estimates clear

regress y_linear c.x1##c.x1 x2 i.z [pweight=pw] if sample_a
estimates store partial_a
dump_estimation, label("Partial sample A weighted linear model")

regress y_linear c.x1##c.x1 x2 i.z mediator [pweight=pw] if sample_b
estimates store partial_b
dump_estimation, label("Partial sample B weighted linear model")

capture noisily suest partial_a partial_b
local rc = _rc
display as result "SUEST PARTIAL-OVERLAP COMMON-WEIGHT RETURN CODE = `rc'"
if `rc' == 0 {
    dump_estimation, label("SUEST: partial overlap with common pweight")
}

regress y_linear c.x1##c.x1 x2 i.z [pweight=pw] if sample_a, vce(robust)
estimates store partial_a_robust
regress y_linear c.x1##c.x1 x2 i.z mediator [pweight=pw] if sample_b, vce(robust)
estimates store partial_b_robust

capture noisily mecompare x1, models(partial_a_robust partial_b_robust) amount(1)
local rc = _rc
display as result "MECOMPARE PARTIAL-OVERLAP COMMON-WEIGHT RETURN CODE = `rc'"
if `rc' == 0 {
    dump_postestimation, label("MECOMPARE: partial overlap with common pweight, x1 + 1")
}

display as text _newline "=============================================================================="
display as text "4. DISJOINT SAMPLES WITH ONE COMMON PWEIGHT"
display as text "=============================================================================="

estimates clear

regress y_linear c.x1##c.x1 x2 i.z [pweight=pw] if sample_left
estimates store disjoint_left
dump_estimation, label("Disjoint left sample weighted linear model")

regress y_linear c.x1##c.x1 x2 i.z mediator [pweight=pw] if sample_right
estimates store disjoint_right
dump_estimation, label("Disjoint right sample weighted linear model")

capture noisily suest disjoint_left disjoint_right
local rc = _rc
display as result "SUEST DISJOINT COMMON-WEIGHT RETURN CODE = `rc'"
if `rc' == 0 {
    dump_estimation, label("SUEST: disjoint samples with common pweight")
}

regress y_linear c.x1##c.x1 x2 i.z [pweight=pw] if sample_left, vce(robust)
estimates store disjoint_left_robust
regress y_linear c.x1##c.x1 x2 i.z mediator [pweight=pw] if sample_right, vce(robust)
estimates store disjoint_right_robust

capture noisily mecompare x1, models(disjoint_left_robust disjoint_right_robust) amount(1)
local rc = _rc
display as result "MECOMPARE DISJOINT COMMON-WEIGHT RETURN CODE = `rc'"
if `rc' == 0 {
    dump_postestimation, label("MECOMPARE: disjoint samples with common pweight, x1 + 1")
}

display as text _newline "=============================================================================="
display as text "5. WEIGHT-COMPATIBILITY EXPERIMENTS"
display as text "=============================================================================="

display as text _newline "5A. SAME SAMPLE; DIFFERENT VARIABLE NAMES BUT IDENTICAL VALUES"
estimates clear
regress y_linear x1 x2 [pweight=pw]
estimates store same_name_a
regress y_linear x1 x2 mediator [pweight=pw_copy]
estimates store same_name_b
capture noisily suest same_name_a same_name_b
local rc = _rc
display as result "SUEST IDENTICAL VALUES / DIFFERENT NAMES RETURN CODE = `rc'"
if `rc' == 0 {
    dump_estimation, label("SUEST: identical weight values stored under different variable names")
}

display as text _newline "5B. PARTIAL OVERLAP; SAME WEIGHTS ON OVERLAP, DIFFERENT OFF OVERLAP"
estimates clear
regress y_linear x1 x2 [pweight=pw] if sample_a
estimates store overlap_weight_a
regress y_linear x1 x2 mediator [pweight=pw_partial_b] if sample_b
estimates store overlap_weight_b
capture noisily suest overlap_weight_a overlap_weight_b
local rc = _rc
display as result "SUEST SAME ON OVERLAP / DIFFERENT OFF OVERLAP RETURN CODE = `rc'"
if `rc' == 0 {
    dump_estimation, label("SUEST: same weights on overlap but different on unique observations")
}

display as text _newline "5C. DISJOINT SAMPLES; ENTIRELY DIFFERENT WEIGHT VARIABLES AND VALUES"
estimates clear
regress y_linear x1 x2 [pweight=pw_left] if sample_left
estimates store different_left
regress y_linear x1 x2 mediator [pweight=pw_right] if sample_right
estimates store different_right
capture noisily suest different_left different_right
local rc = _rc
display as result "SUEST DISJOINT / DIFFERENT WEIGHTS RETURN CODE = `rc'"
if `rc' == 0 {
    dump_estimation, label("SUEST: disjoint samples with entirely different pweights")
}

display as text _newline "5D. SAME SAMPLE; PWEIGHTS DIFFER ONLY BY A CONSTANT MULTIPLIER"
estimates clear
regress y_linear x1 x2 [pweight=pw]
estimates store scale_a
dump_estimation, label("Weight-scale model using pw")
regress y_linear x1 x2 [pweight=pw10]
estimates store scale_b
dump_estimation, label("Weight-scale model using 10*pw")
capture noisily suest scale_a scale_b
local rc = _rc
display as result "SUEST DIFFERENT WEIGHT SCALE RETURN CODE = `rc'"
if `rc' == 0 {
    dump_estimation, label("SUEST: same sample with weights differing by factor of 10")
}

display as text _newline "=============================================================================="
display as text "6. FINAL FILES"
display as text "=============================================================================="
display as text "Upload both files:"
display as result "  suest_pweight_benchmark.log"
display as result "  suest_pweight_benchmark_data.dta"
display as text "Completed: " c(current_date) " " c(current_time)
display as text "=============================================================================="

log close pwbench
