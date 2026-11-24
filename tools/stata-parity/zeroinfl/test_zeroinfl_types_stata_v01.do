* test_zeroinfl_types_stata_v01.do
* Checks suest's zero-inflated prediction types against Stata's predict()
* statistics after zip and zinb:
*   R type = "response"  <->  predict(n)
*   R type = "count"     <->  predict(ir)
*   R type = "zero"      <->  predict(pr)      (probability of a structural zero)
*   R type = "prob"      <->  predict(pr(k))   (k = 0, 1, 2 here)
* Part 1: zip and zinb with vce(robust), average predictions and dydx(ment).
* Part 2: official suest of the two models, the ZINB - ZIP difference in the
*         probability of a structural zero, by expression().
* Lines start ME| (margins: case, term, estimate, SE), FIT| (convergence, N).
* File (same folder): zeroinfl_biochem.dta (pscl's bioChemists data).
* R values: zeroinfl_r_reference_v01.csv (zeroinfl_types_reference_v01.R).
version 16
clear all
set more off
set linesize 255
capture log close
log using test_zeroinfl_types_stata_v01.log, text replace

* Print the last margins results: estimate and standard error per term
capture program drop dumpme
program define dumpme
    args case
    tempname t
    matrix `t' = r(table)
    local names : colfullnames `t'
    local k = colsof(`t')
    forvalues i = 1/`k' {
        local ri : word `i' of `names'
        display as text "ME|`case'|`ri'|" %21.0g `t'[1, `i'] "|" %21.0g `t'[2, `i']
    }
end

use zeroinfl_biochem.dta, clear
summarize art female kid5 ment

* Part 1: each model on its own, robust covariance
foreach cmd in zip zinb {
    local M = upper("`cmd'")
    `cmd' art i.female kid5 ment, inflate(ment) vce(robust) iterate(100)
    display as text "FIT|`M'|converged=" e(converged) "|N=" e(N)
    foreach s in n ir pr {
        capture noisily margins, predict(`s')
        if !_rc dumpme mean_`s'|`M'
        capture noisily margins, dydx(ment) predict(`s')
        if !_rc dumpme dydx_`s'|`M'
    }
    forvalues k = 0/2 {
        capture noisily margins, predict(pr(`k'))
        if !_rc dumpme mean_prk|`M'::`k'
        capture noisily margins, dydx(ment) predict(pr(`k'))
        if !_rc dumpme dydx_prk|`M'::`k'
    }
}

* Part 2: suest of the two models; equations are ZIP art, ZIP inflate,
* ZINB art, ZINB inflate, ZINB lnalpha
zip art i.female kid5 ment, inflate(ment) iterate(100)
estimates store ZIP
zinb art i.female kid5 ment, inflate(ment) iterate(100)
estimates store ZINB
capture noisily suest ZIP ZINB
if !_rc {
    display as text "FIT|suest|N=" e(N)
    local eqs : coleq e(b)
    display as text "FIT|suest|equations=" "`: list uniq eqs'"
    local diff invlogit(predict(equation(#4) xb)) - invlogit(predict(equation(#2) xb))
    capture noisily margins, expression(`diff')
    if !_rc dumpme diff_mean_pr
    capture noisily margins, dydx(ment) expression(`diff')
    if !_rc dumpme diff_dydx_pr
}

log close
