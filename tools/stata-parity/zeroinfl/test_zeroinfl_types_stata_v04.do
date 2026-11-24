* test_zeroinfl_types_stata_v04.do
* The count-part mean, R type = "count", as Stata's exp() of the count
* equation, named for the outcome (art), and with predict(xb) beside it.
* Lines start ME| (margins: case, term, estimate, SE), FIT| (convergence, N).
* File (same folder): zeroinfl_biochem.dta (from v01).
* R values: zeroinfl_r_reference_v02.csv.
version 16
clear all
set more off
set linesize 255
capture log close
log using test_zeroinfl_types_stata_v04.log, text replace

foreach cmd in zip zinb margins {
    which `cmd'
}

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

foreach cmd in zip zinb {
    local M = upper("`cmd'")
    `cmd' art i.female kid5 ment, inflate(ment) vce(robust) iterate(100)
    display as text "FIT|`M'|converged=" e(converged) "|N=" e(N)
    capture noisily margins, expression(exp(predict(equation(art))))
    if !_rc dumpme mean_count|`M'
    capture noisily margins, dydx(ment) expression(exp(predict(equation(art))))
    if !_rc dumpme dydx_count|`M'
    capture noisily margins, expression(exp(predict(xb)))
    if !_rc dumpme mean_count_xb|`M'
    capture noisily margins, dydx(ment) expression(exp(predict(xb)))
    if !_rc dumpme dydx_count_xb|`M'
}

log close
