* test_suest_r_stata_v22.do
* Follow-up to round 21: four checks of single models fit natively in Stata.
*   1. xtpoisson, re: robust and clustered standard errors of each model alone,
*      and suest2 with one and with both models (round 21: suest2's standard
*      errors differ from R's by 0.3% to 12%)
*   2. hetprobit: margins, dydx(x) after each model alone with vce(robust)
*      (round 21: suest2 gave a standard error 0.4% below R's for model A)
*   3. ivprobit: the predicted probability and the effect of x after each model
*      alone with vce(robust) (suest2 gives only the linear prediction)
*   4. mecompare after xtreg, fe (round 21: no effects could be estimated)
* Results print as SE|check|name|estimate|std. err. lines.
* File (same folder): v20_fixture.dta.
version 16
clear all
set more off
set linesize 255
set maxiter 50
capture log close
log using test_suest_r_stata_v22.log, text replace

* showse CHECK SOURCE   prints each estimate and standard error of SOURCE(b)
* and SOURCE(V), where SOURCE is e or r
capture program drop showse
program define showse
    args check src
    tempname b V
    matrix `b' = `src'(b)
    matrix `V' = `src'(V)
    local names : colfullnames `b'
    local k = colsof(`b')
    forvalues i = 1/`k' {
        local n : word `i' of `names'
        display as text "SE|`check'|`n'|" %21.0g `b'[1, `i'] "|" %21.0g sqrt(`V'[`i', `i'])
    }
end

foreach cmd in suest2 suest2_cleanup mecompare xtpoisson hetprobit ivprobit {
    capture noisily which `cmd'
}
display "maxiter = " c(maxiter)

use v20_fixture.dta, clear
xtset id t
local A x z f2 f3
local B x f2 f3
local it iterate(50)

* 1. xtpoisson, re
xtpoisson yc `A', re vce(robust) `it'
showse xtp_A_robust e
xtpoisson yc `A', re vce(cluster id) `it'
showse xtp_A_cluster e
xtpoisson yc `B', re vce(robust) `it'
showse xtp_B_robust e
quietly xtpoisson yc `A', re `it'
estimates store A
quietly xtpoisson yc `B', re `it'
estimates store B
capture noisily suest2 A
if !_rc showse xtp_suest2_A e
capture noisily suest2 A B
if !_rc showse xtp_suest2_AB e
capture noisily suest2_cleanup, force

* 2. hetprobit
hetprobit yh `A', het(w) vce(robust) `it'
showse het_A_robust e
margins, dydx(x)
showse het_A_dydx r
hetprobit yh `B', het(w) vce(robust) `it'
margins, dydx(x)
showse het_B_dydx r

* 3. ivprobit (maximum likelihood)
foreach m in A B {
    ivprobit yivb ``m'' (endog = instr1 instr2), vce(robust) `it'
    showse iv_`m'_robust e
    foreach opt in "pr" "pr fixedasf" "pr noasf" {
        local tag = subinstr("`opt'", " ", "_", .)
        capture noisily margins, predict(`opt')
        if !_rc showse iv_`m'_mean_`tag' r
        capture noisily margins, dydx(x) predict(`opt')
        if !_rc showse iv_`m'_dydx_`tag' r
    }
}

* 4. mecompare after xtreg, fe
quietly xtreg y `A', fe
estimates store A
quietly xtreg y `B', fe
estimates store B
capture noisily mecompare x, models(A B)
capture noisily mecompare x, models(A B) uncentered
capture noisily suest2_cleanup, force

log close
