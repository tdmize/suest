* test_suest_r_stata_v18.do
* R suest vs Stata, round 18 (Batch 4): Gompertz and generalized gamma
* survival models (streg) and bivariate probit with different predictors in
* each equation (biprobit).
*   - Pairs of each model through suest2, from conventional fits.
*   - Each model's own vce(robust) fit, which should equal R's block for it.
*   - Average predictions with margins: median survival time after streg
*     (Stata's default prediction) and Pr(y1 = 1, y2 = 1) after biprobit.
* Every iterative model has iterate(50), and set maxiter 50 caps anything
* else. A FIT| line reports each model's sample and convergence; INFO| lines
* report log likelihoods.
* File (same folder): batch4_fixture_v01.dta. R values:
* batch4_r_reference_v01.csv (make_v18.R).
version 16
clear all
set more off
set linesize 255
set maxiter 50
capture log close
log using test_suest_r_stata_v18.log, text replace

* Print a coefficient vector and its covariance one element per line
capture program drop dumpbv
program define dumpbv
    args case b V
    local names : colfullnames `b'
    local k = colsof(`b')
    forvalues i = 1/`k' {
        local ri : word `i' of `names'
        display as text "EV|`case'|b|`ri'||" %21.0g `b'[1, `i']
    }
    forvalues i = 1/`k' {
        local ri : word `i' of `names'
        forvalues j = 1/`k' {
            local cj : word `j' of `names'
            display as text "EV|`case'|V|`ri'|`cj'|" %21.0g `V'[`i', `j']
        }
    }
end

* e(b) and e(V)
capture program drop dumpev
program define dumpev
    args case
    capture confirm matrix e(b)
    if _rc {
        display as error "EV|`case'|no results"
        exit
    }
    tempname b V
    matrix `b' = e(b)
    matrix `V' = e(V)
    dumpbv `case' `b' `V'
end

* fitstore NAME CASE: command   fits quietly, stores the results as NAME for
* suest2, and reports the sample and convergence
capture program drop fitstore
program define fitstore
    gettoken name 0 : 0
    gettoken case 0 : 0, parse(":")
    gettoken colon 0 : 0, parse(":")
    capture quietly `0'
    if _rc {
        display as error "FIT|`case'|`name'|failed with r(" _rc ")"
        ereturn clear
        exit
    }
    estimates store `name'
    display as text "FIT|`case'|`name'|converged=" e(converged) "|N=" e(N) "|ll=" %21.0g e(ll)
end

* combine CASE NAMES   runs suest2 on the stored models and prints the result
capture program drop combine
program define combine
    gettoken case 0 : 0
    capture noisily suest2 `0'
    if !_rc dumpev `case'
    else display as error "EV|`case'|suest2 failed with r(" _rc ")"
    capture noisily suest2_cleanup, force
end

* native CASE: command   one model's own fit (vce(robust)), its log
* likelihood, and its e(b) and e(V)
capture program drop native
program define native
    gettoken case 0 : 0, parse(":")
    gettoken colon 0 : 0, parse(":")
    capture quietly `0'
    if _rc {
        display as error "EV|`case'|failed with r(" _rc ")"
        ereturn clear
        exit
    }
    display as text "FIT|`case'|converged=" e(converged) "|N=" e(N) "|iterations=" e(ic)
    display as text "INFO|`case'|ll|" %21.0g e(ll)
    display as text "INFO|`case'|N|" e(N)
    if "`e(N_fail)'" != "" display as text "INFO|`case'|N_fail|" e(N_fail)
    dumpev `case'
end

* margdump CASE: margins ...   runs margins after the current fit and prints
* r(b) and r(V)
capture program drop margdump
program define margdump
    gettoken case 0 : 0, parse(":")
    gettoken colon 0 : 0, parse(":")
    capture noisily `0'
    if _rc {
        display as error "EV|`case'|margins failed with r(" _rc ")"
        exit
    }
    tempname b V
    matrix `b' = r(b)
    matrix `V' = r(V)
    dumpbv `case' `b' `V'
end

foreach cmd in suest2 suest2_cleanup stset streg biprobit probit margins {
    capture noisily which `cmd'
}
display "maxiter = " c(maxiter)

use batch4_fixture_v01.dta, clear
local opt iterate(50) nrtolerance(1e-10)

*-------------------------------------------------------------------------------
* 1. Gompertz: without and with a covariate on gamma (ancillary())
*-------------------------------------------------------------------------------
stset tm, failure(ev)
fitstore A gomp_pair: streg x z, distribution(gompertz) `opt'
fitstore B gomp_pair: streg x z, distribution(gompertz) ancillary(z) `opt'
combine gomp_pair A B
native gomp_pair_native_A: streg x z, distribution(gompertz) vce(robust) `opt'
margdump gomp_A_margins: margins, predict(median time)
native gomp_pair_native_B: streg x z, distribution(gompertz) ancillary(z) vce(robust) `opt'

*-------------------------------------------------------------------------------
* 2. Gompertz, bivariate probit, and probit in one system
*-------------------------------------------------------------------------------
fitstore A three_mixed: streg x z, distribution(gompertz) `opt'
fitstore B three_mixed: biprobit (y1 = x z) (y2 = x w), `opt'
fitstore C three_mixed: probit yp x w, `opt'
combine three_mixed A B C

*-------------------------------------------------------------------------------
* 3. Generalized gamma: without and with covariates on lnsigma (ancillary())
*    and kappa (anc2())
*-------------------------------------------------------------------------------
stset tg, failure(eg)
fitstore A ggamma_pair: streg x z, distribution(ggamma) `opt'
fitstore B ggamma_pair: streg x z, distribution(ggamma) ancillary(z) anc2(b) `opt'
combine ggamma_pair A B
native ggamma_pair_native_A: streg x z, distribution(ggamma) vce(robust) `opt'
native ggamma_pair_native_B: streg x z, distribution(ggamma) ancillary(z) anc2(b) vce(robust) `opt'
margdump ggamma_B_margins: margins, predict(median time)

*-------------------------------------------------------------------------------
* 4. Gompertz with delayed entry (entry > 0 for 40% of the rows)
*-------------------------------------------------------------------------------
stset tm, failure(ev) enter(time entry)
fitstore A entry_pair: streg x z, distribution(gompertz) `opt'
fitstore B entry_pair: streg x z, distribution(gompertz) ancillary(z) `opt'
combine entry_pair A B
native entry_pair_native_A: streg x z, distribution(gompertz) vce(robust) `opt'
native entry_pair_native_B: streg x z, distribution(gompertz) ancillary(z) vce(robust) `opt'

*-------------------------------------------------------------------------------
* 5. Bivariate probit with different predictors in each equation
*-------------------------------------------------------------------------------
fitstore A biprobit_pair: biprobit (y1 = x z) (y2 = x w), `opt'
fitstore B biprobit_pair: biprobit (y1 = x) (y2 = x w z), `opt'
combine biprobit_pair A B
native biprobit_pair_native_A: biprobit (y1 = x z) (y2 = x w), vce(robust) `opt'
margdump biprobit_A_margins: margins, predict(p11)
native biprobit_pair_native_B: biprobit (y1 = x) (y2 = x w z), vce(robust) `opt'

log close
