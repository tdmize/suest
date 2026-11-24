* test_suest_r_stata_v09.do
* R suest vs Stata, round 9: the new survey models (Batch 1). Each system is
* fit with svy: and combined with suest2; R fits the same models with svyolr(),
* svyVGAM::svy_vglm(), svysurvreg(), and svyivreg(). Each model's native
* svy: e(V) is printed too, so R's blocks can be checked without suest2.
* Every iterative model has iterate(50), and set maxiter 50 caps anything
* else. A FIT| line reports convergence for each model.
*   1-2: ordered logit and ordered probit pairs
*   3:   multinomial logit pair
*   4:   negative binomial pair
*   5:   generalized ordered logit pair (gologit2; needs ssc install gologit2)
*   6-8: tobit, interval regression, and Weibull (time metric) pairs
*   9:   two-stage least squares pair
*   10:  six different models in one system
*   11:  ordered logit + multinomial logit on a subpopulation
* File (same folder): survey_models_fixture_v01.dta.
* R values: survey_models_r_reference_v01.csv (survey_models_reference_v01.R).
version 16
clear all
set more off
set linesize 255
set maxiter 50
capture log close
log using test_suest_r_stata_v09.log, text replace

* Print e(b) and e(V) one element per line, at full precision
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

* fitstore NAME CASE: command   fits quietly, stores, reports convergence,
* and prints the native e(b) and e(V) as case CASE
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
    display as text "FIT|`case'|`name'|converged=" e(converged) "|N=" e(N)
    dumpev `case'_native_`name'
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

foreach cmd in suest2 suest2_cleanup gologit2 {
    capture noisily which `cmd'
}
display "maxiter = " c(maxiter)

use survey_models_fixture_v01.dta, clear
svyset psu [pweight = w], strata(strata)
stset t, failure(ev)

* 1. Ordered logit pair
fitstore A ologit_pair: svy: ologit o x z, iterate(50)
fitstore B ologit_pair: svy: ologit o x, iterate(50)
combine ologit_pair A B

* 2. Ordered probit pair
fitstore A oprobit_pair: svy: oprobit o x z, iterate(50)
fitstore B oprobit_pair: svy: oprobit o x, iterate(50)
combine oprobit_pair A B

* 3. Multinomial logit pair (base outcome 1, as R's refLevel = 1)
fitstore A mlogit_pair: svy: mlogit m x z, baseoutcome(1) iterate(50)
fitstore B mlogit_pair: svy: mlogit m x, baseoutcome(1) iterate(50)
combine mlogit_pair A B

* 4. Negative binomial pair (mean dispersion, Stata's default)
fitstore A nbreg_pair: svy: nbreg cnt x z, iterate(50)
fitstore B nbreg_pair: svy: nbreg cnt x, iterate(50)
combine nbreg_pair A B

* 5. Generalized ordered logit: all effects free + z parallel
fitstore A gologit_pair: svy: gologit2 o x z, iterate(50)
fitstore B gologit_pair: svy: gologit2 o x z, pl(z) iterate(50)
combine gologit_pair A B

* 6. Tobit pair, left-censored at 0
fitstore A tobit_pair: svy: tobit yc x z, ll(0) iterate(50)
fitstore B tobit_pair: svy: tobit yc x, ll(0) iterate(50)
combine tobit_pair A B

* 7. Interval regression pair (lo missing = left-censored, hi missing =
*    right-censored, lo == hi = exact)
fitstore A intreg_pair: svy: intreg lo hi x z, iterate(50)
fitstore B intreg_pair: svy: intreg lo hi x, iterate(50)
combine intreg_pair A B

* 8. Weibull pair, accelerated failure-time metric
fitstore A weibull_pair: svy: streg x z, distribution(weibull) time iterate(50)
fitstore B weibull_pair: svy: streg x, distribution(weibull) time iterate(50)
combine weibull_pair A B

* 9. Two-stage least squares pair
fitstore A ivreg_pair: svy: ivregress 2sls yiv z (endo = inst)
fitstore B ivreg_pair: svy: ivregress 2sls yiv (endo = inst)
combine ivreg_pair A B

* 10. Six different models in one system
fitstore A mixed_six: svy: regress yiv x z
fitstore B mixed_six: svy: logit yb x z, iterate(50)
fitstore C mixed_six: svy: ologit o x z, iterate(50)
fitstore D mixed_six: svy: mlogit m x z, baseoutcome(1) iterate(50)
fitstore E mixed_six: svy: nbreg cnt x z, iterate(50)
fitstore F mixed_six: svy: tobit yc x z, ll(0) iterate(50)
combine mixed_six A B C D E F

* 11. Ordered logit + multinomial logit on the subpopulation dom == 1
fitstore A domain_pair: svy, subpop(dom): ologit o x z, iterate(50)
fitstore B domain_pair: svy, subpop(dom): mlogit m x z, baseoutcome(1) iterate(50)
combine domain_pair A B

log close
