* test_suest_r_stata_v17.do
* R suest vs Stata, round 17: the four population-averaged systems through
* suest2. As in the earlier GEE benchmark (suest_r_gee_crosslang_benchmark.do),
* the models are conventional xtreg, xtlogit, and xtpoisson fits with pa, and
* suest2 supplies the joint robust covariance (round 16 used vce(robust),
* which suest2 refuses for population-averaged models).
* Every iterative model has iterate(50), and set maxiter 50 caps anything
* else. A FIT| line reports each model's sample.
* File (same folder): batch3_fixture_v01.dta. R values:
* batch3_r_reference_v02.csv (make_v14.R, which also wrote the matrices below).
version 16
clear all
set more off
set linesize 255
set maxiter 50
capture log close
log using test_suest_r_stata_v17.log, text replace

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

* ngroups: e(N_g) as one number in r(ng) (a matrix after multilevel commands)
capture program drop ngroups
program define ngroups, rclass
    tempname G
    local ng .
    capture local ng = e(N_g)
    capture matrix `G' = e(N_g)
    if !_rc local ng = `G'[1, 1]
    return local ng "`ng'"
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
    ngroups
    display as text "FIT|`case'|`name'|converged=" e(converged) "|N=" e(N) "|N_g=`r(ng)'"
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

foreach cmd in suest2 suest2_cleanup xtreg xtlogit xtpoisson xtgee {
    capture noisily which `cmd'
}
display "maxiter = " c(maxiter)

use batch3_fixture_v01.dta, clear
xtset id time

* Working correlations estimated by R (geepack), one per model
matrix R_ar1_A = (1, 0.63451033658162737, 0.40260336722893003, 0.25545599804932489, 0.1620894713040727 \ 0.63451033658162737, 1, 0.63451033658162737, 0.40260336722893003, 0.25545599804932489 \ 0.40260336722893003, 0.63451033658162737, 1, 0.63451033658162737, 0.40260336722893003 \ 0.25545599804932489, 0.40260336722893003, 0.63451033658162737, 1, 0.63451033658162737 \ 0.1620894713040727, 0.25545599804932489, 0.40260336722893003, 0.63451033658162737, 1)
matrix R_ar1_B = (1, 0.59912285089421524, 0.35894819046361209, 0.21505406319387904, 0.12884380343710153 \ 0.59912285089421524, 1, 0.59912285089421524, 0.35894819046361209, 0.21505406319387904 \ 0.35894819046361209, 0.59912285089421524, 1, 0.59912285089421524, 0.35894819046361209 \ 0.21505406319387904, 0.35894819046361209, 0.59912285089421524, 1, 0.59912285089421524 \ 0.12884380343710153, 0.21505406319387904, 0.35894819046361209, 0.59912285089421524, 1)
matrix R_un_A = (1, 0.59122020064386749, 0.34727748864490326, 0.24460651421158977, 0.26588484105976901 \ 0.59122020064386749, 1, 0.54368118817358879, 0.38960492015421716, 0.3652685185582486 \ 0.34727748864490326, 0.54368118817358879, 1, 0.5491416013068342, 0.44215331541017433 \ 0.24460651421158977, 0.38960492015421716, 0.5491416013068342, 1, 0.66201615654943113 \ 0.26588484105976901, 0.3652685185582486, 0.44215331541017433, 0.66201615654943113, 1)
matrix R_un_B = (1, 0.57435244035993749, 0.3414670432372755, 0.27369174122190787, 0.23690419185309941 \ 0.57435244035993749, 1, 0.4574492867544574, 0.29900437141486608, 0.28229659003146224 \ 0.3414670432372755, 0.4574492867544574, 1, 0.50100906796880729, 0.42814224854456762 \ 0.27369174122190787, 0.29900437141486608, 0.50100906796880729, 1, 0.64420416708421246 \ 0.23690419185309941, 0.28229659003146224, 0.42814224854456762, 0.64420416708421246, 1)
matrix R_mixed_A = (1, 0.13377544512283135, 0.017895869717811663, 0.0023940279373604535, 0.00032026215295688851 \ 0.13377544512283135, 1, 0.13377544512283135, 0.017895869717811663, 0.0023940279373604535 \ 0.017895869717811663, 0.13377544512283135, 1, 0.13377544512283135, 0.017895869717811663 \ 0.0023940279373604535, 0.017895869717811663, 0.13377544512283135, 1, 0.13377544512283135 \ 0.00032026215295688851, 0.0023940279373604535, 0.017895869717811663, 0.13377544512283135, 1)
matrix R_mixed_B = (1, 0.10699477446500912, 0.081776805181793019, -0.10000646511737721, 0.00092395361014768697 \ 0.10699477446500912, 1, 0.052163792404524865, 0.089548005190274282, -0.0502066897913339 \ 0.081776805181793019, 0.052163792404524865, 1, -0.033212453389608126, -0.097086620921867539 \ -0.10000646511737721, 0.089548005190274282, -0.033212453389608126, 1, 0.27170968567965176 \ 0.00092395361014768697, -0.0502066897913339, -0.097086620921867539, 0.27170968567965176, 1)
matrix R_short_A = (1, 0.61747776045534886, 0.3677423758993793, 0.23726945131204752, 0.3087041928693508 \ 0.61747776045534886, 1, 0.57356038264722686, 0.39586412064297988, 0.35192684911332106 \ 0.3677423758993793, 0.57356038264722686, 1, 0.58079168487262312, 0.51611716405428321 \ 0.23726945131204752, 0.39586412064297988, 0.58079168487262312, 1, 0.55172934361476522 \ 0.3087041928693508, 0.35192684911332106, 0.51611716405428321, 0.55172934361476522, 1)
matrix R_short_B = (1, 0.10425790186116769, 0.010869710100492874, 0.0011332531689165293, 0.00011815059766875681 \ 0.10425790186116769, 1, 0.10425790186116769, 0.010869710100492874, 0.0011332531689165293 \ 0.010869710100492874, 0.10425790186116769, 1, 0.10425790186116769, 0.010869710100492874 \ 0.0011332531689165293, 0.010869710100492874, 0.10425790186116769, 1, 0.10425790186116769 \ 0.00011815059766875681, 0.0011332531689165293, 0.010869710100492874, 0.10425790186116769, 1)

local tol tolerance(1e-9) iterate(50)

* 1. Linear, ar1
fitstore A gee_ar1: xtreg yg x z, pa corr(fixed R_ar1_A) `tol'
fitstore B gee_ar1: xtreg yg x, pa corr(fixed R_ar1_B) `tol'
combine gee_ar1 A B

* 2. Linear, unstructured
fitstore A gee_un: xtreg yg x z, pa corr(fixed R_un_A) `tol'
fitstore B gee_un: xtreg yg x, pa corr(fixed R_un_B) `tol'
combine gee_un A B

* 3. Logit with ar1 and Poisson with unstructured
fitstore A gee_mixed: xtlogit yb x z, pa corr(fixed R_mixed_A) `tol'
fitstore B gee_mixed: xtpoisson yc x z, pa corr(fixed R_mixed_B) `tol'
combine gee_mixed A B

* 4. Panels that end early (short == 1): linear unstructured, logit ar1
fitstore A gee_short: xtreg yg x z if short, pa corr(fixed R_short_A) `tol'
fitstore B gee_short: xtlogit yb x z if short, pa corr(fixed R_short_B) `tol'
combine gee_short A B

log close
