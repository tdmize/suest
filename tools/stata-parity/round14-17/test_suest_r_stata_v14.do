* test_suest_r_stata_v14.do
* R suest vs Stata, round 14: the Batch 3 models (round 13 stopped at the
* menbreg FIT line, where e(N_g) is a matrix; now read with ngroups).
*   Part 1: constant-dispersion negative binomial (nbreg and menbreg with
*           dispersion(constant)).
*   Part 2: population-averaged models (xtgee) with ar1 and unstructured
*           working correlations. Stata and R's geepack estimate these
*           correlations differently, so each model is fit twice: with the
*           correlation matrix R estimated, through corr(fixed), which should
*           reproduce R exactly; and with Stata's own corr(ar1) or
*           corr(unstructured), to see how far the two estimates are apart.
*           A last check shows what Stata does with gaps inside panels.
*   Part 3: generalized ordered models (gologit2) with link(cauchit) and
*           link(loglog), and one system mixing three kinds of model.
* Every iterative model has iterate(50), and set maxiter 50 caps anything
* else. A FIT| line reports each model's sample and convergence.
* File (same folder): batch3_fixture_v01.dta. R values:
* batch3_r_reference_v02.csv (make_v14.R, which also wrote the matrices below).
version 16
clear all
set more off
set linesize 255
set maxiter 50
capture log close
log using test_suest_r_stata_v14.log, text replace

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

* Print xtgee's working correlation e(R) and scale e(phi)
capture program drop dumpr
program define dumpr
    args case
    capture confirm matrix e(R)
    if _rc exit
    tempname R
    matrix `R' = e(R)
    local r = rowsof(`R')
    local c = colsof(`R')
    forvalues i = 1/`r' {
        forvalues j = 1/`c' {
            display as text "R|`case'|`i'|`j'|" %21.0g `R'[`i', `j']
        }
    }
    display as text "INFO|`case'|phi|" %21.0g e(phi)
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

* native CASE: command   fits quietly and prints e(b), e(V), the log
* likelihood, and for xtgee the working correlation, as case CASE
capture program drop native
program define native
    gettoken case 0 : 0, parse(":")
    gettoken colon 0 : 0, parse(":")
    capture quietly `0'
    if _rc {
        display as error "FIT|`case'|failed with r(" _rc ")"
        ereturn clear
        exit
    }
    ngroups
    local ng `r(ng)'
    display as text "FIT|`case'|converged=" e(converged) "|N=" e(N) "|N_g=`ng'" ///
        "|clusters=" e(N_clust)
    display as text "INFO|`case'|ll|" %21.0g e(ll)
    display as text "INFO|`case'|N|" %21.0g e(N)
    display as text "INFO|`case'|N_g|`ng'"
    dumpev `case'
    dumpr `case'
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

foreach cmd in suest2 suest2_cleanup gologit2 nbreg menbreg xtgee {
    capture noisily which `cmd'
}
display "maxiter = " c(maxiter)

use batch3_fixture_v01.dta, clear

* ---------------------------------------------------------------------------
* Part 1. Constant-dispersion negative binomial
* ---------------------------------------------------------------------------
* 1. A pair with constant dispersion: Var(y) = mu*(1 + delta)
fitstore A nb1_pair: nbreg cnt1 x z, dispersion(constant) iterate(50)
fitstore B nb1_pair: nbreg cnt1 x, dispersion(constant) iterate(50)
combine nb1_pair A B
native nb1_pair_native_A: nbreg cnt1 x z, dispersion(constant) vce(robust) iterate(50)
native nb1_pair_native_B: nbreg cnt1 x, dispersion(constant) vce(robust) iterate(50)

* 2. Constant dispersion, mean dispersion (the default), and Poisson
fitstore A nb_three: nbreg cnt1 x z, dispersion(constant) iterate(50)
fitstore B nb_three: nbreg cnt1 x z, dispersion(mean) iterate(50)
fitstore C nb_three: poisson cnt1 x z, iterate(50)
combine nb_three A B C

* 3. Random intercept for id, constant dispersion. R fits this with the
*    Laplace approximation, so Stata uses intmethod(laplace), alone, with
*    errors clustered on id.
native nb1_ri_native_A: menbreg cnt1 x z || id:, dispersion(constant) intmethod(laplace) vce(cluster id) iterate(50)

* ---------------------------------------------------------------------------
* Part 2. Population-averaged models (xtgee)
* ---------------------------------------------------------------------------
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

* 4. Linear, ar1
fitstore A gee_ar1: xtgee yg x z, family(gaussian) link(identity) corr(fixed R_ar1_A) `tol'
fitstore B gee_ar1: xtgee yg x, family(gaussian) link(identity) corr(fixed R_ar1_B) `tol'
combine gee_ar1 A B
native gee_ar1_fixed_A: xtgee yg x z, family(gaussian) link(identity) corr(fixed R_ar1_A) vce(robust) `tol'
native gee_ar1_fixed_B: xtgee yg x, family(gaussian) link(identity) corr(fixed R_ar1_B) vce(robust) `tol'
native gee_ar1_own_A: xtgee yg x z, family(gaussian) link(identity) corr(ar1) vce(robust) `tol'
native gee_ar1_own_B: xtgee yg x, family(gaussian) link(identity) corr(ar1) vce(robust) `tol'

* 5. Linear, unstructured
fitstore A gee_un: xtgee yg x z, family(gaussian) link(identity) corr(fixed R_un_A) `tol'
fitstore B gee_un: xtgee yg x, family(gaussian) link(identity) corr(fixed R_un_B) `tol'
combine gee_un A B
native gee_un_fixed_A: xtgee yg x z, family(gaussian) link(identity) corr(fixed R_un_A) vce(robust) `tol'
native gee_un_fixed_B: xtgee yg x, family(gaussian) link(identity) corr(fixed R_un_B) vce(robust) `tol'
native gee_un_own_A: xtgee yg x z, family(gaussian) link(identity) corr(unstructured) vce(robust) `tol'
native gee_un_own_B: xtgee yg x, family(gaussian) link(identity) corr(unstructured) vce(robust) `tol'

* 6. Logit with ar1 and Poisson with unstructured, in one system
fitstore A gee_mixed: xtgee yb x z, family(binomial) link(logit) corr(fixed R_mixed_A) `tol'
fitstore B gee_mixed: xtgee yc x z, family(poisson) link(log) corr(fixed R_mixed_B) `tol'
combine gee_mixed A B
native gee_mixed_fixed_A: xtgee yb x z, family(binomial) link(logit) corr(fixed R_mixed_A) vce(robust) `tol'
native gee_mixed_fixed_B: xtgee yc x z, family(poisson) link(log) corr(fixed R_mixed_B) vce(robust) `tol'
native gee_mixed_own_A: xtgee yb x z, family(binomial) link(logit) corr(ar1) vce(robust) `tol'
native gee_mixed_own_B: xtgee yc x z, family(poisson) link(log) corr(unstructured) vce(robust) `tol'

* 7. Panels that end early (short == 1): linear unstructured, logit ar1
fitstore A gee_short: xtgee yg x z if short, family(gaussian) link(identity) corr(fixed R_short_A) `tol'
fitstore B gee_short: xtgee yb x z if short, family(binomial) link(logit) corr(fixed R_short_B) `tol'
combine gee_short A B
native gee_short_fixed_A: xtgee yg x z if short, family(gaussian) link(identity) corr(fixed R_short_A) vce(robust) `tol'
native gee_short_fixed_B: xtgee yb x z if short, family(binomial) link(logit) corr(fixed R_short_B) vce(robust) `tol'
native gee_short_own_A: xtgee yg x z if short, family(gaussian) link(identity) corr(unstructured) vce(robust) `tol'
native gee_short_own_B: xtgee yb x z if short, family(binomial) link(logit) corr(ar1) vce(robust) `tol'

* 8. Gaps inside panels (gap == 1): without force, Stata leaves out panels
*    whose times are not equally spaced; with force, it uses the order of the
*    observations. R's geeglm keeps every panel.
native gap_waves_own_A: xtgee yg x z if gap, family(gaussian) link(identity) corr(ar1) vce(robust) `tol'
native gap_order_own_A: xtgee yg x z if gap, family(gaussian) link(identity) corr(ar1) force vce(robust) `tol'
native gap_un_own_A: xtgee yg x z if gap, family(gaussian) link(identity) corr(unstructured) vce(robust) `tol'

* ---------------------------------------------------------------------------
* Part 3. Generalized ordered models (gologit2)
* ---------------------------------------------------------------------------
* 9. Cauchit link; z parallel in model A
fitstore A cauchit_pair: gologit2 oc x z, link(cauchit) pl(z) iterate(50)
fitstore B cauchit_pair: gologit2 oc x, link(cauchit) iterate(50)
combine cauchit_pair A B
native cauchit_pair_native_A: gologit2 oc x z, link(cauchit) pl(z) robust iterate(50)
native cauchit_pair_native_B: gologit2 oc x, link(cauchit) robust iterate(50)

* 10. Log-log link; z parallel in model A
fitstore A loglog_pair: gologit2 og x z, link(loglog) pl(z) iterate(50)
fitstore B loglog_pair: gologit2 og x, link(loglog) iterate(50)
combine loglog_pair A B
native loglog_pair_native_A: gologit2 og x z, link(loglog) pl(z) robust iterate(50)
native loglog_pair_native_B: gologit2 og x, link(loglog) robust iterate(50)

* 11. Cauchit gologit2 + constant-dispersion nbreg + logit in one system
fitstore A three_mixed: gologit2 oc x z, link(cauchit) pl(z) iterate(50)
fitstore B three_mixed: nbreg cnt1 x z, dispersion(constant) iterate(50)
fitstore C three_mixed: logit yb x z, iterate(50)
combine three_mixed A B C

log close
