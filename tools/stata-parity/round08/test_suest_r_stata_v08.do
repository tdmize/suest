* test_suest_r_stata_v08.do
* R suest vs Stata, round 8. Every iterative model has iterate(50), and
* set maxiter 50 caps anything else, so nothing can run forever. A model that
* stops at 50 iterations prints converged=0 on its FIT| line and is skipped.
*   Part 1: pooled marginal effects after multiple imputation (mecompare;
*           v05 stopped on decimals(8))
*   Part 2: binary random-intercept pairs with 25-point adaptive quadrature
*           (R: glmer nAGQ = 25), so both programs use near-exact integrals
*   Part 3: ordinary + random-intercept logit, 25-point quadrature
*   Part 4: Poisson random-intercept pair, 25-point quadrature
*   Part 5: random-intercept ordered logit/probit, 25-point quadrature
*           (v07 Part 5 used mcaghermite, which suest2 refuses)
*   Part 6: Gaussian random slopes, mixed with covariance(unstructured)
*   Part 7: Gamma-log random intercept, meglm with Laplace (native e(V))
*   Part 8: cloglog pair, gsem with Laplace and vce(oim) (diagnostic)
* Files (same folder): mi_fixture_v01.dta, families_fixture_v01.dta,
*   combined_fixture_v01.dta, clmm_fixture_v01.dta, v08_fixture_v01.dta.
* R values: v08_r_reference_v01.csv and mi_r_reference_v01.csv.
version 16
clear all
set more off
set linesize 255
set maxiter 50
capture log close
log using test_suest_r_stata_v08.log, text replace

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

* Print every r() matrix one element per line, at full precision
capture program drop dumpr
program define dumpr
    args case
    local mats : r(matrices)
    foreach m of local mats {
        tempname M
        matrix `M' = r(`m')
        local rn : rowfullnames `M'
        local cn : colfullnames `M'
        forvalues i = 1/`=rowsof(`M')' {
            local ri : word `i' of `rn'
            forvalues j = 1/`=colsof(`M')' {
                local cj : word `j' of `cn'
                display as text "RM|`case'|`m'|`ri'|`cj'|" %21.0g `M'[`i', `j']
            }
        }
    }
end

* fitstore NAME: command   fits quietly, stores, and reports convergence
capture program drop fitstore
program define fitstore
    gettoken name 0 : 0, parse(":")
    gettoken colon 0 : 0, parse(":")
    capture quietly `0'
    if _rc {
        display as error "FIT|`name'|failed with r(" _rc ")"
        ereturn clear
        exit
    }
    estimates store `name'
    display as text "FIT|`name'|converged=" e(converged) "|iterations=" e(ic) "|N=" e(N)
end

foreach cmd in suest2 mecompare mimrgns {
    capture noisily which `cmd'
}
display "maxiter = " c(maxiter)

* ---------------------------------------------------------------------------
* Part 1. Multiple imputation: pooled marginal effects of x (R rows 6-8 of
*         mi_r_reference_v01.csv: centered unit change for mi1 and mi2, and
*         their difference)
* ---------------------------------------------------------------------------
import delimited using mi_r_reference_v01.csv, varnames(1) stringcols(1) asdouble clear
format estimate std_error %20.0g
list, noobs clean

capture noisily {
    use mi_fixture_v01.dta, clear
    mi import flong, m(imp_m) id(imp_id) imputed(x z) clear
    mi estimate, post: logit y c.x, vce(robust) iterate(50)
    estimates store mi1
    mi estimate, post: logit y c.x c.z, vce(robust) iterate(50)
    estimates store mi2
    capture noisily mecompare x, models(mi1 mi2) decimals(7)
    if !_rc {
        return list
        dumpr mi_mecompare
    }
    capture noisily suest2_cleanup, force
}
if _rc display as error "Part 1 stopped with r(" _rc ")"

* ---------------------------------------------------------------------------
* Part 2. Binary random-intercept pairs, 25-point adaptive quadrature
* ---------------------------------------------------------------------------
use families_fixture_v01.dta, clear
local aghq intmethod(mvaghermite) intpoints(25) iterate(50)

foreach cmd in melogit meprobit mecloglog {
    estimates clear
    fitstore A: `cmd' yb1 x z || id:, `aghq'
    fitstore B: `cmd' yb2 x || id:, `aghq'
    capture noisily suest2 A B
    if !_rc dumpev `cmd'_aghq25
    capture noisily suest2_cleanup, force
}

estimates clear
fitstore A: probit yb1 x z, iterate(50)
fitstore B: meprobit yb2 x || id:, `aghq'
capture noisily suest2 A B
if !_rc dumpev probit_ordinary_aghq25
capture noisily suest2_cleanup, force

* ---------------------------------------------------------------------------
* Part 3. Ordinary + random-intercept logit, 25-point adaptive quadrature
* ---------------------------------------------------------------------------
use combined_fixture_v01.dta, clear
estimates clear
fitstore A: logit y1 x z, iterate(50)
fitstore B: melogit y2 x || id:, `aghq'
capture noisily suest2 A B
if !_rc dumpev logit_melogit_aghq25
capture noisily suest2_cleanup, force

estimates clear
quietly regress yc x z
estimates store A
fitstore B: logit y1 x z, iterate(50)
fitstore C: melogit y2 x || id:, `aghq'
capture noisily suest2 A B C
if !_rc dumpev regress_logit_melogit_aghq25
capture noisily suest2_cleanup, force

* ---------------------------------------------------------------------------
* Part 4. Poisson random-intercept pair, 25-point adaptive quadrature
* ---------------------------------------------------------------------------
use v08_fixture_v01.dta, clear
estimates clear
fitstore A: mepoisson yp1 x z || id:, `aghq'
fitstore B: mepoisson yp2 x || id:, `aghq'
capture noisily suest2 A B
if !_rc dumpev poisson_aghq25
capture noisily suest2_cleanup, force

* ---------------------------------------------------------------------------
* Part 5. Random-intercept ordered logit/probit, 25-point adaptive quadrature
* ---------------------------------------------------------------------------
use clmm_fixture_v01.dta, clear
estimates clear
fitstore A: meologit y x z || id:, `aghq'
fitstore B: meoprobit y x || id:, `aghq'
capture noisily suest2 A B
if !_rc dumpev ologit_oprobit_aghq25
capture noisily suest2_cleanup, force

estimates clear
fitstore A: meologit y x z || id:, `aghq'
fitstore B: logit yb x, iterate(50)
capture noisily suest2 A B
if !_rc dumpev ologit_logit_aghq25
capture noisily suest2_cleanup, force

* ---------------------------------------------------------------------------
* Part 6. Gaussian random intercept and slope (exact in both programs)
* ---------------------------------------------------------------------------
use v08_fixture_v01.dta, clear
estimates clear
fitstore A: mixed ys1 x z || id: x, covariance(unstructured) mle iterate(50)
fitstore B: mixed ys2 x || id: x, covariance(unstructured) mle iterate(50)
capture noisily suest2 A B
if !_rc dumpev mixed_slopes
capture noisily suest2_cleanup, force

* ---------------------------------------------------------------------------
* Part 7. Gamma-log random intercept with Laplace integration, one model at a
*         time (native e(V); R: glmmTMB model-based covariance)
* ---------------------------------------------------------------------------
use families_fixture_v01.dta, clear
estimates clear
fitstore G1: meglm yg1 x z || id:, family(gamma) link(log) intmethod(laplace) iterate(50)
dumpev gamma_laplace_yg1
fitstore G2: meglm yg2 x || id:, family(gamma) link(log) intmethod(laplace) iterate(50)
dumpev gamma_laplace_yg2

* ---------------------------------------------------------------------------
* Part 8. Diagnostic: cloglog pair with gsem, Laplace, model-based vce(oim)
* ---------------------------------------------------------------------------
estimates clear
capture noisily gsem (yb1 <- x z M1[id], cloglog) (yb2 <- x M2[id], cloglog), ///
    cov(M1[id]*M2[id]@0) intmethod(laplace) iterate(50)
if !_rc {
    display as text "FIT|gsem|converged=" e(converged) "|iterations=" e(ic)
    dumpev cloglog_laplace_oim
}

log close
