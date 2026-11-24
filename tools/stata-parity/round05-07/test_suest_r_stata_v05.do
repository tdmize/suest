* R suest vs Stata, two new R features:
*   Part 1: pooled marginal effects after multiple imputation (suest2 + mecompare)
*   Part 2: ordinary models combined with multilevel models (suest2 1.1.0; gsem)
* Files (same folder): mi_fixture_v01.dta, mi_r_reference_v01.csv, combined_fixture_v01.dta
* Part 1 needs mimrgns (ssc install mimrgns).
version 16
clear all
set more off
set linesize 255
capture log close
log using test_suest_r_stata_v05.log, text replace

* Print e(b) and e(V) one element per line, at full precision
capture program drop dumpev
program define dumpev
    args case
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

foreach cmd in suest2 mecompare mimrgns {
    capture noisily which `cmd'
}

* ---------------------------------------------------------------------------
* Part 1. Multiple imputation: five imputations made in R (mice), flong form
* ---------------------------------------------------------------------------
import delimited using mi_r_reference_v01.csv, varnames(1) stringcols(1) asdouble clear
format estimate std_error %20.0g
list, noobs clean

capture noisily {
    use mi_fixture_v01.dta, clear
    mi import flong, m(imp_m) id(imp_id) imputed(x z) clear
    mi describe

    mi estimate, post: logit y c.x, vce(robust)
    estimates store mi1
    mi estimate, post: logit y c.x c.z, vce(robust)
    estimates store mi2

    * 1a. Pooled joint coefficients and covariance (R rows 1-5; Cov(mi1 x, mi2 x) = 0.011480098876)
    capture noisily suest2 mi1 mi2
    if !_rc dumpev mi_suest2
    capture noisily suest2_cleanup, force

    * 1b. Pooled marginal effects of x (centered unit change) and their difference (R rows 6-8)
    capture noisily mecompare x, models(mi1 mi2) decimals(8) commands
    if !_rc dumpev mi_mecompare
    capture noisily suest2_cleanup, force
}
if _rc display as error "Part 1 stopped with r(" _rc ")"

* ---------------------------------------------------------------------------
* Part 2. Ordinary + multilevel systems, clustered on id (80 groups)
* ---------------------------------------------------------------------------
estimates clear
use combined_fixture_v01.dta, clear

* 2a. regress + mixed
quietly regress yc x z
estimates store A
quietly mixed yc2 x || id:, mle
estimates store B
capture noisily suest2 A B
if !_rc dumpev regress_mixed
capture noisily suest2_cleanup, force

* 2b. logit on a subsample + mixed
quietly logit y1 x z if sub == 1
estimates store A
quietly mixed yc2 x || id:, mle
estimates store B
capture noisily suest2 A B
if !_rc dumpev logit_mixed_sub
capture noisily suest2_cleanup, force

* 2c. logit + melogit: joint gsem with Laplace integration (the method R uses)
capture noisily gsem (y1 <- x z, logit) (y2 <- x M1[id], logit), intmethod(laplace) vce(cluster id)
if !_rc dumpev gsem_logit_melogit

* 2d. logit + melogit through suest2 (12-point adaptive quadrature; close to, not equal to, Laplace)
quietly logit y1 x z
estimates store A
quietly melogit y2 x || id:, intmethod(mvaghermite) intpoints(12)
estimates store B
capture noisily suest2 A B
if !_rc dumpev suest2_logit_melogit
capture noisily suest2_cleanup, force

* 2e. regress + logit + melogit, three models
quietly regress yc x z
estimates store A
quietly logit y1 x z
estimates store B
quietly melogit y2 x || id:, intmethod(mvaghermite) intpoints(12)
estimates store C
capture noisily suest2 A B C
if !_rc dumpev suest2_regress_logit_melogit
capture noisily suest2_cleanup, force

* 2f. Native clustered fits: each model's block in 2a-2b should equal these
regress yc x z, vce(cluster id)
logit y1 x z if sub == 1, vce(cluster id)
capture noisily mixed yc2 x || id:, mle vce(cluster id)

log close
