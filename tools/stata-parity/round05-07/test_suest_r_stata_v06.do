* R suest vs Stata for two new R features:
*   Part 1: glmmTMB random-intercept probit, cloglog, Gaussian, and Gamma-log.
*     gsem with Laplace integration is the exact comparison (glmmTMB uses
*     Laplace); suest2 with 12-point adaptive quadrature is close but not equal,
*     except for the Gaussian models (mixed), which are exact.
*   Part 2: lme4. lmer equals mixed; glmer(nAGQ = 7) equals 7-point
*     mode-curvature adaptive quadrature, intmethod(mcaghermite) intpoints(7).
* File (same folder): families_fixture_v01.dta. R values: families_r_reference_v01.csv.
version 16
clear all
set more off
set linesize 255
capture log close
log using test_suest_r_stata_v06.log, text replace

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

capture noisily which suest2
use families_fixture_v01.dta, clear
local aghq intmethod(mvaghermite) intpoints(12)

* 1. Probit pair
capture noisily gsem (yb1 <- x z M1[id], probit) (yb2 <- x M2[id], probit), ///
    cov(M1[id]*M2[id]@0) intmethod(laplace) vce(cluster id)
if !_rc dumpev probit_gsem
quietly meprobit yb1 x z || id:, `aghq'
estimates store A
quietly meprobit yb2 x || id:, `aghq'
estimates store B
capture noisily suest2 A B
if !_rc dumpev probit_suest2
capture noisily suest2_cleanup, force

* 2. Cloglog pair
capture noisily gsem (yb1 <- x z M1[id], cloglog) (yb2 <- x M2[id], cloglog), ///
    cov(M1[id]*M2[id]@0) intmethod(laplace) vce(cluster id)
if !_rc dumpev cloglog_gsem
quietly mecloglog yb1 x z || id:, `aghq'
estimates store A
quietly mecloglog yb2 x || id:, `aghq'
estimates store B
capture noisily suest2 A B
if !_rc dumpev cloglog_suest2
capture noisily suest2_cleanup, force

* 3. Gaussian pair (exact in both)
quietly mixed yc1 x z || id:, mle
estimates store A
quietly mixed yc2 x || id:, mle
estimates store B
capture noisily suest2 A B
if !_rc dumpev gaussian_suest2
capture noisily suest2_cleanup, force

* 4. Gamma (log link) pair
capture noisily gsem (yg1 <- x z M1[id], family(gamma) link(log)) ///
    (yg2 <- x M2[id], family(gamma) link(log)), ///
    cov(M1[id]*M2[id]@0) intmethod(laplace) vce(cluster id)
if !_rc dumpev gamma_gsem
quietly meglm yg1 x z || id:, family(gamma) link(log) `aghq'
estimates store A
quietly meglm yg2 x || id:, family(gamma) link(log) `aghq'
estimates store B
capture noisily suest2 A B
if !_rc dumpev gamma_suest2
capture noisily suest2_cleanup, force

* 5. Ordinary probit + random-intercept probit
capture noisily gsem (yb1 <- x z, probit) (yb2 <- x M2[id], probit), ///
    intmethod(laplace) vce(cluster id)
if !_rc dumpev probit_ordinary_gsem
quietly probit yb1 x z
estimates store A
quietly meprobit yb2 x || id:, `aghq'
estimates store B
capture noisily suest2 A B
if !_rc dumpev probit_ordinary_suest2
capture noisily suest2_cleanup, force

* 6. Gaussian + Gamma random-intercept models
capture noisily gsem (yc1 <- x z M1[id]) (yg2 <- x M2[id], family(gamma) link(log)), ///
    cov(M1[id]*M2[id]@0) intmethod(laplace) vce(cluster id)
if !_rc dumpev gaussian_gamma_gsem
quietly mixed yc1 x z || id:, mle
estimates store A
quietly meglm yg2 x || id:, family(gamma) link(log) `aghq'
estimates store B
capture noisily suest2 A B
if !_rc dumpev gaussian_gamma_suest2
capture noisily suest2_cleanup, force

* ---------------------------------------------------------------------------
* Part 2. lme4 counterparts
* ---------------------------------------------------------------------------
local mcagh intmethod(mcaghermite) intpoints(7)

* 7. Logit pair, 7-point mode-curvature adaptive quadrature (R: glmer nAGQ = 7)
quietly melogit yb1 x z || id:, `mcagh'
estimates store A
quietly melogit yb2 x || id:, `mcagh'
estimates store B
capture noisily suest2 A B
if !_rc dumpev logit_mcagh7_suest2
capture noisily suest2_cleanup, force

* 8. Probit pair, same quadrature
quietly meprobit yb1 x z || id:, `mcagh'
estimates store A
quietly meprobit yb2 x || id:, `mcagh'
estimates store B
capture noisily suest2 A B
if !_rc dumpev probit_mcagh7_suest2
capture noisily suest2_cleanup, force

* 9. mixed + ordinary logit (R: lmer + glm)
quietly mixed yc1 x z || id:, mle
estimates store A
quietly logit yb2 x
estimates store B
capture noisily suest2 A B
if !_rc dumpev mixed_logit_suest2
capture noisily suest2_cleanup, force

log close
