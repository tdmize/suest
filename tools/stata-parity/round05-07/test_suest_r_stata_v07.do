* R suest vs Stata for five new R features:
*   Part 1: survey systems that mix linear, binary, and count models (svy: with suest2)
*   Part 2: fixed-effects Poisson (xtpoisson, fe with suest2; R: fixest::fepois)
*   Part 3: generalized ordered logit (gologit2 with suest2; R: VGAM::vglm)
*   Part 4: Heckman selection, ML (heckman with suest2; R: sampleSelection)
*   Part 5: random-intercept ordered logit/probit (meologit/meoprobit with
*     7-point mode-curvature adaptive quadrature; R: ordinal::clmm, nAGQ = 7)
* Files (same folder): survey_fixture_v01.dta, fepois_fixture_v01.dta,
*   gologit_fixture_v01.dta, heckman_fixture_v01.dta, clmm_fixture_v01.dta.
* Part 3 needs gologit2 (ssc install gologit2).
version 16
clear all
set more off
set linesize 255
capture log close
log using test_suest_r_stata_v07.log, text replace

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
capture noisily which gologit2
* ---------------------------------------------------------------------------
* Part 1. Survey systems
* ---------------------------------------------------------------------------
use survey_fixture_v01.dta, clear
svyset psu [pweight = w], strata(strata)

* 1. svy: regress + svy: logit + svy: poisson
quietly svy: regress y x z
estimates store A
quietly svy: logit yb x z
estimates store B
quietly svy: poisson count x
estimates store C
capture noisily suest2 A B C
if !_rc dumpev regress_logit_poisson
capture noisily suest2_cleanup, force

* 2. Logit on a domain (sub == 1) + Poisson on everyone
quietly svy, subpop(sub): logit yb x
estimates store A
dumpev logit_sub_native
quietly svy: poisson count x
estimates store B
capture noisily suest2 A B
if !_rc dumpev logit_sub_poisson
capture noisily suest2_cleanup, force

* 3. svy: regress + svy: cloglog (R keeps svyglm's expected-information
*    covariance for cloglog; Stata uses observed information)
quietly svy: regress y x z
estimates store A
quietly svy: cloglog yb x z
estimates store B
capture noisily suest2 A B
if !_rc dumpev regress_cloglog
capture noisily suest2_cleanup, force

* ---------------------------------------------------------------------------
* Part 2. Fixed-effects Poisson
* ---------------------------------------------------------------------------
use fepois_fixture_v01.dta, clear
xtset id time
quietly xtpoisson y x z, fe
estimates store A
quietly xtpoisson y2 x, fe
estimates store B
capture noisily suest2 A B
if !_rc dumpev fepois
capture noisily suest2_cleanup, force

* Native robust fits: each model's block above should equal these
xtpoisson y x z, fe vce(robust)
xtpoisson y2 x, fe vce(robust)

* ---------------------------------------------------------------------------
* Part 3. Generalized ordered logit
* ---------------------------------------------------------------------------
use gologit_fixture_v01.dta, clear
capture noisily {
    quietly gologit2 y x z
    estimates store A
    quietly gologit2 y x z, pl(z)
    estimates store B
    capture noisily suest2 A B
    if !_rc dumpev gologit_full_partial
    capture noisily suest2_cleanup, force

    quietly gologit2 y x z, pl(z)
    estimates store A
    quietly logit yb x z
    estimates store B
    capture noisily suest2 A B
    if !_rc dumpev gologit_logit
    capture noisily suest2_cleanup, force
}
if _rc display as error "Part 3 stopped with r(" _rc ")"

* ---------------------------------------------------------------------------
* Part 4. Heckman selection (maximum likelihood)
* ---------------------------------------------------------------------------
use heckman_fixture_v01.dta, clear
quietly heckman y x z, select(s = x w)
estimates store A
quietly regress yc x z
estimates store B
capture noisily suest2 A B
if !_rc dumpev heckman_regress
capture noisily suest2_cleanup, force
heckman y x z, select(s = x w) vce(robust)

* ---------------------------------------------------------------------------
* Part 5. Random-intercept ordered logit and probit
* ---------------------------------------------------------------------------
use clmm_fixture_v01.dta, clear
local mcagh intmethod(mcaghermite) intpoints(7)
quietly meologit y x z || id:, `mcagh'
estimates store A
quietly meoprobit y x || id:, `mcagh'
estimates store B
capture noisily suest2 A B
if !_rc dumpev ologit_oprobit
capture noisily suest2_cleanup, force

quietly meologit y x z || id:, `mcagh'
estimates store A
quietly logit yb x
estimates store B
capture noisily suest2 A B
if !_rc dumpev ologit_logit
capture noisily suest2_cleanup, force

log close
