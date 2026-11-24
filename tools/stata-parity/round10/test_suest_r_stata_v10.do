* test_suest_r_stata_v10.do
* R suest vs Stata, round 10: random-effects structures with several random
* slopes and with classes nested in schools (Batch 2), plus one survey check.
* Every iterative model has iterate(50), and set maxiter 50 caps anything
* else. A FIT| line reports convergence for each model.
*   Part 1: linear models (mixed), exact in both programs, combined with suest2.
*           Each model is also fit alone with vce(cluster school).
*   Part 2: nonlinear models. R fits these with the Laplace approximation, and
*           suest2 accepts only adaptive quadrature, so each model is fit alone
*           with intmethod(laplace) vce(cluster school); R's block for that
*           model should equal it.
*   Part 3: two joint gsem fits with intmethod(laplace) vce(cluster school),
*           which check the cross-model covariances of nonlinear models.
*   Part 4: survey complementary log-log pair (svy: cloglog), which R now
*           computes with observed information, as Stata does.
* Files (same folder): rs_fixture_v01.dta and survey_models_fixture_v01.dta
* (from round 9). R values: rs_r_reference_v01.csv (rs_reference_v01.R).
version 16
clear all
set more off
set linesize 255
set maxiter 50
capture log close
log using test_suest_r_stata_v10.log, text replace

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

* fitstore NAME CASE: command   fits quietly, stores the results as NAME for
* suest2, and reports convergence
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
end

* native CASE: command   fits quietly and prints e(b) and e(V) as case CASE
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
    display as text "FIT|`case'|converged=" e(converged) "|N=" e(N) "|clusters=" e(N_clust)
    dumpev `case'
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

foreach cmd in suest2 suest2_cleanup {
    capture noisily which `cmd'
}
display "maxiter = " c(maxiter)

use rs_fixture_v01.dta, clear

* ---------------------------------------------------------------------------
* Part 1. Linear models
* ---------------------------------------------------------------------------
* 1. Intercept and two slopes for school; intercept and one slope for school
*    with an intercept for classes nested in schools
fitstore A mixed_pair: mixed yc x1 x2 || school: x1 x2, covariance(unstructured) mle iterate(50)
fitstore B mixed_pair: mixed yc2 x1 x2 || school: x1, covariance(unstructured) || class:, mle iterate(50)
combine mixed_pair A B
native mixed_pair_cluster_A: mixed yc x1 x2 || school: x1 x2, covariance(unstructured) mle vce(cluster school) iterate(50)
native mixed_pair_cluster_B: mixed yc2 x1 x2 || school: x1, covariance(unstructured) || class:, mle vce(cluster school) iterate(50)

* 2. The nested linear model with an ordinary logit
fitstore A mixed_logit: mixed yc2 x1 x2 || school: x1, covariance(unstructured) || class:, mle iterate(50)
fitstore B mixed_logit: logit yb x1 x2, iterate(50)
combine mixed_logit A B
native mixed_logit_cluster_B: logit yb x1 x2, vce(cluster school) iterate(50)

* ---------------------------------------------------------------------------
* Part 2. Nonlinear models, one at a time, Laplace with clustered errors
* ---------------------------------------------------------------------------
local lap intmethod(laplace) vce(cluster school) iterate(50)

* 3. Logit: slope for school, classes nested in schools
native melogit_nested_cluster_A: melogit yb x1 x2 || school: x1, covariance(unstructured) || class:, `lap'

* 4. Poisson: intercept and two slopes for school
native mepoisson_slopes_cluster_A: mepoisson count x1 x2 || school: x1 x2, covariance(unstructured) `lap'

* 5. Probit: slope for school, classes nested in schools
native meprobit_nested_cluster_A: meprobit yb x1 x2 || school: x1, covariance(unstructured) || class:, `lap'

* 6. Poisson: intercept and two slopes for school, independent (diagonal)
native mepoisson_diag_cluster_A: mepoisson count x1 x2 || school: x1 x2, covariance(independent) `lap'

* 7. Complementary log-log: slope for school, classes nested in schools
native mecloglog_nested_cluster_A: mecloglog yb x1 || school: x1, covariance(unstructured) || class:, `lap'

* 8. Negative binomial (mean dispersion): slope for school, nested classes
native menbreg_nested_cluster_A: menbreg nb x1 || school: x1, covariance(unstructured) || class:, dispersion(mean) `lap'

* 9. Ordered logit: intercept and two slopes for school
native meologit_slopes_cluster_A: meologit yo x1 x2 || school: x1 x2, covariance(unstructured) `lap'

* 10. Ordered logit: slope for school, classes nested in schools
native meologit_nested_cluster_A: meologit yo x1 x2 || school: x1, covariance(unstructured) || class:, `lap'

* 11. Ordered probit: slope for school
native meoprobit_slopes_cluster_A: meoprobit yo x1 x2 || school: x1, covariance(unstructured) `lap'

* ---------------------------------------------------------------------------
* Part 3. Joint gsem fits (cross-model covariances)
* ---------------------------------------------------------------------------
* 12. Logit (slope for school, nested classes) + Poisson (two slopes for
*     school). Random effects of different models are uncorrelated.
native gsem_logit_poisson: gsem ///
    (yb <- x1 x2 M1[school]@1 c.x1#M2[school]@1 M3[school>class]@1, logit) ///
    (count <- x1 x2 N1[school]@1 c.x1#N2[school]@1 c.x2#N3[school]@1, poisson), ///
    covariance(M1[school]*M2[school] N1[school]*N2[school] N1[school]*N3[school] ///
        N2[school]*N3[school] M1[school]*N1[school]@0 M1[school]*N2[school]@0 ///
        M1[school]*N3[school]@0 M2[school]*N1[school]@0 M2[school]*N2[school]@0 ///
        M2[school]*N3[school]@0) ///
    intmethod(laplace) vce(cluster school) iterate(50)

* 13. Ordered logit (slope for school, nested classes) + ordinary logit
native gsem_ologit_logit: gsem ///
    (yo <- x1 x2 M1[school]@1 c.x1#M2[school]@1 M3[school>class]@1, ologit) ///
    (yb <- x1 x2, logit), ///
    covariance(M1[school]*M2[school]) ///
    intmethod(laplace) vce(cluster school) iterate(50)

* ---------------------------------------------------------------------------
* Part 4. Survey complementary log-log pair
* ---------------------------------------------------------------------------
use survey_models_fixture_v01.dta, clear
svyset psu [pweight = w], strata(strata)
fitstore A svy_cloglog_pair: svy: cloglog yb x z, iterate(50)
fitstore B svy_cloglog_pair: svy: cloglog yb x, iterate(50)
combine svy_cloglog_pair A B
native svy_cloglog_pair_cluster_A: svy: cloglog yb x z, iterate(50)
native svy_cloglog_pair_cluster_B: svy: cloglog yb x, iterate(50)

log close
