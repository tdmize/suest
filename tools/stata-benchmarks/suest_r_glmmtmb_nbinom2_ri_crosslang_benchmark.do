version 16
clear all
set more off
set linesize 255

capture log close _all
log using suest_r_glmmtmb_nbinom2_ri_crosslang_benchmark.log, replace text name(glmmnb2)
di as result "BENCHMARK_FILE_REVISION=1"
di as txt "SUEST R GLMMTMB RANDOM-INTERCEPT NB2 CROSS-LANGUAGE BENCHMARK"
di as txt "Stata version: `c(stata_version)'"
di as txt "NB2: variance = mu + alpha*mu^2; R log_phi = -ln(alpha)."

import delimited using suest_r_glmmtmb_nbinom2_ri_crosslang_benchmark.csv, clear varnames(1)

* Primary component gate: all six fits use the same Laplace likelihood as R.
foreach outcome in y1 y2 y1_partial y2_partial y1_left y2_right {
    di as txt "=== LAPLACE COMPONENT: `outcome' ==="
    quietly menbreg `outcome' x z || id:, dispersion(mean) intmethod(laplace)
    assert e(converged) == 1
    di as txt "N=" e(N) " LL=" %21.15g e(ll)
    matrix list e(N_g), format(%21.15g)
    matrix list e(b), format(%21.15g)
    matrix list e(V), format(%21.15g)
    estimates store nb2_`outcome'
    if inlist("`outcome'", "y1", "y2") {
        di as txt "MARGINAL MEAN: `outcome'"
        quietly margins, predict(mu marginal)
        matrix list r(b), format(%21.15g)
        matrix list r(V), format(%21.15g)
        di as txt "MARGINAL SLOPE X: `outcome'"
        quietly margins, dydx(x) predict(mu marginal)
        matrix list r(b), format(%21.15g)
        matrix list r(V), format(%21.15g)
    }
}

* Primary joint gate: zero covariance between latent effects keeps the fitted
* likelihood separable; robust scores retain covariance across equations.
di as txt "=== BALANCED JOINT GSEM, LAPLACE ==="
quietly gsem (y1 <- x z M1[id], family(nbinomial mean) link(log)) ///
    (y2 <- x z M2[id], family(nbinomial mean) link(log)), ///
    cov(M1[id]*M2[id]@0) intmethod(laplace) vce(robust)
assert e(converged) == 1
capture noisily di as txt "N=" e(N) " N_clust=" e(N_clust) " LL=" %21.15g e(ll)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
di as result "GSEM_LAPLACE_BALANCED_COMPLETE=1"

di as txt "=== PARTIAL OVERLAP, HIGHER-CLUSTER GSEM LAPLACE ==="
quietly gsem (y1_partial <- x z M1[id], family(nbinomial mean) link(log)) ///
    (y2_partial <- x z M2[id], family(nbinomial mean) link(log)), ///
    cov(M1[id]*M2[id]@0) intmethod(laplace) vce(cluster higher)
assert e(converged) == 1
capture noisily di as txt "N=" e(N) " N_clust=" e(N_clust) " LL=" %21.15g e(ll)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
di as result "GSEM_LAPLACE_PARTIAL_COMPLETE=1"

di as txt "=== DISJOINT-GROUP GSEM, LAPLACE ==="
quietly gsem (y1_left <- x z M1[id], family(nbinomial mean) link(log)) ///
    (y2_right <- x z M2[id], family(nbinomial mean) link(log)), ///
    cov(M1[id]*M2[id]@0) intmethod(laplace) vce(robust)
assert e(converged) == 1
capture noisily di as txt "N=" e(N) " N_clust=" e(N_clust) " LL=" %21.15g e(ll)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
di as result "GSEM_LAPLACE_DISJOINT_COMPLETE=1"
di as result "NBINOM2_PRIMARY_LAPLACE_BENCHMARK_COMPLETE=1"

* Secondary approximation comparison. Unsupported suest2 routes are logged
* without interrupting the completed primary menbreg/gsem gate.
capture which suest2
local has_suest2 = (_rc == 0)
if `has_suest2' {
    di as txt "=== OPTIONAL SUEST2 LAPLACE CHECK ==="
    capture noisily suest2 nb2_y1 nb2_y2
    di as result "SUEST2_LAPLACE_RC=" _rc
    capture noisily suest2_cleanup, force
    foreach case in balanced partial disjoint {
        local first y1
        local second y2
        local options ""
        if "`case'" == "partial" {
            local first y1_partial
            local second y2_partial
            local options ", cluster(higher)"
        }
        if "`case'" == "disjoint" {
            local first y1_left
            local second y2_right
        }
        di as txt "=== OPTIONAL ADAPTIVE 12-POINT SUEST2: `case' ==="
        quietly menbreg `first' x z || id:, dispersion(mean) intmethod(mvaghermite) intpoints(12)
        assert e(converged) == 1
        estimates store nb2_adaptive1
        quietly menbreg `second' x z || id:, dispersion(mean) intmethod(mvaghermite) intpoints(12)
        assert e(converged) == 1
        estimates store nb2_adaptive2
        capture noisily suest2 nb2_adaptive1 nb2_adaptive2 `options'
        local result = _rc
        di as result "SUEST2_ADAPTIVE_RC=" `result'
        if `result' == 0 {
            capture noisily di as txt "N=" e(N) " N_clust=" e(N_clust)
            matrix list e(b), format(%21.15g)
            matrix list e(V), format(%21.15g)
        }
        capture noisily suest2_cleanup, force
    }
}
else {
    di as txt "SUEST2_NOT_INSTALLED: optional comparison skipped."
}

di as result "GLMMTMB_NBINOM2_RI_CROSSLANG_BENCHMARK_COMPLETE=1"
log close glmmnb2
