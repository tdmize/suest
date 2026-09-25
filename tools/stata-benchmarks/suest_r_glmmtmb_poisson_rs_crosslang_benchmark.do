version 16
clear all
set more off
set linesize 255
capture log close _all
log using suest_r_glmmtmb_poisson_rs_crosslang_benchmark.log, replace text name(glmmrs)
di as result "BENCHMARK_FILE_REVISION=1"
di as txt "SUEST R GLMMTMB POISSON RANDOM-INTERCEPT AND RANDOM-SLOPE BENCHMARK"
di as txt "Stata version: `c(stata_version)'"
import delimited using suest_r_glmmtmb_poisson_rs_crosslang_benchmark.csv, clear varnames(1)

* Primary components: identical Laplace likelihood and freely correlated effects.
foreach outcome in y1 y2 y1_partial y2_partial y1_left y2_right {
    di as txt "=== LAPLACE COMPONENT: `outcome' ==="
    quietly mepoisson `outcome' x z || id: x, covariance(unstructured) intmethod(laplace)
    assert e(converged) == 1
    di as txt "N=" e(N) " LL=" %21.15g e(ll)
    matrix list e(N_g), format(%21.15g)
    matrix list e(b), format(%21.15g)
    matrix list e(V), format(%21.15g)
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
di as result "POISSON_RS_COMPONENT_BENCHMARK_COMPLETE=1"

* Within each equation, intercept/slope covariance is free. Across equations,
* all four latent covariances are zero so the likelihood remains separable.
* Robust group scores still retain cross-equation sampling covariance.
foreach case in balanced partial disjoint {
    local first y1
    local second y2
    local vce "vce(robust)"
    if "`case'" == "partial" {
        local first y1_partial
        local second y2_partial
        local vce "vce(cluster higher)"
    }
    if "`case'" == "disjoint" {
        local first y1_left
        local second y2_right
    }
    di as txt "=== JOINT GSEM LAPLACE: `case' ==="
    quietly gsem (`first' <- x z I1[id]@1 c.x#S1[id]@1, family(poisson) link(log)) ///
        (`second' <- x z I2[id]@1 c.x#S2[id]@1, family(poisson) link(log)), ///
        cov(I1[id]*S1[id] I2[id]*S2[id] ///
            I1[id]*I2[id]@0 I1[id]*S2[id]@0 S1[id]*I2[id]@0 S1[id]*S2[id]@0) ///
        intmethod(laplace) `vce'
    assert e(converged) == 1
    di as txt "N=" e(N) " N_clust=" e(N_clust) " LL=" %21.15g e(ll)
    matrix list e(b), format(%21.15g)
    matrix list e(V), format(%21.15g)
    di as result "GSEM_LAPLACE_`case'_COMPLETE=1"
}
di as result "GLMMTMB_POISSON_RS_CROSSLANG_BENCHMARK_COMPLETE=1"
log close glmmrs
