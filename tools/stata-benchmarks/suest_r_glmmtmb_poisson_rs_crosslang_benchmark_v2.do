version 16
clear all
set more off
set linesize 255
capture log close _all
log using suest_r_glmmtmb_poisson_rs_crosslang_benchmark_v2.log, replace text name(glmmrs)
di as result "BENCHMARK_FILE_REVISION=2"
di as txt "POISSON RANDOM-SLOPE LAPLACE BENCHMARK: INITIALIZATION DIAGNOSTICS"
di as txt "Stata version: `c(stata_version)'"
di as txt "Same data and model as revision 1; convergence tolerances are unchanged."
di as txt "Each Laplace attempt is capped at 200 iterations; seed fits at 100."
import delimited using suest_r_glmmtmb_poisson_rs_crosslang_benchmark.csv, clear varnames(1)
isid id time
assert _N == 800

* All attempts show their iteration logs. A failed cold Laplace fit is retried
* from a Stata adaptive-quadrature fit of the SAME model and sample. Quadrature
* is only an initializer; only a converged final Laplace fit counts below.
* No R coefficients are supplied. Clear e() before every attempt so that a
* failed command cannot leave an earlier fit masquerading as its own result.
program define rs_fit, eclass
    syntax, First(name) [Second(name) Cluster(name)]
    local points 7
    if "`second'" == "" {
        local model "mepoisson `first' x z || id: x, covariance(unstructured)"
    }
    else {
        local points 3
        local vce "vce(robust)"
        if "`cluster'" != "" local vce "vce(cluster `cluster')"
        local model "gsem (`first' <- x z I1[id]@1 c.x#S1[id]@1, family(poisson) link(log)) (`second' <- x z I2[id]@1 c.x#S2[id]@1, family(poisson) link(log)), cov(I1[id]*S1[id] I2[id]*S2[id] I1[id]*I2[id]@0 I1[id]*S2[id]@0 S1[id]*I2[id]@0 S1[id]*S2[id]@0) `vce'"
    }
    local ok 0
    local route "none"
    di as txt "--- COLD LAPLACE ATTEMPT ---"
    ereturn clear
    capture noisily `model' intmethod(laplace) iterate(200) log showtolerance
    local rc = _rc
    if `rc' == 1 exit 1
    di as result "COLD_LAPLACE_RC=" `rc'
    if `rc' == 0 {
        if e(converged) == 1 {
            local ok 1
            local route "cold"
        }
    }
    if !`ok' {
        di as txt "COLD ATTEMPT DIAGNOSTICS ONLY (not accepted benchmark estimates)"
        capture noisily di as txt "CONVERGED=" e(converged) " LL=" %21.15g e(ll)
        capture noisily matrix list e(b), format(%21.15g)
        capture noisily matrix list e(gradient), format(%21.15g)
        di as txt "--- ADAPTIVE INITIALIZER ONLY: `points' POINTS ---"
        ereturn clear
        capture noisily `model' intmethod(mvaghermite) intpoints(`points') iterate(100) log showtolerance
        local seed_rc = _rc
        if `seed_rc' == 1 exit 1
        di as result "ADAPTIVE_INITIALIZER_RC=" `seed_rc'
        local seed_ok 0
        if `seed_rc' == 0 {
            if e(converged) == 1 local seed_ok 1
        }
        di as result "ADAPTIVE_INITIALIZER_CONVERGED=" `seed_ok'
        capture noisily matrix list e(b), format(%21.15g)
        if `seed_ok' {
            tempname initial
            matrix `initial' = e(b)
            di as txt "--- WARM LAPLACE ATTEMPT ---"
            ereturn clear
            capture noisily `model' intmethod(laplace) from(`initial') iterate(200) log showtolerance
            local warm_rc = _rc
            if `warm_rc' == 1 exit 1
            di as result "WARM_LAPLACE_RC=" `warm_rc'
            if `warm_rc' == 0 {
                if e(converged) == 1 {
                    local ok 1
                    local route "adaptive_start"
                }
            }
            if !`ok' {
                di as txt "WARM ATTEMPT DIAGNOSTICS ONLY (not accepted benchmark estimates)"
                capture noisily di as txt "CONVERGED=" e(converged) " LL=" %21.15g e(ll)
                capture noisily matrix list e(b), format(%21.15g)
                capture noisily matrix list e(gradient), format(%21.15g)
            }
        }
    }
    di as result "FINAL_LAPLACE_ACCEPTED=" `ok' " START_ROUTE=`route'"
    * Do not expose an adaptive fit as the final model if a retry was skipped.
    if !`ok' ereturn clear
    ereturn scalar rs_success = `ok'
end

local component_count 0
local margin_count 0
foreach outcome in y1 y2 y1_partial y2_partial y1_left y2_right {
    di as txt "=== FIT ATTEMPTS: `outcome' ==="
    rs_fit, first(`outcome')
    if e(rs_success) == 1 {
        local ++component_count
        di as txt "=== LAPLACE COMPONENT: `outcome' ==="
        di as txt "N=" e(N) " LL=" %21.15g e(ll)
        matrix list e(N_g), format(%21.15g)
        matrix list e(b), format(%21.15g)
        matrix list e(V), format(%21.15g)
        if inlist("`outcome'", "y1", "y2") {
            di as txt "MARGINAL MEAN: `outcome'"
            capture noisily margins, predict(mu marginal)
            local mean_rc = _rc
            di as result "MARGINAL_MEAN_RC=" `mean_rc'
            if `mean_rc' == 1 exit 1
            if `mean_rc' == 0 {
                matrix list r(b), format(%21.15g)
                matrix list r(V), format(%21.15g)
                local mean_ok = rowsof(r(b)) == 1 & colsof(r(b)) == 1 & rowsof(r(V)) == 1 & colsof(r(V)) == 1
                if `mean_ok' {
                    local mean_ok = !missing(el(r(b), 1, 1), el(r(V), 1, 1)) & el(r(V), 1, 1) >= 0
                }
                di as result "MARGINAL_MEAN_ACCEPTED=" `mean_ok'
                if `mean_ok' local ++margin_count
            }
            di as txt "MARGINAL SLOPE X: `outcome'"
            capture noisily margins, dydx(x) predict(mu marginal)
            local slope_rc = _rc
            di as result "MARGINAL_SLOPE_RC=" `slope_rc'
            if `slope_rc' == 1 exit 1
            if `slope_rc' == 0 {
                matrix list r(b), format(%21.15g)
                matrix list r(V), format(%21.15g)
                local slope_ok = rowsof(r(b)) == 1 & colsof(r(b)) == 1 & rowsof(r(V)) == 1 & colsof(r(V)) == 1
                if `slope_ok' {
                    local slope_ok = !missing(el(r(b), 1, 1), el(r(V), 1, 1)) & el(r(V), 1, 1) >= 0
                }
                di as result "MARGINAL_SLOPE_ACCEPTED=" `slope_ok'
                if `slope_ok' local ++margin_count
            }
        }
    }
    else {
        di as error "COMPONENT_LAPLACE_FAILED: `outcome'; continuing to next case."
    }
}
if `component_count' == 6 di as result "POISSON_RS_COMPONENT_BENCHMARK_COMPLETE=1"

* Within each equation, intercept/slope covariance is free. Across equations,
* all four latent covariances are fixed to zero, exactly as in revision 1.
local joint_count 0
foreach case in balanced partial disjoint {
    local first y1
    local second y2
    local options ""
    if "`case'" == "partial" {
        local first y1_partial
        local second y2_partial
        local options "cluster(higher)"
    }
    if "`case'" == "disjoint" {
        local first y1_left
        local second y2_right
    }
    di as txt "=== JOINT FIT ATTEMPTS: `case' ==="
    rs_fit, first(`first') second(`second') `options'
    if e(rs_success) == 1 {
        local ++joint_count
        di as txt "=== JOINT GSEM LAPLACE: `case' ==="
        di as txt "N=" e(N) " N_clust=" e(N_clust) " LL=" %21.15g e(ll)
        matrix list e(b), format(%21.15g)
        matrix list e(V), format(%21.15g)
        di as result "GSEM_LAPLACE_`case'_COMPLETE=1"
    }
    else {
        di as error "JOINT_LAPLACE_FAILED: `case'; continuing to next case."
    }
}
di as result "CONVERGED_LAPLACE_COMPONENTS=" `component_count' "/6"
di as result "CONVERGED_LAPLACE_JOINT_SYSTEMS=" `joint_count' "/3"
di as result "SUCCESSFUL_NATIVE_MARGINS=" `margin_count' "/4"
if `component_count' == 6 & `joint_count' == 3 & `margin_count' == 4 {
    di as result "GLMMTMB_POISSON_RS_CROSSLANG_BENCHMARK_COMPLETE=1"
}
else {
    di as error "GLMMTMB_POISSON_RS_CROSSLANG_BENCHMARK_INCOMPLETE=1"
}
di as result "POISSON_RS_DIAGNOSTIC_RUN_FINISHED=1"
log close glmmrs
