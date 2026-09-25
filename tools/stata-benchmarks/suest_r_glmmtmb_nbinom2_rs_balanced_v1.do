version 16
clear all
set more off
set linesize 255
capture log close _all
log using suest_r_glmmtmb_nbinom2_rs_balanced_v1.log, replace text name(nb2rs)
di as result "BENCHMARK_FILE_REVISION=1"
di as txt "NB2 RANDOM SLOPES: TWO BALANCED COMPONENTS AND ONE JOINT FIT ONLY"
di as txt "R seeds are checked at zero iterations; every fit is capped at 15 iterations."
di as txt "Stop at the first failed check. No automatic restarts or quadrature runs."
di as txt "Stata version: `c(stata_version)'"
import delimited using suest_r_glmmtmb_nbinom2_rs_crosslang_benchmark.csv, clear varnames(1)
isid id time
assert _N == 1000
assert !missing(y1, y2, x, z, id)
do suest_r_glmmtmb_nbinom2_rs_start_values.do
matrix nb2rs_components = J(2, 7, .)
matrix rownames nb2rs_components = y1 y2
matrix colnames nb2rs_components = x z _cons lnalpha var_slope var_intercept covariance
scalar nb2rs_target_ll = 0

* Match Stata's own parameter stripes, allowing either covariance orientation.
mata:
real rowvector nb2rs_indices(string matrix stripe, string matrix wanted) {
    real scalar i
    real rowvector answer, hit
    string scalar reverse
    answer = J(1, rows(wanted), .)
    for (i = 1; i <= rows(wanted); i++) {
        hit = selectindex((stripe[,1] :== wanted[i,1]) :& (stripe[,2] :== wanted[i,2]))'
        if (length(hit) == 0) {
            reverse = ""
            if (wanted[i,2] == "cov(x[id],_cons[id])") reverse = "cov(_cons[id],x[id])"
            if (wanted[i,2] == "cov(I1[id],S1[id])") reverse = "cov(S1[id],I1[id])"
            if (wanted[i,2] == "cov(I2[id],S2[id])") reverse = "cov(S2[id],I2[id])"
            if (reverse != "") hit = selectindex((stripe[,1] :== wanted[i,1]) :& (stripe[,2] :== reverse))'
        }
        if (length(hit) != 1) {
            errprintf("Unmatched parameter: %s:%s\n", wanted[i,1], wanted[i,2])
            _error(3498)
        }
        answer[i] = hit[1]
    }
    return(answer)
}
string matrix nb2rs_component_stripes(string scalar outcome) {
    return((outcome, "x" \ outcome, "z" \ outcome, "_cons" \ ///
        "/", "lnalpha" \ "/", "var(x[id])" \ "/", "var(_cons[id])" \ "/", "cov(x[id],_cons[id])"))
}
void nb2rs_component_start(string scalar outcome, real scalar row) {
    string matrix stripe
    real rowvector index, b
    real matrix seed
    stripe = st_matrixcolstripe("nb2rs_start")
    index = nb2rs_indices(stripe, nb2rs_component_stripes(outcome))
    b = st_matrix("nb2rs_start")
    if (cols(b) != 7) _error(3498)
    seed = st_matrix("nb2rs_r_starts")
    b[1,index] = seed[row,]
    st_matrix("nb2rs_start", b)
    st_matrixcolstripe("nb2rs_start", stripe)
}
void nb2rs_component_check() {
    real matrix v, ng
    v = st_matrix("e(V)"); ng = st_matrix("e(N_g)")
    st_numscalar("nb2rs_ok", rows(v) == 7 & cols(v) == 7 & !hasmissing(v) & ng[1,1] == 100)
}
void nb2rs_component_record(string scalar outcome, real scalar row) {
    real rowvector index, b
    real matrix components
    index = nb2rs_indices(st_matrixcolstripe("e(b)"), nb2rs_component_stripes(outcome))
    components = st_matrix("nb2rs_components"); b = st_matrix("e(b)")
    components[row,] = b[1,index]
    st_matrix("nb2rs_components", components)
}
void nb2rs_margin_check() {
    real matrix b, v
    b = st_matrix("r(b)"); v = st_matrix("r(V)")
    st_numscalar("nb2rs_margin_ok", rows(b) == 1 & cols(b) == 1 & !hasmissing(b) & ///
        rows(v) == 1 & cols(v) == 1 & !hasmissing(v) & v[1,1] >= 0)
}
end

local row 0
foreach outcome in y1 y2 {
    local ++row
    local model "menbreg `outcome' x z || id: x, covariance(unstructured) dispersion(mean) intmethod(laplace)"
    di as txt "=== COMPONENT PARAMETER TEMPLATE: `outcome' ==="
    `model' noestimate startvalues(zero)
    matrix nb2rs_start = e(b)
    mata: nb2rs_component_start("`outcome'", `row')
    matrix list nb2rs_start, format(%21.15g)
    di as txt "=== COMPONENT ZERO-STEP CHECK: `outcome' ==="
    ereturn clear
    capture noisily `model' from(nb2rs_start) startvalues(zero) iterate(0) log
    local eval_rc = _rc
    di as result "START_EVALUATION_RC=" `eval_rc'
    if !inlist(`eval_rc', 0, 430) {
        log close nb2rs
        exit `eval_rc'
    }
    scalar nb2rs_gap = e(ll)-scalar(nb2rs_r_`outcome'_ll)
    mata: st_numscalar("nb2rs_change", max(abs(st_matrix("e(b)")-st_matrix("nb2rs_start"))))
    di as txt "INITIAL_LOGLIK_GAP=" %21.15g scalar(nb2rs_gap)
    di as txt "MAX_START_PARAMETER_CHANGE=" %21.15g scalar(nb2rs_change)
    if missing(scalar(nb2rs_gap), scalar(nb2rs_change)) | abs(scalar(nb2rs_gap)) > 1e-4 | scalar(nb2rs_change) > 1e-10 {
        di as error "COMPONENT_START_CHECK_FAILED: return this log."
        log close nb2rs
        exit 459
    }
    di as txt "=== LAPLACE COMPONENT: `outcome' ==="
    ereturn clear
    capture noisily `model' from(nb2rs_start) startvalues(zero) iterate(15) log showtolerance
    local fit_rc = _rc
    di as result "COMPONENT_FIT_RC=" `fit_rc'
    if `fit_rc' != 0 {
        di as error "COMPONENT_FIT_FAILED: return this log."
        log close nb2rs
        exit `fit_rc'
    }
    di as txt "CONVERGED=" e(converged) " N=" e(N) " LL=" %21.15g e(ll)
    matrix list e(N_g), format(%21.15g)
    matrix list e(b), format(%21.15g)
    matrix list e(V), format(%21.15g)
    mata: nb2rs_component_check()
    if e(converged) != 1 | e(N) != 1000 | missing(e(ll)) | abs(e(ll)-scalar(nb2rs_r_`outcome'_ll)) > 1e-4 | scalar(nb2rs_ok) != 1 {
        di as error "COMPONENT_CHECK_FAILED: return this log."
        log close nb2rs
        exit 459
    }
    mata: nb2rs_component_record("`outcome'", `row')
    scalar nb2rs_target_ll = scalar(nb2rs_target_ll)+e(ll)
    estimates store nb2rs_`outcome'
    di as txt "MARGINAL MEAN: `outcome'"
    margins, predict(mu marginal)
    matrix list r(b), format(%21.15g)
    matrix list r(V), format(%21.15g)
    mata: nb2rs_margin_check()
    if scalar(nb2rs_margin_ok) != 1 {
        di as error "COMPONENT_MEAN_FAILED: return this log."
        log close nb2rs
        exit 459
    }
    di as txt "MARGINAL SLOPE X: `outcome'"
    margins, dydx(x) predict(mu marginal)
    matrix list r(b), format(%21.15g)
    matrix list r(V), format(%21.15g)
    mata: nb2rs_margin_check()
    if scalar(nb2rs_margin_ok) != 1 {
        di as error "COMPONENT_SLOPE_FAILED: return this log."
        log close nb2rs
        exit 459
    }
    di as result "NB2_RS_COMPONENT_`outcome'_COMPLETE=1"
}

* Start the joint fit from actual Stata component estimates. Both within-
* equation intercept/slope covariances are free; all cross-equation latent
* covariances are fixed at zero. The robust VCE retains outcome dependence.
local model "gsem (y1 <- x z I1[id]@1 c.x#S1[id]@1, family(nbinomial mean) link(log)) (y2 <- x z I2[id]@1 c.x#S2[id]@1, family(nbinomial mean) link(log)), cov(I1[id]*S1[id] I2[id]*S2[id] I1[id]*I2[id]@0 I1[id]*S2[id]@0 S1[id]*I2[id]@0 S1[id]*S2[id]@0) intmethod(laplace) adaptopts(tolerance(1e-12)) vce(robust)"
di as txt "=== JOINT PARAMETER TEMPLATE (NO ESTIMATION) ==="
`model' noestimate startvalues(zero)
matrix nb2rs_start = e(b)
mata:
nb2rs_s = st_matrixcolstripe("nb2rs_start")
nb2rs_w = ("y1", "x" \ "y1", "z" \ "y1", "_cons" \ "/y1", "lnalpha" \ ///
    "/", "var(S1[id])" \ "/", "var(I1[id])" \ "/", "cov(I1[id],S1[id])" \ ///
    "y2", "x" \ "y2", "z" \ "y2", "_cons" \ "/y2", "lnalpha" \ ///
    "/", "var(S2[id])" \ "/", "var(I2[id])" \ "/", "cov(I2[id],S2[id])")
nb2rs_ix = nb2rs_indices(nb2rs_s, nb2rs_w)
nb2rs_c = st_matrix("nb2rs_components")
nb2rs_b = st_matrix("nb2rs_start")
nb2rs_b[1,nb2rs_ix] = (nb2rs_c[1,], nb2rs_c[2,])
st_matrix("nb2rs_start", nb2rs_b)
st_matrixcolstripe("nb2rs_start", nb2rs_s)
end
di as txt "EXPECTED_JOINT_LOGLIK=" %21.15g scalar(nb2rs_target_ll)
matrix list nb2rs_start, format(%21.15g)
di as txt "=== JOINT AT COMPONENT ESTIMATES (ZERO ITERATIONS) ==="
ereturn clear
capture noisily `model' from(nb2rs_start) startvalues(zero) iterate(0) log
local eval_rc = _rc
di as result "START_EVALUATION_RC=" `eval_rc'
if !inlist(`eval_rc', 0, 430) {
    log close nb2rs
    exit `eval_rc'
}
scalar nb2rs_gap = e(ll)-scalar(nb2rs_target_ll)
mata: st_numscalar("nb2rs_change", max(abs(st_matrix("e(b)")-st_matrix("nb2rs_start"))))
di as txt "INITIAL_LOGLIK_GAP=" %21.15g scalar(nb2rs_gap)
di as txt "MAX_START_PARAMETER_CHANGE=" %21.15g scalar(nb2rs_change)
if missing(scalar(nb2rs_gap), scalar(nb2rs_change)) | abs(scalar(nb2rs_gap)) > 1e-4 | scalar(nb2rs_change) > 1e-10 {
    di as error "JOINT_START_CHECK_FAILED: return this log."
    log close nb2rs
    exit 459
}
di as txt "=== BALANCED JOINT GSEM LAPLACE (MAXIMUM 15 ITERATIONS) ==="
ereturn clear
capture noisily `model' from(nb2rs_start) startvalues(zero) iterate(15) log showtolerance
local fit_rc = _rc
di as result "BALANCED_FIT_RC=" `fit_rc'
if `fit_rc' != 0 {
    di as error "JOINT_FIT_FAILED: return this log."
    log close nb2rs
    exit `fit_rc'
}
di as txt "CONVERGED=" e(converged) " N=" e(N) " N_clust=" e(N_clust) " LL=" %21.15g e(ll)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
matrix list e(V_modelbased), format(%21.15g)
capture noisily matrix list e(gradient), format(%21.15g)
mata:
nb2rs_v = st_matrix("e(V)"); nb2rs_vm = st_matrix("e(V_modelbased)")
st_numscalar("nb2rs_ok", rows(nb2rs_v) == cols(nb2rs_b) & cols(nb2rs_v) == cols(nb2rs_b) & !hasmissing(nb2rs_v) & ///
    rows(nb2rs_vm) == cols(nb2rs_b) & cols(nb2rs_vm) == cols(nb2rs_b) & !hasmissing(nb2rs_vm))
end
if e(converged) != 1 | e(N) != 1000 | e(N_clust) != 100 | missing(e(ll)) | abs(e(ll)-scalar(nb2rs_target_ll)) > 1e-4 | scalar(nb2rs_ok) != 1 {
    di as error "JOINT_CHECK_FAILED: return this log."
    log close nb2rs
    exit 459
}
di as result "NB2_RS_BALANCED_COMPLETE=1"
di as txt "Completion confirms successful fits; cross-language agreement will be checked from this log."
log close nb2rs
