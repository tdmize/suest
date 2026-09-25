version 16
clear all
set more off
set linesize 255
capture log close _all
log using suest_r_glmmtmb_logit_rs_disjoint_v1.log, replace text name(logitrs)
di as result "BENCHMARK_FILE_REVISION=1"
di as txt "LOGIT RANDOM SLOPES: DISJOINT ODD/EVEN GROUPS, 100 UNION CLUSTERS; THREE FITS ONLY"
di as txt "R seeds are checked at zero iterations; full and preliminary fixed-effects fits are capped at 15 iterations."
di as txt "Stop at the first failed check. No automatic restarts or additional fitting runs."
di as txt "Stata version: `c(stata_version)'"
import delimited using suest_r_glmmtmb_logit_rs_crosslang_benchmark.csv, clear varnames(1)
isid id time
assert _N == 1200
assert !missing(x, z, id, higher)
quietly count if !missing(y1_left)
assert r(N) == 600
quietly count if !missing(y2_right)
assert r(N) == 600
quietly count if !missing(y1_left, y2_right)
assert r(N) == 0
quietly count if !missing(y1_left) | !missing(y2_right)
assert r(N) == 1200
assert missing(y1_left) == (mod(id, 2) == 0)
assert missing(y2_right) == (mod(id, 2) == 1)
do suest_r_glmmtmb_logit_rs_disjoint_starts.do
matrix logitrs_components = J(2, 6, .)
matrix rownames logitrs_components = y1_left y2_right
matrix colnames logitrs_components = x z _cons var_slope var_intercept covariance
scalar logitrs_target_ll = 0

* Match Stata's own parameter stripes, allowing either covariance orientation.
mata:
real rowvector logitrs_indices(string matrix stripe, string matrix wanted) {
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
string matrix logitrs_component_stripes(string scalar outcome) {
    return((outcome, "x" \ outcome, "z" \ outcome, "_cons" \ ///
        "/", "var(x[id])" \ "/", "var(_cons[id])" \ "/", "cov(x[id],_cons[id])"))
}
void logitrs_component_start(string scalar outcome, real scalar row) {
    string matrix stripe
    real rowvector index, b
    real matrix seed
    stripe = st_matrixcolstripe("logitrs_start")
    index = logitrs_indices(stripe, logitrs_component_stripes(outcome))
    b = st_matrix("logitrs_start")
    if (cols(b) != 6) _error(3498)
    seed = st_matrix("logitrs_r_starts")
    b[1,index] = seed[row,]
    st_matrix("logitrs_start", b)
    st_matrixcolstripe("logitrs_start", stripe)
}
void logitrs_component_check() {
    real matrix v, ng
    v = st_matrix("e(V)"); ng = st_matrix("e(N_g)")
    st_numscalar("logitrs_ok", rows(v) == 6 & cols(v) == 6 & !hasmissing(v) & ng[1,1] == 50)
}
void logitrs_component_record(string scalar outcome, real scalar row) {
    real rowvector index, b
    real matrix components
    index = logitrs_indices(st_matrixcolstripe("e(b)"), logitrs_component_stripes(outcome))
    components = st_matrix("logitrs_components"); b = st_matrix("e(b)")
    components[row,] = b[1,index]
    st_matrix("logitrs_components", components)
}
void logitrs_margin_check() {
    real matrix b, v
    b = st_matrix("r(b)"); v = st_matrix("r(V)")
    st_numscalar("logitrs_margin_ok", rows(b) == 1 & cols(b) == 1 & !hasmissing(b) & ///
        rows(v) == 1 & cols(v) == 1 & !hasmissing(v) & v[1,1] >= 0)
}
end

local row 0
foreach outcome in y1_left y2_right {
    local ++row
    local expected_n = 600
    local model "melogit `outcome' x z || id: x, covariance(unstructured) intmethod(laplace)"
    di as txt "=== COMPONENT PARAMETER TEMPLATE: `outcome' ==="
    `model' noestimate startvalues(fixedonly, iterate(15))
    matrix logitrs_start = e(b)
    mata: logitrs_component_start("`outcome'", `row')
    matrix list logitrs_start, format(%21.15g)
    di as txt "=== COMPONENT ZERO-STEP CHECK: `outcome' ==="
    ereturn clear
    capture noisily `model' from(logitrs_start) startvalues(fixedonly, iterate(15)) iterate(0) log
    local eval_rc = _rc
    di as result "START_EVALUATION_RC=" `eval_rc'
    if !inlist(`eval_rc', 0, 430) {
        log close logitrs
        exit `eval_rc'
    }
    scalar logitrs_gap = e(ll)-scalar(logitrs_r_`outcome'_ll)
    mata: st_numscalar("logitrs_change", max(abs(st_matrix("e(b)")-st_matrix("logitrs_start"))))
    di as txt "INITIAL_LOGLIK_GAP=" %21.15g scalar(logitrs_gap)
    di as txt "MAX_START_PARAMETER_CHANGE=" %21.15g scalar(logitrs_change)
    if missing(scalar(logitrs_gap), scalar(logitrs_change)) | abs(scalar(logitrs_gap)) > 1e-4 | scalar(logitrs_change) > 1e-10 {
        di as error "COMPONENT_START_CHECK_FAILED: return this log."
        log close logitrs
        exit 459
    }
    di as txt "=== LAPLACE COMPONENT: `outcome' ==="
    ereturn clear
    capture noisily `model' from(logitrs_start) startvalues(fixedonly, iterate(15)) iterate(15) log showtolerance
    local fit_rc = _rc
    di as result "COMPONENT_FIT_RC=" `fit_rc'
    if `fit_rc' != 0 {
        di as error "COMPONENT_FIT_FAILED: return this log."
        log close logitrs
        exit `fit_rc'
    }
    di as txt "CONVERGED=" e(converged) " N=" e(N) " LL=" %21.15g e(ll)
    matrix list e(N_g), format(%21.15g)
    matrix list e(b), format(%21.15g)
    matrix list e(V), format(%21.15g)
    mata: logitrs_component_check()
    if e(converged) != 1 | e(N) != `expected_n' | missing(e(ll)) | abs(e(ll)-scalar(logitrs_r_`outcome'_ll)) > 1e-4 | scalar(logitrs_ok) != 1 {
        di as error "COMPONENT_CHECK_FAILED: return this log."
        log close logitrs
        exit 459
    }
    mata: logitrs_component_record("`outcome'", `row')
    scalar logitrs_target_ll = scalar(logitrs_target_ll)+e(ll)
    estimates store logitrs_`outcome'
    di as txt "MARGINAL MEAN: `outcome'"
    margins, predict(mu marginal intpoints(40))
    matrix list r(b), format(%21.15g)
    matrix list r(V), format(%21.15g)
    mata: logitrs_margin_check()
    if scalar(logitrs_margin_ok) != 1 {
        di as error "COMPONENT_MEAN_FAILED: return this log."
        log close logitrs
        exit 459
    }
    di as txt "MARGINAL SLOPE X: `outcome'"
    margins, dydx(x) predict(mu marginal intpoints(40))
    matrix list r(b), format(%21.15g)
    matrix list r(V), format(%21.15g)
    mata: logitrs_margin_check()
    if scalar(logitrs_margin_ok) != 1 {
        di as error "COMPONENT_SLOPE_FAILED: return this log."
        log close logitrs
        exit 459
    }
    di as result "LOGIT_RS_COMPONENT_`outcome'_COMPLETE=1"
}

* Start the joint fit from actual Stata component estimates. Both within-
* equation intercept/slope covariances are free; all cross-equation latent
* covariances are fixed at zero. The robust VCE retains outcome dependence.
local model "gsem (y1_left <- x z I1[id]@1 c.x#S1[id]@1, family(bernoulli) link(logit)) (y2_right <- x z I2[id]@1 c.x#S2[id]@1, family(bernoulli) link(logit)), cov(I1[id]*S1[id] I2[id]*S2[id] I1[id]*I2[id]@0 I1[id]*S2[id]@0 S1[id]*I2[id]@0 S1[id]*S2[id]@0) intmethod(laplace) adaptopts(tolerance(1e-12)) vce(robust)"
di as txt "=== JOINT PARAMETER TEMPLATE (NO ESTIMATION) ==="
`model' noestimate startvalues(zero)
matrix logitrs_start = e(b)
mata:
logitrs_s = st_matrixcolstripe("logitrs_start")
logitrs_w = ("y1_left", "x" \ "y1_left", "z" \ "y1_left", "_cons" \ ///
    "/", "var(S1[id])" \ "/", "var(I1[id])" \ "/", "cov(I1[id],S1[id])" \ ///
    "y2_right", "x" \ "y2_right", "z" \ "y2_right", "_cons" \ ///
    "/", "var(S2[id])" \ "/", "var(I2[id])" \ "/", "cov(I2[id],S2[id])")
logitrs_ix = logitrs_indices(logitrs_s, logitrs_w)
logitrs_c = st_matrix("logitrs_components")
logitrs_b = st_matrix("logitrs_start")
if (cols(logitrs_b) != 16) _error(3498)
logitrs_b[1,logitrs_ix] = (logitrs_c[1,], logitrs_c[2,])
st_matrix("logitrs_start", logitrs_b)
st_matrixcolstripe("logitrs_start", logitrs_s)
end
di as txt "EXPECTED_JOINT_LOGLIK=" %21.15g scalar(logitrs_target_ll)
matrix list logitrs_start, format(%21.15g)
di as txt "=== JOINT AT COMPONENT ESTIMATES (ZERO ITERATIONS) ==="
ereturn clear
capture noisily `model' from(logitrs_start) startvalues(zero) iterate(0) log
local eval_rc = _rc
di as result "START_EVALUATION_RC=" `eval_rc'
if !inlist(`eval_rc', 0, 430) {
    log close logitrs
    exit `eval_rc'
}
scalar logitrs_gap = e(ll)-scalar(logitrs_target_ll)
mata: st_numscalar("logitrs_change", max(abs(st_matrix("e(b)")-st_matrix("logitrs_start"))))
di as txt "INITIAL_LOGLIK_GAP=" %21.15g scalar(logitrs_gap)
di as txt "MAX_START_PARAMETER_CHANGE=" %21.15g scalar(logitrs_change)
if missing(scalar(logitrs_gap), scalar(logitrs_change)) | abs(scalar(logitrs_gap)) > 1e-4 | scalar(logitrs_change) > 1e-10 {
    di as error "JOINT_START_CHECK_FAILED: return this log."
    log close logitrs
    exit 459
}
di as txt "=== DISJOINT JOINT GSEM LAPLACE (MAXIMUM 15 ITERATIONS) ==="
ereturn clear
capture noisily `model' from(logitrs_start) startvalues(zero) iterate(15) log showtolerance
local fit_rc = _rc
di as result "DISJOINT_FIT_RC=" `fit_rc'
if `fit_rc' != 0 {
    di as error "JOINT_FIT_FAILED: return this log."
    log close logitrs
    exit `fit_rc'
}
di as txt "CONVERGED=" e(converged) " N=" e(N) " N_clust=" e(N_clust) " LL=" %21.15g e(ll)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
matrix list e(V_modelbased), format(%21.15g)
capture noisily matrix list e(gradient), format(%21.15g)
mata:
logitrs_v = st_matrix("e(V)"); logitrs_vm = st_matrix("e(V_modelbased)")
st_numscalar("logitrs_ok", rows(logitrs_v) == cols(logitrs_b) & cols(logitrs_v) == cols(logitrs_b) & !hasmissing(logitrs_v) & ///
    rows(logitrs_vm) == cols(logitrs_b) & cols(logitrs_vm) == cols(logitrs_b) & !hasmissing(logitrs_vm))
end
if e(converged) != 1 | e(N) != 1200 | e(N_clust) != 100 | missing(e(ll)) | abs(e(ll)-scalar(logitrs_target_ll)) > 1e-4 | scalar(logitrs_ok) != 1 {
    di as error "JOINT_CHECK_FAILED: return this log."
    log close logitrs
    exit 459
}
di as result "LOGIT_RS_DISJOINT_COMPLETE=1"
di as txt "Completion confirms successful fits; cross-language agreement will be checked from this log."
log close logitrs
