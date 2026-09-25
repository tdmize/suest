version 16
clear all
set more off
set linesize 255
capture log close _all
log using suest_r_glmmtmb_poisson_rs_balanced_v3.log, replace text name(glmmrs)
di as result "BENCHMARK_FILE_REVISION=3"
di as txt "BALANCED JOINT MODEL ONLY; RETURNED STATA COMPONENT STARTS"
di as txt "Stata version: `c(stata_version)'"
import delimited using suest_r_glmmtmb_poisson_rs_crosslang_benchmark.csv, clear varnames(1)
isid id time
assert _N == 800
assert !missing(y1, y2, x, z, id)

* Exact successful Laplace component estimates from the returned v2 log.
* Columns: beta_x, beta_z, intercept, var_slope, var_intercept, covariance.
matrix rs_components = (.268437223537171, -.172164480056694, .581914673223932, .176558869264552, .245752485470852, .0908034866545292 \ ///
    .234543117800827, .136938795446797, .419193228426063, .152194025941785, .256778047093543, .0220419636638305)
matrix rownames rs_components = y1 y2
matrix colnames rs_components = x z _cons var_slope var_intercept covariance
scalar rs_target_ll = -1440.49068695319 -1371.93813135237
di as txt "EXPECTED_JOINT_LOGLIK=" %21.15g scalar(rs_target_ll)
matrix list rs_components, format(%21.15g)

* Same joint likelihood and robust VCE as v2. No component re-estimation.
local model "gsem (y1 <- x z I1[id]@1 c.x#S1[id]@1, family(poisson) link(log)) (y2 <- x z I2[id]@1 c.x#S2[id]@1, family(poisson) link(log)), cov(I1[id]*S1[id] I2[id]*S2[id] I1[id]*I2[id]@0 I1[id]*S2[id]@0 S1[id]*I2[id]@0 S1[id]*S2[id]@0) intmethod(laplace) vce(robust)"

* Obtain Stata's own parameter stripes, then fill all 12 free parameters by
* name. A missing/ambiguous match stops here, before any optimization.
di as txt "=== JOINT PARAMETER TEMPLATE (NO ESTIMATION) ==="
`model' noestimate
matrix rs_start = e(b)
mata:
rs_b = st_matrix("rs_start")
rs_s = st_matrixcolstripe("rs_start")
rs_c = st_matrix("rs_components")
rs_wanted = ("y1", "x" \ "y1", "z" \ "y1", "_cons" \ ///
    "/", "var(S1[id])" \ "/", "var(I1[id])" \ "/", "cov(I1[id],S1[id])" \ ///
    "y2", "x" \ "y2", "z" \ "y2", "_cons" \ ///
    "/", "var(S2[id])" \ "/", "var(I2[id])" \ "/", "cov(I2[id],S2[id])")
rs_values = (rs_c[1,], rs_c[2,])
for (rs_i = 1; rs_i <= 12; rs_i++) {
    rs_hit = selectindex((rs_s[,1] :== rs_wanted[rs_i,1]) :& (rs_s[,2] :== rs_wanted[rs_i,2]))
    if (length(rs_hit) == 0 & (rs_i == 6 | rs_i == 12)) {
        rs_reverse = rs_i == 6 ? "cov(S1[id],I1[id])" : "cov(S2[id],I2[id])"
        rs_hit = selectindex((rs_s[,1] :== "/") :& (rs_s[,2] :== rs_reverse))
    }
    if (length(rs_hit) != 1) {
        errprintf("Unmatched parameter: %s:%s\n", rs_wanted[rs_i,1], rs_wanted[rs_i,2])
        _error(3498)
    }
    rs_b[1,rs_hit[1]] = rs_values[rs_i]
}
st_matrix("rs_start", rs_b)
st_matrixcolstripe("rs_start", rs_s)
end
matrix list rs_start, format(%21.15g)

* Evaluate once without optimizing. Agreement with the sum of component
* likelihoods checks the mapping/specification before a longer fit is allowed.
di as txt "=== JOINT AT COMPONENT ESTIMATES (ZERO ITERATIONS) ==="
ereturn clear
capture noisily `model' from(rs_start) startvalues(zero) iterate(0) log
local eval_rc = _rc
di as result "START_EVALUATION_RC=" `eval_rc'
if !inlist(`eval_rc', 0, 430) {
    log close glmmrs
    exit `eval_rc'
}
scalar rs_initial_ll = e(ll)
di as txt "INITIAL_JOINT_LOGLIK=" %21.15g scalar(rs_initial_ll)
di as txt "INITIAL_LOGLIK_GAP=" %21.15g (scalar(rs_initial_ll)-scalar(rs_target_ll))
capture noisily matrix list e(b), format(%21.15g)
capture noisily matrix list e(gradient), format(%21.15g)
if missing(scalar(rs_initial_ll)) | abs(scalar(rs_initial_ll)-scalar(rs_target_ll)) > 1e-4 {
    di as error "START_LIKELIHOOD_MISMATCH: stopped before optimization. Return this log."
    log close glmmrs
    exit 459
}
* A zero-step evaluation must preserve the complete supplied parameter vector.
mata: st_numscalar("rs_start_change", max(abs(st_matrix("e(b)")-st_matrix("rs_start"))))
di as txt "MAX_START_PARAMETER_CHANGE=" %21.15g scalar(rs_start_change)
if missing(scalar(rs_start_change)) | scalar(rs_start_change) > 1e-10 {
    di as error "START_PARAMETERS_CHANGED: stopped before optimization. Return this log."
    log close glmmrs
    exit 459
}

di as txt "=== BALANCED JOINT LAPLACE FROM COMPONENT FITS (MAXIMUM 15 ITERATIONS) ==="
ereturn clear
capture noisily `model' from(rs_start) startvalues(zero) iterate(15) log showtolerance
local fit_rc = _rc
di as result "BALANCED_FIT_RC=" `fit_rc'
capture noisily di as txt "CONVERGED=" e(converged) " N=" e(N) " N_clust=" e(N_clust) " LL=" %21.15g e(ll)
capture noisily matrix list e(b), format(%21.15g)
capture noisily matrix list e(V), format(%21.15g)
capture noisily matrix list e(gradient), format(%21.15g)
local ok 0
scalar rs_vcov_ok = 0
if `fit_rc' == 0 {
    if e(converged) == 1 & e(N) == 800 & e(N_clust) == 100 & abs(e(ll)-scalar(rs_target_ll)) < 1e-4 local ok 1
    if `ok' {
        mata: rs_v = st_matrix("e(V)"); st_numscalar("rs_vcov_ok", rows(rs_v) == cols(rs_b) & cols(rs_v) == cols(rs_b) & !hasmissing(rs_v))
        local ok = scalar(rs_vcov_ok)
    }
}
di as result "FINITE_JOINT_COVARIANCE=" scalar(rs_vcov_ok)
if `ok' {
    di as result "POISSON_RS_BALANCED_WARMSTART_COMPLETE=1"
}
else {
    di as error "POISSON_RS_BALANCED_WARMSTART_INCOMPLETE=1"
}
di as result "POISSON_RS_FOCUSED_RUN_FINISHED=1"
log close glmmrs
