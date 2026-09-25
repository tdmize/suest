version 16
clear all
set more off
set linesize 255
capture log close _all
log using suest_r_glmmtmb_poisson_rs_remaining_v5.log, replace text name(glmmrs)
di as result "BENCHMARK_FILE_REVISION=5"
di as txt "PARTIAL AND DISJOINT JOINT MODELS; SUCCESSFUL STATA COMPONENT STARTS"
di as txt "Same Laplace likelihood and tighter internal tolerance as balanced v4."
di as txt "Maximum 15 optimization iterations per fit; stop at the first failure."
di as txt "Stata version: `c(stata_version)'"
import delimited using suest_r_glmmtmb_poisson_rs_crosslang_benchmark.csv, clear varnames(1)
isid id time
assert _N == 800
assert !missing(x, z, id, higher)
* Verify the supplied sample patterns before fitting either joint model.
quietly count if !missing(y1_partial)
assert r(N) == 775
quietly count if !missing(y2_partial)
assert r(N) == 780
quietly count if !missing(y1_partial, y2_partial)
assert r(N) == 755
quietly count if !missing(y1_left)
assert r(N) == 400
quietly count if !missing(y2_right)
assert r(N) == 400
quietly count if !missing(y1_left, y2_right)
assert r(N) == 0
bysort id (higher): assert higher[1] == higher[_N]

di as txt "=== SAMPLE CASE: partial ==="
* Exact successful Laplace component estimates from the returned v2 log.
* Columns: beta_x, beta_z, intercept, var_slope, var_intercept, covariance.
matrix rs_components = (.271284848596774, -.184986424866489, .567619713906397, .169394546915135, .257748081547828, .0898870925491518 \ ///
    .228126752202568, .136553411904341, .422160987385958, .153560362896159, .248196290451108, .0259134654180823)
matrix rownames rs_components = y1_partial y2_partial
matrix colnames rs_components = x z _cons var_slope var_intercept covariance
scalar rs_target_ll = -1389.75007946091 -1338.53800664237
di as txt "EXPECTED_JOINT_LOGLIK=" %21.15g scalar(rs_target_ll)
matrix list rs_components, format(%21.15g)

* Same joint likelihood and robust VCE as v2. No component re-estimation.
local model "gsem (y1_partial <- x z I1[id]@1 c.x#S1[id]@1, family(poisson) link(log)) (y2_partial <- x z I2[id]@1 c.x#S2[id]@1, family(poisson) link(log)), cov(I1[id]*S1[id] I2[id]*S2[id] I1[id]*I2[id]@0 I1[id]*S2[id]@0 S1[id]*I2[id]@0 S1[id]*S2[id]@0) intmethod(laplace) adaptopts(tolerance(1e-12)) vce(cluster higher)"

* Obtain Stata's own parameter stripes, then fill all 12 free parameters by
* name. A missing/ambiguous match stops here, before any optimization.
di as txt "=== JOINT PARAMETER TEMPLATE (NO ESTIMATION) ==="
`model' noestimate
matrix rs_start = e(b)
mata:
rs_b = st_matrix("rs_start")
rs_s = st_matrixcolstripe("rs_start")
rs_c = st_matrix("rs_components")
rs_wanted = ("y1_partial", "x" \ "y1_partial", "z" \ "y1_partial", "_cons" \ ///
    "/", "var(S1[id])" \ "/", "var(I1[id])" \ "/", "cov(I1[id],S1[id])" \ ///
    "y2_partial", "x" \ "y2_partial", "z" \ "y2_partial", "_cons" \ ///
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

di as txt "=== PARTIAL JOINT LAPLACE FROM COMPONENT FITS (MAXIMUM 15 ITERATIONS) ==="
ereturn clear
capture noisily `model' from(rs_start) startvalues(zero) iterate(15) log showtolerance
local fit_rc = _rc
if `fit_rc' == 1 {
    log close glmmrs
    exit 1
}
di as result "PARTIAL_FIT_RC=" `fit_rc'
capture noisily di as txt "CONVERGED=" e(converged) " N=" e(N) " N_clust=" e(N_clust) " LL=" %21.15g e(ll)
capture noisily matrix list e(b), format(%21.15g)
capture noisily matrix list e(V), format(%21.15g)
capture noisily matrix list e(V_modelbased), format(%21.15g)
local modelbased_rc = _rc
di as result "MODELBASED_COVARIANCE_RC=" `modelbased_rc'
capture noisily matrix list e(gradient), format(%21.15g)
local ok 0
scalar rs_vcov_ok = 0
if `fit_rc' == 0 {
    if e(converged) == 1 & e(N) == 800 & e(N_clust) == 20 & abs(e(ll)-scalar(rs_target_ll)) < 1e-4 local ok 1
    if `ok' {
        mata: rs_v = st_matrix("e(V)"); st_numscalar("rs_vcov_ok", rows(rs_v) == cols(rs_b) & cols(rs_v) == cols(rs_b) & !hasmissing(rs_v))
        local ok = scalar(rs_vcov_ok) & `modelbased_rc' == 0
        if `ok' {
            mata: rs_vm = st_matrix("e(V_modelbased)"); st_numscalar("rs_modelbased_ok", rows(rs_vm) == cols(rs_b) & cols(rs_vm) == cols(rs_b) & !hasmissing(rs_vm))
            local ok = scalar(rs_modelbased_ok)
        }
    }
}
di as result "FINITE_JOINT_COVARIANCE=" scalar(rs_vcov_ok)
if `ok' {
    di as result "POISSON_RS_PARTIAL_PRECISION_COMPLETE=1"
}
else {
    di as error "POISSON_RS_PARTIAL_PRECISION_INCOMPLETE=1"
    log close glmmrs
    exit 459
}

di as txt "=== SAMPLE CASE: disjoint ==="
* Exact successful Laplace component estimates from the returned v2 log.
* Columns: beta_x, beta_z, intercept, var_slope, var_intercept, covariance.
matrix rs_components = (.281949035745495, -.217119752503643, .509650418862606, .118231574542104, .371005189235335, .109728146698855 \ ///
    .271810868264066, .178733355393067, .395582030044842, .140694744398986, .194020980686681, .023551965266474)
matrix rownames rs_components = y1_left y2_right
matrix colnames rs_components = x z _cons var_slope var_intercept covariance
scalar rs_target_ll = -709.596417255838 -673.491682963962
di as txt "EXPECTED_JOINT_LOGLIK=" %21.15g scalar(rs_target_ll)
matrix list rs_components, format(%21.15g)

* Same joint likelihood and robust VCE as v2. No component re-estimation.
local model "gsem (y1_left <- x z I1[id]@1 c.x#S1[id]@1, family(poisson) link(log)) (y2_right <- x z I2[id]@1 c.x#S2[id]@1, family(poisson) link(log)), cov(I1[id]*S1[id] I2[id]*S2[id] I1[id]*I2[id]@0 I1[id]*S2[id]@0 S1[id]*I2[id]@0 S1[id]*S2[id]@0) intmethod(laplace) adaptopts(tolerance(1e-12)) vce(robust)"

* Obtain Stata's own parameter stripes, then fill all 12 free parameters by
* name. A missing/ambiguous match stops here, before any optimization.
di as txt "=== JOINT PARAMETER TEMPLATE (NO ESTIMATION) ==="
`model' noestimate
matrix rs_start = e(b)
mata:
rs_b = st_matrix("rs_start")
rs_s = st_matrixcolstripe("rs_start")
rs_c = st_matrix("rs_components")
rs_wanted = ("y1_left", "x" \ "y1_left", "z" \ "y1_left", "_cons" \ ///
    "/", "var(S1[id])" \ "/", "var(I1[id])" \ "/", "cov(I1[id],S1[id])" \ ///
    "y2_right", "x" \ "y2_right", "z" \ "y2_right", "_cons" \ ///
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

di as txt "=== DISJOINT JOINT LAPLACE FROM COMPONENT FITS (MAXIMUM 15 ITERATIONS) ==="
ereturn clear
capture noisily `model' from(rs_start) startvalues(zero) iterate(15) log showtolerance
local fit_rc = _rc
if `fit_rc' == 1 {
    log close glmmrs
    exit 1
}
di as result "DISJOINT_FIT_RC=" `fit_rc'
capture noisily di as txt "CONVERGED=" e(converged) " N=" e(N) " N_clust=" e(N_clust) " LL=" %21.15g e(ll)
capture noisily matrix list e(b), format(%21.15g)
capture noisily matrix list e(V), format(%21.15g)
capture noisily matrix list e(V_modelbased), format(%21.15g)
local modelbased_rc = _rc
di as result "MODELBASED_COVARIANCE_RC=" `modelbased_rc'
capture noisily matrix list e(gradient), format(%21.15g)
local ok 0
scalar rs_vcov_ok = 0
if `fit_rc' == 0 {
    if e(converged) == 1 & e(N) == 800 & e(N_clust) == 100 & abs(e(ll)-scalar(rs_target_ll)) < 1e-4 local ok 1
    if `ok' {
        mata: rs_v = st_matrix("e(V)"); st_numscalar("rs_vcov_ok", rows(rs_v) == cols(rs_b) & cols(rs_v) == cols(rs_b) & !hasmissing(rs_v))
        local ok = scalar(rs_vcov_ok) & `modelbased_rc' == 0
        if `ok' {
            mata: rs_vm = st_matrix("e(V_modelbased)"); st_numscalar("rs_modelbased_ok", rows(rs_vm) == cols(rs_b) & cols(rs_vm) == cols(rs_b) & !hasmissing(rs_vm))
            local ok = scalar(rs_modelbased_ok)
        }
    }
}
di as result "FINITE_JOINT_COVARIANCE=" scalar(rs_vcov_ok)
if `ok' {
    di as result "POISSON_RS_DISJOINT_PRECISION_COMPLETE=1"
}
else {
    di as error "POISSON_RS_DISJOINT_PRECISION_INCOMPLETE=1"
    log close glmmrs
    exit 459
}
di as result "POISSON_RS_REMAINING_JOINT_MODELS_COMPLETE=1"
di as result "POISSON_RS_FOCUSED_RUN_FINISHED=1"
log close glmmrs
