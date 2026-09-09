version 16
clear all
set more off
set linesize 255
capture log close _all
log using suest_r_survey_binary_benchmark.log, replace text name(surveybinary)
di as result "BENCHMARK_FILE_REVISION=1"
di as txt "Stata version: `c(stata_version)'"
capture noisily which suest2
local available = (_rc == 0)
import delimited using suest_r_survey_binary_benchmark.csv, clear varnames(1) asdouble

* Six deterministic systems. The full design is retained; domains use subpop().
* Cases 1/3/5 use logit; cases 2/4/6 use probit.
forvalues case = 1/6 {
    capture drop dom_a dom_b
    generate byte dom_a = 1
    generate byte dom_b = 1
    local cmd "logit"
    local link "logit"
    local a_rhs "x"
    local b_rhs "x z"
    local a_opts ""
    local b_opts ""
    local stack_vars "a_cons a_x b_cons b_x b_z"
    local stack_opts ""

    if `case' == 2 {
        local cmd "probit"
        local link "probit"
    }
    if `case' == 3 {
        replace dom_a = unit <= 7
        replace dom_b = unit >= 4
        local a_rhs "x i.f"
        local b_rhs "x z i.f"
        local stack_vars "a_cons a_x a_f1 a_f2 b_cons b_x b_z b_f1 b_f2"
    }
    if `case' == 4 {
        replace dom_a = unit <= 5
        replace dom_b = unit > 5
        local cmd "probit"
        local link "probit"
        local a_opts ", offset(off_a)"
        local b_opts ", offset(off_b)"
        local stack_opts "offset(off_stack)"
    }
    if `case' == 5 {
        replace dom_a = mod(psu, 6) == 1
        replace dom_b = mod(psu, 6) == 2
        local a_rhs "x i.f"
        local b_rhs "x z i.f"
        local a_opts ", offset(off_a)"
        local b_opts ", offset(off_b)"
        local stack_vars "a_cons a_x a_f1 a_f2 b_cons b_x b_z b_f1 b_f2"
        local stack_opts "offset(off_stack)"
    }
    if `case' == 6 {
        replace dom_a = strata <= 2
        replace dom_b = strata > 2
        local cmd "probit"
        local link "probit"
        local a_opts ", offset(off_a)"
        local b_opts ", offset(off_b)"
        local stack_opts "offset(off_stack)"
    }

    if `case' == 1 {
        svyset psu [pweight = w], strata(strata) singleunit(missing)
    }
    else {
        svyset psu [pweight = w], strata(strata) fpc(pop) singleunit(missing)
    }

    di as result "SURVEY_BINARY_CASE=`case' LINK=`link'"
    quietly svy, subpop(dom_a): `cmd' yb `a_rhs' `a_opts'
    estimates store survey_a
    di as txt "NATIVE_MODEL=A DF=" e(df_r)
    matrix list e(b), format(%21.15g)
    matrix list e(V), format(%21.15g)

    * Same-sample cases deliberately use the same subpopulation expression.
    if `case' <= 2 {
        quietly svy, subpop(dom_a): `cmd' yb2 `b_rhs' `b_opts'
    }
    else {
        quietly svy, subpop(dom_b): `cmd' yb2 `b_rhs' `b_opts'
    }
    estimates store survey_b
    di as txt "NATIVE_MODEL=B DF=" e(df_r)
    matrix list e(b), format(%21.15g)
    matrix list e(V), format(%21.15g)

    if `available' {
        capture noisily suest2 survey_a survey_b
        local rc = _rc
        di as result "SUEST2_RC=`rc'"
        if `rc' == 0 {
            di as txt "SUEST2_JOINT_MODEL"
            matrix list e(b), format(%21.15g)
            matrix list e(V), format(%21.15g)
        }
        capture noisily suest2_cleanup, force
    }
    else {
        di as result "SUEST2_UNAVAILABLE=1"
    }

    * Independent full-joint reference. With one common link per case, the
    * expanded survey GLM has block-separable estimating equations but common
    * PSU/stratum linearization, yielding the complete cross-model covariance.
    preserve
    expand 2, generate(equation_b)
    generate byte y_stack = cond(equation_b == 0, yb, yb2)
    generate byte stack_subpop = cond(equation_b == 0, dom_a, dom_b)
    generate byte a_cons = equation_b == 0
    generate double a_x = x * a_cons
    generate byte b_cons = equation_b == 1
    generate double b_x = x * b_cons
    generate double b_z = z * b_cons
    if inlist(`case', 3, 5) {
        generate byte a_f1 = (f == 1) * a_cons
        generate byte a_f2 = (f == 2) * a_cons
        generate byte b_f1 = (f == 1) * b_cons
        generate byte b_f2 = (f == 2) * b_cons
    }
    if inlist(`case', 4, 5, 6) {
        generate double off_stack = cond(equation_b == 0, off_a, off_b)
    }
    quietly svy, subpop(stack_subpop): glm y_stack `stack_vars', ///
        family(binomial) link(`link') noconstant `stack_opts'
    di as txt "STACKED_JOINT_MODEL DF=" e(df_r)
    matrix list e(b), format(%21.15g)
    matrix list e(V), format(%21.15g)
    restore
}

di as result "SURVEY_BINARY_BENCHMARK_COMPLETE=1"
log close surveybinary
