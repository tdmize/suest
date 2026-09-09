version 16
clear all
set more off
set linesize 255
capture log close _all
log using suest_r_survey_linear_benchmark_v2.log, replace text name(surveylinear)
di as result "BENCHMARK_FILE_REVISION=2"
di as txt "Stata version: `c(stata_version)'"
capture noisily which suest2
local available = (_rc == 0)
import delimited using suest_r_survey_linear_benchmark.csv, clear varnames(1) asdouble

* Full design is always retained. Domains use subpop(), never an if restriction.
forvalues case = 1/6 {
    capture drop dom_a dom_b
    generate byte dom_a = 1
    generate byte dom_b = 1
    if `case' == 3 {
        replace dom_a = unit <= 7
        replace dom_b = unit >= 4
    }
    if `case' == 4 {
        replace dom_a = unit <= 5
        replace dom_b = unit > 5
    }
    if `case' == 5 {
        replace dom_a = mod(psu, 6) == 1
        replace dom_b = mod(psu, 6) == 2
    }
    if `case' == 6 {
        replace dom_a = strata <= 2
        replace dom_b = strata > 2
    }
    if `case' == 1 {
        svyset psu [pweight = w], strata(strata) singleunit(missing)
    }
    else {
        svyset psu [pweight = w], strata(strata) fpc(pop) singleunit(missing)
    }
    di as result "SURVEY_CASE=`case'"
    quietly svy, subpop(dom_a): regress y x
    estimates store survey_a
    di as txt "NATIVE_MODEL=A DF=" e(df_r)
    matrix list e(b), format(%21.15g)
    matrix list e(V), format(%21.15g)
    * Cases 1-2 deliberately use the same subpopulation expression. This
    * distinguishes suest2's call-level compatibility check from a substantive
    * difference in estimation samples.
    if `case' <= 2 {
        quietly svy, subpop(dom_a): regress y2 x z
    }
    else {
        quietly svy, subpop(dom_b): regress y2 x z
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
            di as txt "JOINT_MODEL"
            matrix list e(b), format(%21.15g)
            matrix list e(V), format(%21.15g)
        }
        capture noisily suest2_cleanup, force
    }
    else {
        di as result "SUEST2_UNAVAILABLE=1"
    }

    * Independent full-joint reference: stack the two equations as one survey
    * regression. Separate intercept/slope columns make its normal equations
    * block diagonal, while PSU totals and stratum centering produce the same
    * cross-equation covariance as the R joint linearization.
    preserve
    expand 2, generate(equation_b)
    generate double y_stack = cond(equation_b == 0, y, y2)
    generate byte stack_subpop = cond(equation_b == 0, dom_a, dom_b)
    generate byte a_cons = equation_b == 0
    generate double a_x = x * a_cons
    generate byte b_cons = equation_b == 1
    generate double b_x = x * b_cons
    generate double b_z = z * b_cons
    quietly svy, subpop(stack_subpop): regress y_stack a_cons a_x b_cons b_x b_z, noconstant
    di as txt "STACKED_JOINT_MODEL DF=" e(df_r)
    matrix list e(b), format(%21.15g)
    matrix list e(V), format(%21.15g)
    restore
}
di as result "SURVEY_LINEAR_BENCHMARK_V2_COMPLETE=1"
log close surveylinear
