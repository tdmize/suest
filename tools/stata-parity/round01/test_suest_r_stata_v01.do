* R suest vs Stata suest2: (1) 31 matched ME helper estimate/SE pairs with
* R references regenerated under the combined-sample N/(N-1) correction;
* (2) joint e(b)/e(V) for model families with no prior Stata comparison and
* for different-sample ordinary systems. Run from the folder holding the inputs.
version 16
clear all
set more off
set linesize 255
capture log close
log using test_suest_r_stata_v01.log, text replace
about
foreach cmd in suest2 meinequality totalme {
    which `cmd'
}
foreach f in matched_fixture.dta r_reference_v04.csv parity_fixture_v01.dta {
    confirm file `f'
}

*-------------------------------------------------------------------------------
* PART 1. Matched ME helper gate (31 pairs)
*-------------------------------------------------------------------------------
use matched_fixture.dta, clear
tempfile me_post
capture postclose __results
postfile __results str40 case int row double estimate std_error using `me_post', replace
capture program drop record_result
program define record_result
    syntax , CASE(string) [ROWS(integer 0)]
    tempname table
    matrix `table' = r(table)
    local n = rowsof(`table')
    if `rows' > 0 {
        if `rows' > `n' exit 498
        local n = `rows'
    }
    forvalues j = 1/`n' {
        if missing(`table'[`j',1], `table'[`j',2]) {
            display as error "Missing estimate or SE in `case', row `j'."
            exit 498
        }
        post __results ("`case'") (`j') (`table'[`j',1]) (`table'[`j',2])
    }
end

capture noisily {
    quietly regress y i.f i.b c.z
    meinequality f
    record_result, case(linear_weighted_shares)
    meinequality f, unweighted
    record_result, case(linear_equal_pairs)

    quietly logit yb i.f i.b c.z, iterate(1000) tolerance(1e-10)
    meinequality f
    record_result, case(logit_inequality)
    quietly probit yb i.f i.b c.z, iterate(1000) tolerance(1e-10)
    meinequality f
    record_result, case(probit_inequality)

    quietly mlogit ym i.f i.b c.z, baseoutcome(1) iterate(1000) tolerance(1e-10)
    meinequality f
    record_result, case(multinomial_inequality)
    meinequality f, unweighted
    record_result, case(multinomial_equal_pairs)
    foreach v in f b z {
        totalme `v', centered
        record_result, case(multinomial_total_`v')
    }
    quietly ologit yo i.f i.b c.z, iterate(1000) tolerance(1e-10)
    meinequality f
    record_result, case(ordered_inequality)
    foreach v in f b z {
        totalme `v', centered
        record_result, case(ordered_total_`v')
    }

    quietly mlogit ym i.f i.b c.z [fw=w], baseoutcome(1) iterate(1000) tolerance(1e-10)
    totalme f
    record_result, case(weighted_multinomial_total)

    * Subgroup s has no value labels: totalme 1.8.2 should accept it (1.8.1 gave r(198)).
    quietly mlogit ym i.f##i.s i.b c.z, baseoutcome(1) iterate(1000) tolerance(1e-10)
    estimates store SubgroupFit
    capture noisily totalme f, over(s)
    display as text "Unlabeled subgroup return code: " _rc
    quietly estimates restore SubgroupFit
    label define ME_subgroup 1 "First" 2 "Second", replace
    label values s ME_subgroup
    totalme f, over(s)
    record_result, case(subgroup_total) rows(2)

    quietly mlogit ym i.f i.b c.z if id <= 300, baseoutcome(1) vce(robust) iterate(1000) tolerance(1e-10)
    estimates store First
    quietly mlogit ym i.f i.b c.z if id > 300, baseoutcome(1) vce(robust) iterate(1000) tolerance(1e-10)
    estimates store Second
    totalme f, models(First Second) groups
    record_result, case(disjoint_total)

    quietly regress y i.f i.b c.z, vce(robust)
    estimates store Base
    quietly regress y i.f i.b c.z i.s, vce(robust)
    estimates store Full
    meinequality f, models(Base Full)
    record_result, case(same_sample_inequality)
}
local part1_rc = _rc
postclose __results
local me_failed = .
if `part1_rc' == 0 {
    use `me_post', clear
    export delimited using stata_results_v04.csv, replace
    tempfile reference
    preserve
    import delimited using r_reference_v04.csv, varnames(1) stringcols(1) asdouble clear
    save `reference'
    restore
    merge 1:1 case row using `reference'
    generate double estimate_error = estimate - r_estimate
    generate double se_error = std_error - r_std_error
    generate byte pass = _merge == 3 & !missing(estimate_error, se_error) & ///
        abs(estimate_error) <= 2e-5 + 1e-4 * abs(r_estimate) & ///
        abs(se_error) <= 2e-5 + 1e-4 * abs(r_std_error)
    format estimate std_error r_estimate r_std_error estimate_error se_error %21.17g
    sort case row
    export delimited using matched_comparison_v04.csv, replace
    list case row estimate r_estimate std_error r_std_error pass, noobs sep(0)
    quietly count if pass == 0
    local me_failed = r(N)
    display as text "PART 1: " _N - `me_failed' " of " _N " pairs passed (31 expected)."
}
else display as error "PART 1 stopped with r(`part1_rc'); see above."

*-------------------------------------------------------------------------------
* PART 2. Joint e(b) and e(V) for families without a Stata comparison
*-------------------------------------------------------------------------------
use parity_fixture_v01.dta, clear
tempfile par_post
capture postclose __par
postfile __par str24 case str4 kind str60 row str60 col double value using `par_post', replace

capture program drop run_suest2
program define run_suest2
    args case fitrc
    if `fitrc' {
        display as error "`case': model fit failed, r(`fitrc')"
        post __par ("`case'") ("rc") ("fit") ("") (`fitrc')
        exit
    }
    capture noisily suest2 A B
    if _rc {
        display as error "`case': suest2 failed, r(" _rc ")"
        post __par ("`case'") ("rc") ("suest2") ("") (_rc)
        exit
    }
    tempname b V
    matrix `b' = e(b)
    matrix `V' = e(V)
    local names : colfullnames `b'
    local k = colsof(`b')
    forvalues i = 1/`k' {
        local ni : word `i' of `names'
        post __par ("`case'") ("b") ("`ni'") ("") (`b'[1,`i'])
        forvalues j = 1/`k' {
            local nj : word `j' of `names'
            post __par ("`case'") ("V") ("`ni'") ("`nj'") (`V'[`i',`j'])
        }
    }
    post __par ("`case'") ("N") ("union") ("") (e(N))
    display as text "`case': stored " `k' " parameters, N = " e(N)
    capture suest2_cleanup, force
end

capture program drop fit_pair
program define fit_pair
    args case cmdA cmdB
    capture estimates drop A B
    capture noisily {
        quietly `cmdA'
        estimates store A
        quietly `cmdB'
        estimates store B
    }
    run_suest2 `case' `=_rc'
end

fit_pair cloglog        "cloglog yb x1 x2, vce(robust)"                           "cloglog yb x1 x2 x3, vce(robust)"
fit_pair gamma_log      "glm ygam x1 x2, family(gamma) link(log) vce(robust)"     "glm ygam x1 x2 x3, family(gamma) link(log) vce(robust)"
fit_pair fracreg_logit  "fracreg logit yf x1 x2, vce(robust)"                     "fracreg logit yf x1 x2 x3, vce(robust)"
fit_pair fracreg_probit "fracreg probit yf x1 x2, vce(robust)"                    "fracreg probit yf x1 x2 x3, vce(robust)"
fit_pair zip            "zip yz x1 x2, inflate(x2) vce(robust)"                   "zip yz x1 x2 x3, inflate(x2) vce(robust)"
fit_pair zinb           "zinb yz x1 x2, inflate(x2) vce(robust)"                  "zinb yz x1 x2 x3, inflate(x2) vce(robust)"
fit_pair tobit          "tobit ytob x1 x2, ll(0) vce(robust)"                     "tobit ytob x1 x2 x3, ll(0) vce(robust)"
fit_pair intreg         "intreg lo hi x1 x2, vce(robust)"                         "intreg lo hi x1 x2 x3, vce(robust)"
fit_pair truncreg       "truncreg ytr x1 x2, ll(0) vce(robust)"                   "truncreg ytr x1 x2 x3, ll(0) vce(robust)"
fit_pair betareg        "betareg ybeta x1 x2, vce(robust)"                        "betareg ybeta x1 x2 x3, vce(robust)"
fit_pair hetprobit      "hetprobit yh x1 x2, het(x3) vce(robust)"                 "hetprobit yh x1 x2 x3, het(x3) vce(robust)"

* Different samples, no clusters: partial overlap (id <= 500 vs id > 300) and disjoint
fit_pair regress_partial  "regress ylin x1 x2 if id <= 500, vce(robust)"          "regress ylin x1 x2 x3 if id > 300, vce(robust)"
fit_pair regress_disjoint "regress ylin x1 x2 if id <= 400, vce(robust)"          "regress ylin x1 x2 x3 if id > 400, vce(robust)"
fit_pair logit_partial    "logit ylogit x1 x2 if id <= 500, vce(robust)"          "logit ylogit x1 x2 x3 if id > 300, vce(robust)"
fit_pair logit_disjoint   "logit ylogit x1 x2 if id <= 400, vce(robust)"          "logit ylogit x1 x2 x3 if id > 400, vce(robust)"
fit_pair cloglog_partial  "cloglog yb x1 x2 if id <= 500, vce(robust)"            "cloglog yb x1 x2 x3 if id > 300, vce(robust)"
fit_pair tobit_partial    "tobit ytob x1 x2 if id <= 500, ll(0) vce(robust)"      "tobit ytob x1 x2 x3 if id > 300, ll(0) vce(robust)"

* Parametric survival, accelerated failure-time metric
stset time, failure(event)
fit_pair streg_weibull     "streg x1 x2, distribution(weibull) time vce(robust)"      "streg x1 x2 x3, distribution(weibull) time vce(robust)"
fit_pair streg_exponential "streg x1 x2, distribution(exponential) time vce(robust)"  "streg x1 x2 x3, distribution(exponential) time vce(robust)"
fit_pair streg_lognormal   "streg x1 x2, distribution(lognormal) vce(robust)"         "streg x1 x2 x3, distribution(lognormal) vce(robust)"
fit_pair streg_loglogistic "streg x1 x2, distribution(loglogistic) vce(robust)"       "streg x1 x2 x3, distribution(loglogistic) vce(robust)"

postclose __par
use `par_post', clear
export delimited using parity_stata_v01.csv, replace
quietly count if kind == "rc"
local par_failed = r(N)
quietly levelsof case if kind == "N", local(done)
display as text "PART 2: " `: word count `done'' " of 21 systems stored; " `par_failed' " failed."
if `par_failed' list case row value if kind == "rc", noobs sep(0)

display as text _n "Return: test_suest_r_stata_v01.log, matched_comparison_v04.csv, parity_stata_v01.csv"
log close
