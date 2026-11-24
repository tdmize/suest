* test_suest_r_stata_v21.do
* R suest vs Stata, round 21: predictions and marginal effects for each model
* type (round 20 stopped at a Mata error). Changes from v20: the Mata
* function's argument case (a reserved word in Mata) is now cs; xtgee is now
* xtlogit, pa, which suest2 accepts; ivprobit drops the ml option (maximum
* likelihood is the default).
* For each of 38 pairs of models (A: x z f2 f3; B: x f2 f3), suest2 A B,
* then:
*   - margins           the default prediction of each model (MG|case|pred)
*   - margins, dydx(x)  the average marginal effect of x (MG|case|dydx)
*   - mecompare x, models(A B) uncentered: the effect of a one-unit increase
*     in x in each model and the difference (MC|case|change)
*   - for some models, other predict() options (MG|case|pred_<option>)
* plus the joint e(b) and e(V) (EV|case). Each fit prints a FIT| line.
* Every iterative model has iterate(50), and set maxiter 50 caps anything
* else.
* File (same folder): v20_fixture.dta. R values: v20_r_reference.csv
* (make_v20.R). Comparison: compare_v20.R.
version 16
clear all
set more off
set linesize 255
set maxiter 50
capture log close
log using test_suest_r_stata_v21.log, text replace

* Print a row vector and its standard errors, one element per line
capture mata: mata drop dumpvec()
mata:
void dumpvec(string scalar tag, string scalar cs, string scalar stat,
    string scalar bname, string scalar vname)
{
    real matrix b, V
    string matrix s
    real scalar i
    b = st_matrix(bname)
    V = st_matrix(vname)
    s = st_matrixcolstripe(bname)
    for (i = 1; i <= cols(b); i++) {
        printf("%s|%s|%s|%s:%s|%21.0g|%21.0g\n", tag, cs, stat,
            s[i, 1], s[i, 2], b[1, i], sqrt(V[i, i]))
    }
}
end

* Print e(b) and e(V) one element per line, at full precision
capture program drop dumpev
program define dumpev
    args case
    capture confirm matrix e(b)
    if _rc {
        display as error "EV|`case'|no results"
        exit
    }
    tempname b V
    matrix `b' = e(b)
    matrix `V' = e(V)
    local names : colfullnames `b'
    local k = colsof(`b')
    forvalues i = 1/`k' {
        local ri : word `i' of `names'
        display as text "EV|`case'|b|`ri'||" %21.0g `b'[1, `i']
    }
    forvalues i = 1/`k' {
        local ri : word `i' of `names'
        forvalues j = 1/`k' {
            local cj : word `j' of `names'
            display as text "EV|`case'|V|`ri'|`cj'|" %21.0g `V'[`i', `j']
        }
    }
end

* fitstore NAME CASE: command   fits quietly, stores the results as NAME, and
* reports the sample and convergence
capture program drop fitstore
program define fitstore
    gettoken name 0 : 0
    gettoken case 0 : 0, parse(":")
    gettoken colon 0 : 0, parse(":")
    capture quietly `0'
    if _rc {
        display as error "FIT|`case'|`name'|failed with r(" _rc ")"
        ereturn clear
        exit
    }
    estimates store `name'
    display as text "FIT|`case'|`name'|converged=" e(converged) "|N=" e(N)
end

* margdump CASE STAT: margins ...   runs margins and prints every estimate
capture program drop margdump
program define margdump
    gettoken case 0 : 0
    gettoken stat 0 : 0, parse(":")
    gettoken colon 0 : 0, parse(":")
    capture noisily `0'
    if _rc {
        display as error "MG|`case'|`stat'|failed with r(" _rc ")"
        exit
    }
    mata: dumpvec("MG", "`case'", "`stat'", "r(b)", "r(V)")
end

* sysstart CASE   suest2 A B, the joint results, the default prediction of
* each model, and the average marginal effect of x
capture program drop sysstart
program define sysstart
    args case
    capture noisily suest2 A B
    if _rc {
        display as error "EV|`case'|suest2 failed with r(" _rc ")"
        exit
    }
    dumpev `case'
    margdump `case' pred: margins
    margdump `case' dydx: margins, dydx(x)
end

* sysend CASE   mecompare for a one-unit increase in x, then clean up
capture program drop sysend
program define sysend
    args case
    capture noisily mecompare x, models(A B) uncentered
    if _rc display as error "MC|`case'|change|failed with r(" _rc ")"
    else mata: dumpvec("MC", "`case'", "change", "e(b)", "e(V)")
    capture noisily suest2_cleanup, force
end

foreach cmd in suest2 suest2_cleanup mecompare gologit2 margins xtlogit {
    capture noisily which `cmd'
}
display "maxiter = " c(maxiter)

use v20_fixture.dta, clear
xtset id t
local A x z f2 f3
local B x f2 f3
local it iterate(50)

* ---------------------------------------------------------------------------
* Part 1. Single-level models
* ---------------------------------------------------------------------------
fitstore A regress: regress y `A'
fitstore B regress: regress y `B'
sysstart regress
sysend regress

fitstore A logit: logit yb `A', `it'
fitstore B logit: logit yb `B', `it'
sysstart logit
sysend logit

fitstore A ologit: ologit yo `A', `it'
fitstore B ologit: ologit yo `B', `it'
sysstart ologit
sysend ologit

fitstore A oprobit: oprobit yo `A', `it'
fitstore B oprobit: oprobit yo `B', `it'
sysstart oprobit
sysend oprobit

fitstore A mlogit: mlogit ym `A', baseoutcome(1) `it'
fitstore B mlogit: mlogit ym `B', baseoutcome(1) `it'
sysstart mlogit
sysend mlogit

fitstore A poisson: poisson yc `A', `it'
fitstore B poisson: poisson yc `B', `it'
sysstart poisson
sysend poisson

fitstore A nbreg: nbreg ynb `A', `it'
fitstore B nbreg: nbreg ynb `B', `it'
sysstart nbreg
sysend nbreg

fitstore A nbreg_constant: nbreg ynb `A', dispersion(constant) `it'
fitstore B nbreg_constant: nbreg ynb `B', dispersion(constant) `it'
sysstart nbreg_constant
sysend nbreg_constant

fitstore A glm_gamma: glm ypos `A', family(gamma) link(log) `it'
fitstore B glm_gamma: glm ypos `B', family(gamma) link(log) `it'
sysstart glm_gamma
sysend glm_gamma

fitstore A fracreg_logit: fracreg logit yf `A', `it'
fitstore B fracreg_logit: fracreg logit yf `B', `it'
sysstart fracreg_logit
sysend fracreg_logit

fitstore A fracreg_probit: fracreg probit yf `A', `it'
fitstore B fracreg_probit: fracreg probit yf `B', `it'
sysstart fracreg_probit
sysend fracreg_probit

fitstore A tobit: tobit ycens `A', ll(0) `it'
fitstore B tobit: tobit ycens `B', ll(0) `it'
sysstart tobit
sysend tobit

fitstore A intreg: intreg lo hi `A', `it'
fitstore B intreg: intreg lo hi `B', `it'
sysstart intreg
sysend intreg

fitstore A truncreg: truncreg ytr `A', ll(0) `it'
fitstore B truncreg: truncreg ytr `B', ll(0) `it'
sysstart truncreg
sysend truncreg

fitstore A betareg: betareg ybeta `A', `it'
fitstore B betareg: betareg ybeta `B', `it'
sysstart betareg
sysend betareg

fitstore A hetprobit: hetprobit yh `A', het(w) `it'
fitstore B hetprobit: hetprobit yh `B', het(w) `it'
sysstart hetprobit
sysend hetprobit

fitstore A biprobit: biprobit (yb1 = `A') (yb2 = `A'), `it'
fitstore B biprobit: biprobit (yb1 = `B') (yb2 = `B'), `it'
sysstart biprobit
sysend biprobit

fitstore A ivregress: ivregress 2sls yiv `A' (endog = instr1 instr2)
fitstore B ivregress: ivregress 2sls yiv `B' (endog = instr1 instr2)
sysstart ivregress
sysend ivregress

* R predicts the probability from the outcome equation alone, Stata's
* predict(pr fixedasf); Stata's default pr conditions on the first stage's errors
fitstore A ivprobit: ivprobit yivb `A' (endog = instr1 instr2), `it'
fitstore B ivprobit: ivprobit yivb `B' (endog = instr1 instr2), `it'
sysstart ivprobit
margdump ivprobit pred_pr_A: margins, predict(model(A) pr)
margdump ivprobit pred_pr_B: margins, predict(model(B) pr)
margdump ivprobit pred_fixedasf_A: margins, predict(model(A) pr fixedasf)
margdump ivprobit pred_fixedasf_B: margins, predict(model(B) pr fixedasf)
margdump ivprobit dydx_fixedasf_A: margins, dydx(x) predict(model(A) pr fixedasf)
margdump ivprobit dydx_fixedasf_B: margins, dydx(x) predict(model(B) pr fixedasf)
sysend ivprobit

fitstore A heckman: heckman yheck `A', select(sel = x w) `it'
fitstore B heckman: heckman yheck `B', select(sel = x w) `it'
sysstart heckman
sysend heckman

fitstore A gologit2: gologit2 yo `A', npl(x) `it'
fitstore B gologit2: gologit2 yo `B', pl `it'
sysstart gologit2
sysend gologit2

* ---------------------------------------------------------------------------
* Part 2. Panel models (R: plm, nlme, pglm, fixest, geepack)
* ---------------------------------------------------------------------------
fitstore A xtreg_fe: xtreg y `A', fe
fitstore B xtreg_fe: xtreg y `B', fe
sysstart xtreg_fe
sysend xtreg_fe

fitstore A xtreg_re: xtreg y `A', re
fitstore B xtreg_re: xtreg y `B', re
sysstart xtreg_re
sysend xtreg_re

fitstore A xtreg_mle: xtreg y `A', mle `it'
fitstore B xtreg_mle: xtreg y `B', mle `it'
sysstart xtreg_mle
sysend xtreg_mle

* xtlogit and xtprobit, re predict xb by default; pr is the probability
* averaged over the random effect (R's prediction)
* suest2 needs adaptive quadrature here; R's pglm uses nonadaptive quadrature
* with 12 points, so the two differ slightly
local q12 intmethod(mvaghermite) intpoints(12) `it'
fitstore A xtlogit_re: xtlogit yb `A', re `q12'
fitstore B xtlogit_re: xtlogit yb `B', re `q12'
sysstart xtlogit_re
margdump xtlogit_re pred_pr_A: margins, predict(model(A) pr)
margdump xtlogit_re pred_pr_B: margins, predict(model(B) pr)
margdump xtlogit_re dydx_pr_A: margins, dydx(x) predict(model(A) pr)
margdump xtlogit_re dydx_pr_B: margins, dydx(x) predict(model(B) pr)
sysend xtlogit_re

fitstore A xtprobit_re: xtprobit yb `A', re `q12'
fitstore B xtprobit_re: xtprobit yb `B', re `q12'
sysstart xtprobit_re
margdump xtprobit_re pred_pr_A: margins, predict(model(A) pr)
margdump xtprobit_re pred_pr_B: margins, predict(model(B) pr)
margdump xtprobit_re dydx_pr_A: margins, dydx(x) predict(model(A) pr)
margdump xtprobit_re dydx_pr_B: margins, dydx(x) predict(model(B) pr)
sysend xtprobit_re

* xtpoisson predicts xb by default; nu0 is the expected count with the random
* (or fixed) effect at zero
fitstore A xtpoisson_re: xtpoisson yc `A', re `it'
fitstore B xtpoisson_re: xtpoisson yc `B', re `it'
sysstart xtpoisson_re
margdump xtpoisson_re pred_nu0_A: margins, predict(model(A) nu0)
margdump xtpoisson_re pred_nu0_B: margins, predict(model(B) nu0)
margdump xtpoisson_re dydx_nu0_A: margins, dydx(x) predict(model(A) nu0)
margdump xtpoisson_re dydx_nu0_B: margins, dydx(x) predict(model(B) nu0)
sysend xtpoisson_re

fitstore A xtpoisson_fe: xtpoisson yc `A', fe `it'
fitstore B xtpoisson_fe: xtpoisson yc `B', fe `it'
sysstart xtpoisson_fe
margdump xtpoisson_fe pred_nu0_A: margins, predict(model(A) nu0)
margdump xtpoisson_fe pred_nu0_B: margins, predict(model(B) nu0)
sysend xtpoisson_fe

* Population-averaged logit with R's exchangeable working correlation
* (geepack's estimate) fixed, since Stata estimates the correlation a little
* differently; suest2 does not accept xtgee itself
matrix RA = J(6, 6, 0.09535210004788397)
matrix RB = J(6, 6, 0.093034942275636051)
forvalues i = 1/6 {
    matrix RA[`i', `i'] = 1
    matrix RB[`i', `i'] = 1
}
fitstore A xtgee_logit: xtlogit yb `A', pa corr(fixed RA) `it'
fitstore B xtgee_logit: xtlogit yb `B', pa corr(fixed RB) `it'
sysstart xtgee_logit
sysend xtgee_logit

* ---------------------------------------------------------------------------
* Part 3. Multilevel models (R: lme4 with 25-point adaptive quadrature,
*         glmmTMB, ordinal)
* ---------------------------------------------------------------------------
fitstore A mixed: mixed y `A' || id: x, covariance(unstructured) mle `it'
fitstore B mixed: mixed y `B' || id: x, covariance(unstructured) mle `it'
sysstart mixed
sysend mixed

local q25 intmethod(mvaghermite) intpoints(25) `it'
fitstore A melogit: melogit yb `A' || id:, `q25'
fitstore B melogit: melogit yb `B' || id:, `q25'
sysstart melogit
sysend melogit

fitstore A meprobit: meprobit yb `A' || id:, `q25'
fitstore B meprobit: meprobit yb `B' || id:, `q25'
sysstart meprobit
sysend meprobit

fitstore A mepoisson: mepoisson yc `A' || id:, `q25'
fitstore B mepoisson: mepoisson yc `B' || id:, `q25'
sysstart mepoisson
sysend mepoisson

* R's glmmTMB uses the Laplace approximation, which suest2 does not accept, so
* these results are close but not identical
fitstore A menbreg: menbreg yc `A' || id:, `q25'
fitstore B menbreg: menbreg yc `B' || id:, `q25'
sysstart menbreg
sysend menbreg

fitstore A meologit: meologit yo `A' || id:, `q25'
fitstore B meologit: meologit yo `B' || id:, `q25'
sysstart meologit
sysend meologit

* ---------------------------------------------------------------------------
* Part 4. Survey and probability-weighted models
* ---------------------------------------------------------------------------
svyset psu [pweight = wt], strata(strata)
fitstore A svy_logit: svy: logit yb `A', `it'
fitstore B svy_logit: svy: logit yb `B', `it'
sysstart svy_logit
sysend svy_logit

fitstore A svy_ologit: svy: ologit yo `A', `it'
fitstore B svy_ologit: svy: ologit yo `B', `it'
sysstart svy_ologit
sysend svy_ologit

svyset, clear
fitstore A pweight_logit: logit yb `A' [pweight = wt], `it'
fitstore B pweight_logit: logit yb `B' [pweight = wt], `it'
sysstart pweight_logit
sysend pweight_logit

log close
