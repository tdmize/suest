* test_twopart_stata_v03.do
* Two-part models after suest2 and mecompare: reference values for the R
* suest twopart option (MEPS 2004, as on the suest2 and mecompare pages).
* Lines start ME| (margins: case, term, estimate, SE) and FIT| (fit details).
version 16
clear all
set more off
set linesize 255
capture log close
log using test_twopart_stata_v03.log, text replace

foreach cmd in suest2 mecompare {
    which `cmd'
}

capture program drop dumpme
program define dumpme
    args case
    tempname t
    matrix `t' = r(table)
    local names : colfullnames `t'
    local k = colsof(`t')
    forvalues i = 1/`k' {
        local ri : word `i' of `names'
        display as text "ME|`case'|`ri'|" %21.0g `t'[1, `i'] "|" %21.0g `t'[2, `i']
    }
end

use "https://tdmize.github.io/data/data/meps", clear
generate byte any = exp_tot > 0
generate double lnexp = ln(exp_tot) if any
local x "c.age i.female c.lninc c.pcs12 c.mcs12 i.ins_unins i.anylim"

* 1. suest2 page: logit + Poisson, clustered on duid
quietly logit any `x', vce(cluster duid) iterate(50)
display as text "FIT|Any|converged=" e(converged) "|N=" e(N)
estimates store Any
quietly poisson exp_tot `x' if any, vce(cluster duid) iterate(50)
display as text "FIT|Amount|converged=" e(converged) "|N=" e(N)
estimates store Amount
suest2 Any Amount
margins, predict(twopart(Any Amount))
dumpme tp_mean
margins, dydx(age female ins_unins anylim) predict(twopart(Any Amount))
dumpme tp_dydx
margins, dydx(ins_unins) predict(model(Any))
dumpme tp_any
margins, dydx(ins_unins) predict(model(Amount))
dumpme tp_amount
margins if any, dydx(ins_unins) predict(model(Amount))
dumpme tp_amount_ifany
margins, dydx(ins_unins) over(female) predict(twopart(Any Amount))
dumpme tp_over
margins, dydx(ins_unins) over(r.female) predict(twopart(Any Amount))
dumpme tp_overr
predict double yhat if e(sample), twopart(Any Amount)
correlate exp_tot yhat
display as text "FIT|tp_cor|rho=" %21.0g r(rho)

* 2. Other second parts (and a probit first part), clustered on duid
quietly glm exp_tot `x' if any, family(gamma) link(log) vce(cluster duid) iterate(50)
display as text "FIT|Gamma|converged=" e(converged) "|N=" e(N)
estimates store Gamma
quietly nbreg exp_tot `x' if any, vce(cluster duid) iterate(50)
display as text "FIT|NegBin|converged=" e(converged) "|N=" e(N)
estimates store NegBin
quietly regress exp_tot `x' if any, vce(cluster duid)
estimates store Linear
quietly regress lnexp `x' if any, vce(cluster duid)
estimates store LnAmount
quietly probit any `x', vce(cluster duid) iterate(50)
display as text "FIT|Probit|converged=" e(converged) "|N=" e(N)
estimates store Probit

foreach m in Gamma NegBin Linear {
    suest2 Any `m'
    margins, predict(twopart(Any `m'))
    dumpme `m'_mean
    margins, dydx(age ins_unins) predict(twopart(Any `m'))
    dumpme `m'_dydx
}
suest2 Probit Amount
margins, predict(twopart(Probit Amount))
dumpme Probit_mean
margins, dydx(age ins_unins) predict(twopart(Probit Amount))
dumpme Probit_dydx

* 3. Log outcome: Duan's smearing estimate and the normal-theory retransformation
suest2 Any LnAmount
foreach r in duan normal {
    margins, predict(twopart(Any LnAmount) `r')
    dumpme `r'_mean
    margins, dydx(age ins_unins) predict(twopart(Any LnAmount) `r')
    dumpme `r'_dydx
}

* 4. Sampling weights on both parts, clustered on duid
quietly logit any `x' [pw=wtdper], iterate(50)
estimates store AnyW
quietly poisson exp_tot `x' if any [pw=wtdper], iterate(50)
estimates store AmountW
quietly glm exp_tot `x' if any [pw=wtdper], family(gamma) link(log) iterate(50)
estimates store GammaW
foreach m in AmountW GammaW {
    suest2 AnyW `m', cluster(duid)
    display as text "FIT|`m'|wtype=" e(wtype) "|wexp=" e(wexp) "|N_clust=" e(N_clust)
    margins, predict(twopart(AnyW `m'))
    dumpme `m'_mean
    margins, dydx(age ins_unins) predict(twopart(AnyW `m'))
    dumpme `m'_dydx
    margins, predict(twopart(AnyW `m')) noweights
    dumpme `m'_mean_now
    margins, dydx(age ins_unins) predict(twopart(AnyW `m')) noweights
    dumpme `m'_dydx_now
}

* 5. Survey design: weights and clusters from svyset; second part with if any
svyset duid [pw=wtdper]
quietly svy: logit any `x', iterate(50)
estimates store AnyS
quietly svy: poisson exp_tot `x' if any, iterate(50)
estimates store AmountS
suest2 AnyS AmountS
margins, predict(twopart(AnyS AmountS))
dumpme svy_mean
margins, dydx(age ins_unins) predict(twopart(AnyS AmountS))
dumpme svy_dydx

* 6. mecompare page: robust standard errors, four predictors
local x "c.age i.female c.lninc i.ins_unins"
quietly logit any `x', vce(robust) iterate(50)
estimates store Any
quietly glm exp_tot `x' if any, family(gamma) link(log) vce(robust) iterate(50)
estimates store Gamma
quietly poisson exp_tot `x' if any, vce(robust) iterate(50)
estimates store Poisson
quietly poisson exp_tot `x', vce(robust) iterate(50)
estimates store OnePart
mecompare ins_unins age, models(Any Gamma Poisson) twopart(Any Gamma Any Poisson)
mecompare ins_unins age, models(Any Poisson OnePart) twopart(Any Poisson)

log close
