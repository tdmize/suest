* test_suest_r_stata_v19.do
* R suest vs Stata, round 19: median survival time predictions after streg
* (Stata's default prediction) for Weibull, exponential, lognormal, and
* log-logistic models. For each model's own vce(robust) fit: the average
* median survival time and the average effect of x on it (margins), and the
* log likelihood. Every model has iterate(50), and set maxiter 50 caps
* anything else.
* File (same folder): batch4_fixture_v01.dta (from round 18). R values:
* survreg_median_r_reference_v01.csv (make_v19.R).
version 16
clear all
set more off
set linesize 255
set maxiter 50
capture log close
log using test_suest_r_stata_v19.log, text replace

* margdump CASE: margins ...   runs margins and prints its estimate and
* standard error
capture program drop margdump
program define margdump
    gettoken case 0 : 0, parse(":")
    gettoken colon 0 : 0, parse(":")
    capture noisily `0'
    if _rc {
        display as error "M|`case'|margins failed with r(" _rc ")"
        exit
    }
    tempname b V
    matrix `b' = r(b)
    matrix `V' = r(V)
    display as text "M|`case'|" %21.0g `b'[1, 1] "|" %21.0g sqrt(`V'[1, 1])
end

* survfit DIST   fits streg with vce(robust), then prints the log likelihood,
* the average median survival time, and the average effect of x on it
capture program drop survfit
program define survfit
    args dist
    capture noisily streg x z, distribution(`dist') vce(robust) iterate(50) nrtolerance(1e-10)
    if _rc {
        display as error "M|`dist'|streg failed with r(" _rc ")"
        exit
    }
    display as text "FIT|`dist'|converged=" e(converged) "|N=" e(N) "|N_fail=" e(N_fail)
    display as text "M|`dist'_ll|" %21.0g e(ll) "|"
    margdump `dist'_median: margins, predict(median time)
    margdump `dist'_dydx_x: margins, dydx(x) predict(median time)
end

foreach cmd in stset streg margins {
    capture noisily which `cmd'
}
display "maxiter = " c(maxiter)

use batch4_fixture_v01.dta, clear
stset tm, failure(ev)
foreach dist in weibull exponential lognormal loglogistic {
    survfit `dist'
}

log close
