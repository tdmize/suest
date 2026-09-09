version 16
clear all
set more off
capture log close _all
log using suest_r_gee_crosslang_benchmark.log, replace text name(geecross)
di as result "BENCHMARK_FILE_REVISION=1"
di as txt "Stata version: `c(stata_version)'"
which suest2
import delimited using suest_r_gee_crosslang_benchmark.csv, clear varnames(1)
xtset id time

* Store conventional estimates: suest2 supplies the joint robust VCE.
foreach corr in independent exchangeable {
    foreach cmd in xtreg xtlogit xtprobit xtcloglog xtpoisson {
        local outcome binary
        if "`cmd'" == "xtreg" local outcome y
        if "`cmd'" == "xtpoisson" local outcome count
        di as txt "=== GEE CASE: `cmd' `corr' ==="
        quietly `cmd' `outcome' x, pa corr(`corr')
        estimates store gee_base
        quietly `cmd' `outcome' x z, pa corr(`corr')
        estimates store gee_adjusted
        suest2 gee_base gee_adjusted
        di as txt "N=" e(N) " N_clust=" e(N_clust)
        matrix list e(b), format(%21.15g)
        matrix list e(V), format(%21.15g)
        capture noisily suest2_cleanup, force
    }
}
di as txt "=== GEE CASE: MIXED FAMILY, UNEQUAL SAMPLES, HIGHER CLUSTERS ==="
quietly xtreg y x z if id <= 100, pa corr(exchangeable)
estimates store gee_gaussian
quietly xtlogit binary x z if id >= 21 & time != 4, pa corr(exchangeable)
estimates store gee_logit
suest2 gee_gaussian gee_logit, cluster(higher)
di as txt "N=" e(N) " N_clust=" e(N_clust)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force
di as result "GEE_CROSSLANG_BENCHMARK_COMPLETE=1"
log close geecross
