version 16
clear all
set more off

capture log close _all
log using suest_r_panel_ml_crosslang_benchmark.log, replace text name(panelml)

di as result "BENCHMARK_FILE_REVISION=1"
di as txt "SUEST R PANEL ML CROSS-LANGUAGE BENCHMARK"
di as txt "Stata version: `c(stata_version)'"
which suest2

import delimited using suest_r_panel_fe_crosslang_benchmark.csv, clear varnames(1)
xtset id time

di as txt "=== CASE 1: RANDOM-INTERCEPT ML, BALANCED, DEFAULT PANEL CLUSTERS ==="
quietly xtreg y1 x z i.time, mle
estimates store ml1_balanced
quietly xtreg y2 x z i.time, mle
estimates store ml2_balanced
suest2 ml1_balanced ml2_balanced
di as txt "N=" e(N) " N_clust=" e(N_clust) " N_g=" e(N_g)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as txt "=== CASE 2: RANDOM-INTERCEPT ML, UNBALANCED, HIGHER CLUSTERS ==="
quietly xtreg y1_partial x z i.time, mle
estimates store ml1_partial
quietly xtreg y2_partial x z i.time, mle
estimates store ml2_partial
suest2 ml1_partial ml2_partial, cluster(higher)
di as txt "N=" e(N) " N_clust=" e(N_clust) " N_g=" e(N_g)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as result "PANEL_ML_CROSSLANG_BENCHMARK_COMPLETE=1"
log close panelml
