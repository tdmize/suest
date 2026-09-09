version 16
clear all
set more off

capture log close _all
log using suest_r_panel_be_re_crosslang_benchmark.log, replace text name(panelbere)

di as result "BENCHMARK_FILE_REVISION=1"
di as txt "SUEST R PANEL BE/RE CROSS-LANGUAGE BENCHMARK"
di as txt "Stata version: `c(stata_version)'"
which suest2

import delimited using suest_r_panel_fe_crosslang_benchmark.csv, clear varnames(1)
xtset id time

di as txt "=== CASE 1: BETWEEN, BALANCED, DEFAULT PANEL CLUSTERS ==="
quietly xtreg y1 x z, be
estimates store be1_balanced
quietly xtreg y2 x z, be
estimates store be2_balanced
suest2 be1_balanced be2_balanced
di as txt "N=" e(N) " N_clust=" e(N_clust)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as txt "=== CASE 2: BETWEEN, UNBALANCED, HIGHER CLUSTERS ==="
quietly xtreg y1_partial x z, be
estimates store be1_partial
quietly xtreg y2_partial x z, be
estimates store be2_partial
suest2 be1_partial be2_partial, cluster(higher)
di as txt "N=" e(N) " N_clust=" e(N_clust)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as txt "=== CASE 3: RANDOM EFFECTS, BALANCED, DEFAULT PANEL CLUSTERS ==="
quietly xtreg y1 x z i.time, re
estimates store re1_balanced
quietly xtreg y2 x z i.time, re
estimates store re2_balanced
suest2 re1_balanced re2_balanced
di as txt "N=" e(N) " N_clust=" e(N_clust)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as txt "=== CASE 4: RANDOM EFFECTS, UNBALANCED, HIGHER CLUSTERS ==="
quietly xtreg y1_partial x z i.time, re
estimates store re1_partial
quietly xtreg y2_partial x z i.time, re
estimates store re2_partial
suest2 re1_partial re2_partial, cluster(higher)
di as txt "N=" e(N) " N_clust=" e(N_clust)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as result "PANEL_BE_RE_CROSSLANG_BENCHMARK_COMPLETE=1"
log close panelbere
