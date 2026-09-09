version 16
clear all
set more off

capture log close _all
log using suest_r_panel_fe_crosslang_benchmark.log, replace text name(panelcross)

di as result "BENCHMARK_FILE_REVISION=1"
di as txt "SUEST R PANEL FE CROSS-LANGUAGE BENCHMARK"
di as txt "Stata version: `c(stata_version)'"
which suest2

import delimited using suest_r_panel_fe_crosslang_benchmark.csv, clear varnames(1)
xtset id time

di as txt "=== CASE 1: BALANCED PANEL, DEFAULT PANEL CLUSTERS ==="
quietly xtreg y1 x z i.time, fe
estimates store fe1_balanced
quietly xtreg y2 x z i.time, fe
estimates store fe2_balanced
suest2 fe1_balanced fe2_balanced
di as txt "N=" e(N) " N_clust=" e(N_clust)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as txt "=== CASE 2: PARTIALLY OVERLAPPING UNBALANCED PANELS ==="
quietly xtreg y1_partial x z i.time, fe
estimates store fe1_partial
quietly xtreg y2_partial x z i.time, fe
estimates store fe2_partial
suest2 fe1_partial fe2_partial
di as txt "N=" e(N) " N_clust=" e(N_clust)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as txt "=== CASE 3: HIGHER-LEVEL CLUSTERS ==="
quietly xtreg y1_partial x z i.time, fe
estimates store fe1_higher
quietly xtreg y2_partial x z i.time, fe
estimates store fe2_higher
suest2 fe1_higher fe2_higher, cluster(higher)
di as txt "N=" e(N) " N_clust=" e(N_clust)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as result "PANEL_FE_CROSSLANG_BENCHMARK_COMPLETE=1"
log close panelcross
