version 16
clear all
set more off

capture log close _all
log using suest_r_panel_poisson_re_crosslang_benchmark.log, replace text name(panelrepoisson)

di as result "BENCHMARK_FILE_REVISION=1"
di as txt "SUEST R GAMMA RANDOM-EFFECTS PANEL POISSON CROSS-LANGUAGE BENCHMARK"
di as txt "Stata version: `c(stata_version)'"
which suest2

import delimited using suest_r_panel_poisson_re_crosslang_benchmark.csv, clear varnames(1)
xtset id time

di as txt "=== CASE 1A: BALANCED MODEL Y1 ==="
quietly xtpoisson y1 x z, re
di as txt "N=" e(N) " N_g=" e(N_g) " LL=" %21.15g e(ll)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
predict double y1_nu0 if e(sample), nu0
quietly summarize y1_nu0 if e(sample), meanonly
di as result "Y1_NU0_MEAN=" %21.15g r(mean)
quietly margins, dydx(x) predict(nu0)
matrix list r(b), format(%21.15g)
matrix list r(V), format(%21.15g)
estimates store y1_balanced

di as txt "=== CASE 1B: BALANCED MODEL Y2 ==="
quietly xtpoisson y2 x z, re
di as txt "N=" e(N) " N_g=" e(N_g) " LL=" %21.15g e(ll)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
predict double y2_nu0 if e(sample), nu0
quietly summarize y2_nu0 if e(sample), meanonly
di as result "Y2_NU0_MEAN=" %21.15g r(mean)
quietly margins, dydx(x) predict(nu0)
matrix list r(b), format(%21.15g)
matrix list r(V), format(%21.15g)
estimates store y2_balanced

di as txt "=== CASE 1C: BALANCED JOINT SYSTEM, DEFAULT PANEL CLUSTERS ==="
suest2 y1_balanced y2_balanced
di as txt "N=" e(N) " N_clust=" e(N_clust) " N_g=" e(N_g)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as txt "=== CASE 2: PARTIAL OBSERVATION OVERLAP, HIGHER CLUSTERS ==="
quietly xtpoisson y1_partial x z, re
estimates store y1_partial
quietly xtpoisson y2_partial x z, re
estimates store y2_partial
suest2 y1_partial y2_partial, cluster(higher)
di as txt "N=" e(N) " N_clust=" e(N_clust) " N_g=" e(N_g)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as txt "=== CASE 3: DISJOINT PANELS ==="
quietly xtpoisson y1_left x z, re
estimates store y1_left
quietly xtpoisson y2_right x z, re
estimates store y2_right
suest2 y1_left y2_right
di as txt "N=" e(N) " N_clust=" e(N_clust) " N_g=" e(N_g)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as result "PANEL_POISSON_RE_CROSSLANG_BENCHMARK_COMPLETE=1"
log close panelrepoisson
