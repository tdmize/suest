version 16
clear all
set more off

capture log close _all
log using suest_r_panel_logit_re_crosslang_benchmark.log, replace text name(panelrelogit)

di as result "BENCHMARK_FILE_REVISION=2"
di as txt "SUEST R RANDOM-EFFECTS PANEL LOGIT CROSS-LANGUAGE BENCHMARK"
di as txt "Stata version: `c(stata_version)'"
which suest2

import delimited using suest_r_panel_logit_re_crosslang_benchmark.csv, clear varnames(1)
xtset id time

di as txt "=== CASE 1A: BALANCED MODEL Y1, NONADAPTIVE 12-POINT QUADRATURE ==="
quietly xtlogit y1 x z, re intmethod(ghermite) intpoints(12)
di as txt "N=" e(N) " N_g=" e(N_g) " LL=" %21.15g e(ll)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
predict double y1_pr if e(sample), pr
predict double y1_xb if e(sample), xb
quietly summarize y1_pr if e(sample), meanonly
di as result "Y1_PR_MEAN=" %21.15g r(mean)
quietly margins, dydx(x) predict(pr)
matrix list r(b), format(%21.15g)
matrix list r(V), format(%21.15g)
estimates store y1_balanced_nonadaptive

di as txt "=== CASE 1B: BALANCED MODEL Y2, NONADAPTIVE 12-POINT QUADRATURE ==="
quietly xtlogit y2 x z, re intmethod(ghermite) intpoints(12)
di as txt "N=" e(N) " N_g=" e(N_g) " LL=" %21.15g e(ll)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
predict double y2_pr if e(sample), pr
predict double y2_xb if e(sample), xb
quietly summarize y2_pr if e(sample), meanonly
di as result "Y2_PR_MEAN=" %21.15g r(mean)
quietly margins, dydx(x) predict(pr)
matrix list r(b), format(%21.15g)
matrix list r(V), format(%21.15g)
estimates store y2_balanced_nonadaptive

di as txt "=== CASE 1C: BALANCED ADAPTIVE FITS REQUIRED BY SUEST2 ==="
quietly xtlogit y1 x z, re intmethod(mvaghermite) intpoints(12)
di as txt "Y1_ADAPTIVE_N=" e(N) " N_g=" e(N_g) " LL=" %21.15g e(ll)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
estimates store y1_balanced

quietly xtlogit y2 x z, re intmethod(mvaghermite) intpoints(12)
di as txt "Y2_ADAPTIVE_N=" e(N) " N_g=" e(N_g) " LL=" %21.15g e(ll)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
estimates store y2_balanced

di as txt "=== CASE 1D: BALANCED JOINT SYSTEM, DEFAULT PANEL CLUSTERS ==="
suest2 y1_balanced y2_balanced
di as txt "N=" e(N) " N_clust=" e(N_clust) " N_g=" e(N_g)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as txt "=== CASE 2: PARTIAL OBSERVATION OVERLAP, HIGHER CLUSTERS ==="
quietly xtlogit y1_partial x z, re intmethod(mvaghermite) intpoints(12)
estimates store y1_partial
quietly xtlogit y2_partial x z, re intmethod(mvaghermite) intpoints(12)
estimates store y2_partial
suest2 y1_partial y2_partial, cluster(higher)
di as txt "N=" e(N) " N_clust=" e(N_clust) " N_g=" e(N_g)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as txt "=== CASE 3: DISJOINT PANELS ==="
quietly xtlogit y1_left x z, re intmethod(mvaghermite) intpoints(12)
estimates store y1_left
quietly xtlogit y2_right x z, re intmethod(mvaghermite) intpoints(12)
estimates store y2_right
suest2 y1_left y2_right
di as txt "N=" e(N) " N_clust=" e(N_clust) " N_g=" e(N_g)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as result "PANEL_LOGIT_RE_CROSSLANG_BENCHMARK_COMPLETE=1"
log close panelrelogit
