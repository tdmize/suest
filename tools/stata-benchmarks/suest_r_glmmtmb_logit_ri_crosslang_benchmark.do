version 16
clear all
set more off

capture log close _all
log using suest_r_glmmtmb_logit_ri_crosslang_benchmark.log, replace text name(glmmlogit)

di as result "BENCHMARK_FILE_REVISION=3"
di as txt "SUEST R GLMMTMB RANDOM-INTERCEPT LOGIT CROSS-LANGUAGE BENCHMARK"
di as txt "Stata version: `c(stata_version)'"
which suest2

import delimited using suest_r_glmmtmb_logit_ri_crosslang_benchmark.csv, clear varnames(1)

di as txt "=== CASE 1A: BALANCED MODEL Y1, LAPLACE ==="
quietly melogit y1 x z || id:, intmethod(laplace)
di as txt "N=" e(N) " LL=" %21.15g e(ll)
matrix list e(N_g), format(%21.15g)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
predict double y1_mu if e(sample), mu marginal
quietly summarize y1_mu if e(sample), meanonly
di as result "Y1_MU_MARGINAL_MEAN=" %21.15g r(mean)
quietly margins, dydx(x) predict(mu marginal)
matrix list r(b), format(%21.15g)
matrix list r(V), format(%21.15g)
estimates store y1_laplace

di as txt "=== CASE 1B: BALANCED MODEL Y2, LAPLACE ==="
quietly melogit y2 x z || id:, intmethod(laplace)
di as txt "N=" e(N) " LL=" %21.15g e(ll)
matrix list e(N_g), format(%21.15g)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
predict double y2_mu if e(sample), mu marginal
quietly summarize y2_mu if e(sample), meanonly
di as result "Y2_MU_MARGINAL_MEAN=" %21.15g r(mean)
quietly margins, dydx(x) predict(mu marginal)
matrix list r(b), format(%21.15g)
matrix list r(V), format(%21.15g)
estimates store y2_laplace

di as txt "=== CASE 1C: EXPECTED SUEST2 LAPLACE REFUSAL ==="
capture noisily suest2 y1_laplace y2_laplace
di as result "SUEST2_LAPLACE_RC=" _rc
capture noisily suest2_cleanup, force

di as txt "=== CASE 1D: BALANCED ADAPTIVE 12-POINT FITS REQUIRED BY SUEST2 ==="
quietly melogit y1 x z || id:, intmethod(mvaghermite) intpoints(12)
di as txt "Y1_ADAPTIVE_N=" e(N) " LL=" %21.15g e(ll)
matrix list e(N_g), format(%21.15g)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
estimates store y1_adaptive

quietly melogit y2 x z || id:, intmethod(mvaghermite) intpoints(12)
di as txt "Y2_ADAPTIVE_N=" e(N) " LL=" %21.15g e(ll)
matrix list e(N_g), format(%21.15g)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
estimates store y2_adaptive

di as txt "=== CASE 1E: BALANCED JOINT SYSTEM, DEFAULT GROUP CLUSTERS ==="
suest2 y1_adaptive y2_adaptive
di as txt "N=" e(N) " N_clust=" e(N_clust)
capture noisily matrix list e(N_g), format(%21.15g)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as txt "=== CASE 2: PARTIAL OBSERVATION OVERLAP, HIGHER CLUSTERS ==="
quietly melogit y1_partial x z || id:, intmethod(mvaghermite) intpoints(12)
estimates store y1_partial
quietly melogit y2_partial x z || id:, intmethod(mvaghermite) intpoints(12)
estimates store y2_partial
suest2 y1_partial y2_partial, cluster(higher)
di as txt "N=" e(N) " N_clust=" e(N_clust)
capture noisily matrix list e(N_g), format(%21.15g)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as txt "=== CASE 3: DISJOINT GROUPS ==="
quietly melogit y1_left x z || id:, intmethod(mvaghermite) intpoints(12)
estimates store y1_left
quietly melogit y2_right x z || id:, intmethod(mvaghermite) intpoints(12)
estimates store y2_right
suest2 y1_left y2_right
di as txt "N=" e(N) " N_clust=" e(N_clust)
capture noisily matrix list e(N_g), format(%21.15g)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as txt "=== CASE 4A: BALANCED JOINT GSEM, LAPLACE ==="
quietly gsem (y1 <- x z M1[id], family(bernoulli) link(logit)) ///
  (y2 <- x z M2[id], family(bernoulli) link(logit)), ///
  cov(M1[id]*M2[id]@0) intmethod(laplace) vce(robust)
capture noisily di as txt "N=" e(N) " N_clust=" e(N_clust) " LL=" %21.15g e(ll)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
di as result "GSEM_LAPLACE_BALANCED_COMPLETE=1"

di as txt "=== CASE 4B: PARTIAL OBSERVATION OVERLAP, HIGHER-CLUSTER GSEM LAPLACE ==="
quietly gsem (y1_partial <- x z M1[id], family(bernoulli) link(logit)) ///
  (y2_partial <- x z M2[id], family(bernoulli) link(logit)), ///
  cov(M1[id]*M2[id]@0) intmethod(laplace) vce(cluster higher)
capture noisily di as txt "N=" e(N) " N_clust=" e(N_clust) " LL=" %21.15g e(ll)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
di as result "GSEM_LAPLACE_PARTIAL_COMPLETE=1"

di as txt "=== CASE 4C: DISJOINT-GROUP GSEM, LAPLACE ==="
quietly gsem (y1_left <- x z M1[id], family(bernoulli) link(logit)) ///
  (y2_right <- x z M2[id], family(bernoulli) link(logit)), ///
  cov(M1[id]*M2[id]@0) intmethod(laplace) vce(robust)
capture noisily di as txt "N=" e(N) " N_clust=" e(N_clust) " LL=" %21.15g e(ll)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
di as result "GSEM_LAPLACE_DISJOINT_COMPLETE=1"

di as result "GLMMTMB_LOGIT_RI_CROSSLANG_BENCHMARK_COMPLETE=1"
log close glmmlogit
