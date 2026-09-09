version 16
clear all
set more off

capture log close _all
log using suest_r_cluster_crosslang_benchmark.log, replace text name(clustercross)

di as result "BENCHMARK_FILE_REVISION=1"
di as txt "SUEST R CLUSTER CROSS-LANGUAGE BENCHMARK"
di as txt "Stata version: `c(stata_version)'"
which suest2

import delimited using suest_r_cluster_crosslang_benchmark.csv, clear varnames(1)

di as txt "=== CASE 1: IDENTICAL-SAMPLE LM/LOGIT ==="
quietly regress y x z
estimates store c1_lm
quietly logit b x z
estimates store c1_logit
suest2 c1_lm c1_logit, cluster(cluster4)
di as txt "N=" e(N) " N_clust=" e(N_clust)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as txt "=== CASE 2: PARTIAL-OVERLAP LM/LOGIT ==="
quietly regress y x z if id <= 900
estimates store c2_lm
quietly logit b x z if id >= 301
estimates store c2_logit
suest2 c2_lm c2_logit, cluster(cluster4)
di as txt "N=" e(N) " N_clust=" e(N_clust)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as txt "=== CASE 3: DISJOINT ROWS, SHARED CLUSTERS ==="
quietly regress y x z if side == 0
estimates store c3_lm
quietly logit b x z if side == 1
estimates store c3_logit
suest2 c3_lm c3_logit, cluster(cluster2)
di as txt "N=" e(N) " N_clust=" e(N_clust)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as txt "=== CASE 4: IDENTICAL-SAMPLE CLUSTERED 2SLS ==="
quietly ivregress 2sls y_iv x (endogenous = z1)
estimates store c4_iv_base
quietly ivregress 2sls y_iv x (endogenous = z1 z2)
estimates store c4_iv_adjusted
suest2 c4_iv_base c4_iv_adjusted, cluster(cluster4)
di as txt "N=" e(N) " N_clust=" e(N_clust)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as txt "=== CASE 5: PWEIGHTED LM/LOGIT WITH CLUSTERS ==="
quietly regress y x z [pweight=pw]
estimates store c5_lm
quietly logit b x z [pweight=pw]
estimates store c5_logit
suest2 c5_lm c5_logit, cluster(cluster4)
di as txt "N=" e(N) " N_clust=" e(N_clust)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)
capture noisily suest2_cleanup, force

di as result "CLUSTER_CROSSLANG_BENCHMARK_COMPLETE=1"
log close clustercross
