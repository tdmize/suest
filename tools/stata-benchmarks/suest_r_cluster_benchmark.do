version 16
clear all
set more off
set seed 6109

capture log close _all
log using suest_r_cluster_benchmark.log, replace text name(clusterbench)

di as txt "SUEST R CLUSTER BENCHMARK"
di as txt "Stata version: `c(stata_version)'"
which suest2

set obs 1200
gen long id = _n
gen long cluster = ceil(id/4)
gen double cluster_u = rnormal() if mod(id - 1, 4) == 0
sort cluster id
by cluster (id): replace cluster_u = cluster_u[1]
gen double x = rnormal() + .25*cluster_u
gen double z = rnormal()
gen double y = 1 + .45*x - .30*z + .60*cluster_u + rnormal()
gen double pr = invlogit(-.25 + .50*x - .20*z + .45*cluster_u)
gen byte b = runiform() < pr

di as txt "=== CASE 1: IDENTICAL SAMPLES, OFFICIAL SUEST ==="
quietly regress y x z
estimates store c1_lm
quietly logit b x z
estimates store c1_logit
suest c1_lm c1_logit, cluster(cluster)
di as txt "N=" e(N) " N_clust=" e(N_clust) " df_r=" e(df_r)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)

di as txt "=== CASE 1: IDENTICAL SAMPLES, SUEST2 ==="
suest2 c1_lm c1_logit, cluster(cluster)
di as txt "N=" e(N) " N_clust=" e(N_clust) " df_r=" e(df_r)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)

di as txt "=== CASE 2: PARTIAL OVERLAP, OFFICIAL SUEST ==="
quietly regress y x z if id <= 900
estimates store c2_lm
quietly logit b x z if id >= 301
estimates store c2_logit
suest c2_lm c2_logit, cluster(cluster)
di as txt "N=" e(N) " N_clust=" e(N_clust) " df_r=" e(df_r)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)

di as txt "=== CASE 2: PARTIAL OVERLAP, SUEST2 ==="
suest2 c2_lm c2_logit, cluster(cluster)
di as txt "N=" e(N) " N_clust=" e(N_clust) " df_r=" e(df_r)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)

di as txt "=== CASE 3: DISJOINT ROWS, SHARED CLUSTERS, OFFICIAL SUEST ==="
replace cluster = ceil(id/2)
gen byte side = mod(id, 2)
quietly regress y x z if side == 0
estimates store c3_lm
quietly logit b x z if side == 1
estimates store c3_logit
suest c3_lm c3_logit, cluster(cluster)
di as txt "N=" e(N) " N_clust=" e(N_clust) " df_r=" e(df_r)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)

di as txt "=== CASE 3: DISJOINT ROWS, SHARED CLUSTERS, SUEST2 ==="
suest2 c3_lm c3_logit, cluster(cluster)
di as txt "N=" e(N) " N_clust=" e(N_clust) " df_r=" e(df_r)
matrix list e(b), format(%21.15g)
matrix list e(V), format(%21.15g)

di as result "CLUSTER_BENCHMARK_COMPLETE=1"
log close clusterbench
