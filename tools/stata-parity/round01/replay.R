# Round 1 (October 3, 2026): joint coefficients and covariances of 26 systems
# (complementary log-log, Gamma, fractional, zero-inflated, tobit, interval,
# parametric survival, truncated, beta, and heteroskedastic probit pairs, and
# different-sample linear, logit, cloglog, and tobit pairs). Stata's suest2
# results are in parity_stata_v01.csv (exported by test_suest_r_stata_v01.do).
source("helpers.R")
source("parity_reference_v01.R")
source("compare_parity_v02.R")
# Weibull is left out: suest2 1.1.0 used the wrong bread for streg Weibull
# fits (fixed in suest2; round 19 checks Weibull again).
write_gaps(out[out$case != "streg_weibull", ], c("max_b_rel", "loc_V_rel", "all_V_rel"))
