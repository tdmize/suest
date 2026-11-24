# Rounds 14 and 17 (November 8, 2026): constant-dispersion negative binomial
# models (glmmTMB nbinom1; Stata nbreg, dispersion(constant)), alone and with
# the mean-dispersion negative binomial and Poisson; GEE pairs with ar1 and
# unstructured working correlations, which Stata fits with R's correlation
# fixed (round 17: the population-averaged xtreg, xtlogit, and xtpoisson fits
# with suest2); and gologit2 pairs with cauchit and log-log links.
source("helpers.R")
source("reference_v14.R")
source("compare_v14.R")
write_gaps(exact, c("max_b", "max_se_rel", "max_V_rel"))
