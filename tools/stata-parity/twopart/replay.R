# Two-part round (November 23, 2026): suest(twopart = TRUE) against suest2's
# two-part systems on the MEPS data: logit or probit first parts with
# Poisson, Gamma, negative binomial, linear, and log (Duan and normal
# retransformation) second parts; average predictions and effects; subgroups;
# probability weights with and without weighted averaging; and a survey
# design. Needs meps.dta (downloaded unless SUEST_DATA names a folder with it).
source("helpers.R")
source("compare_twopart.R")
write_gaps(cmp, c("rel_est", "rel_se"), check = "check")
