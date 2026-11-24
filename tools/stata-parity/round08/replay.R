# Round 8 (October 5, 2026): random-intercept logit, probit, complementary
# log-log, Poisson, ordered logit and probit models with 25-point adaptive
# quadrature in both programs (R: glmer and clmm; Stata: mvaghermite), alone
# and with ordinary models, and linear random-slope pairs (lmer and glmmTMB;
# Stata mixed).
source("helpers.R")
source("v08_reference_v01.R")
source("compare_v08.R")
write_gaps(out, c("max_b", "max_se_rel", "max_V_rel"))
