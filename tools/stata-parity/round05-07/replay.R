# Rounds 5-7 (October 4-5, 2026): a multiple-imputation system, linear
# random-intercept models (lmer, Stata mixed) alone and with a logit, survey
# systems (linear, logit, and Poisson; a logit on a domain), gologit2 pairs,
# and Heckman with a linear model. Cases where the two programs fit different
# models (Laplace against adaptive quadrature, for example) are left out.
source("helpers.R")
source("reference_v05_v07.R")
source("compare_v05_v07.R")
write_gaps(res, c("max_b", "max_se_rel", "max_V_rel"))
