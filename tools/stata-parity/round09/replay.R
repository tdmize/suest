# Round 9 (November 6, 2026): survey models through svy: with suest2 (ordered
# logit and probit, multinomial logit, negative binomial, generalized ordered
# logit, tobit, interval regression, Weibull, instrumental variables, a
# six-model system, and models on a domain), the joint system and each model's
# own svy: fit.
source("helpers.R")
source("survey_models_reference_v01.R")
source("compare_v09.R")
write_gaps(res, c("max_b", "max_se_rel", "max_V_rel"))
