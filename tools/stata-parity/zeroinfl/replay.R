# Zero-inflated rounds 1 and 4 (November 8, 2026): average predictions and
# effects of ment for each prediction type of zero-inflated Poisson and
# negative binomial models (overall mean, count part, probability of a
# structural zero, Pr(y = k)), and the cross-model differences.
source("helpers.R")
source("zeroinfl_reference.R")
source("compare_zeroinfl.R")
write_gaps(m, c("abs_est", "rel_se"))
