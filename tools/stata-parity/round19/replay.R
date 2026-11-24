# Round 19 (November 8, 2026): median survival time predictions and the
# average effect of x on them for Weibull, exponential, lognormal, and
# log-logistic models (survreg; Stata streg, then margins).
source("helpers.R")
source("make_v19.R")
source("compare_v19.R")
write_gaps(x[!grepl("_ll$", x$case), ], c("rel_diff", "se_rel"))
write_gaps(x[grepl("_ll$", x$case), ], "rel_diff", append = TRUE)
