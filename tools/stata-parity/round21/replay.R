# Rounds 21 and 22 (November 24, 2026): average predictions, average marginal
# effects of x, and mecompare's change for a one-unit increase in x with its
# cross-model difference, for 38 pairs of models through suest2 (round 20
# stopped at an error in the do-file). Standard errors are not stored for
# xtpoisson, re (Stata repeats the panel's dispersion score on every
# observation) or for hetprobit's effects (Stata's margins differs from a
# by-hand delta-method calculation that R matches). After ivprobit, suest2
# gives only the linear prediction, so the probability and its effect are
# checked against round 22's ivprobit, vce(robust) fits with margins,
# predict(pr fixedasf). The fixed-effects Poisson predictions are left out:
# Stata's leave out the fixed effects.
source("helpers.R")
source("make_v20.R")
source("compare_v21.R")
x$check <- paste(x$case, x$stata_stat, x$r_model, x$r_outcome)
drop <- x$case == "ivprobit" | (x$case == "xtpoisson_fe" & grepl("^pred", x$stata_stat))
no_se <- x$case == "xtpoisson_re" | (x$case == "hetprobit" & x$stata_stat == "dydx")
write_gaps(x[!drop, ], "est_diff", check = "check")
write_gaps(x[!drop & !no_se, ], "se_rel", check = "check", append = TRUE)

# ivprobit: round 22's SE|iv_<model>_<mean or dydx>_pr_fixedasf lines
se <- grep("^SE\\|iv_[AB]_(mean|dydx)_pr_fixedasf\\|", readLines("test_suest_r_stata_v22.log"),
  value = TRUE)
p <- strsplit(se, "|", fixed = TRUE)
iv <- data.frame(model = substr(sapply(p, `[`, 2), 4, 4),
  stat = ifelse(grepl("_mean_", sapply(p, `[`, 2)), "pred", "dydx"),
  stata_est = as.numeric(sapply(p, `[`, 4)), stata_se = as.numeric(sapply(p, `[`, 5)))
ref <- read.csv("v20_r_reference.csv")
r <- ref[ref$case == "ivprobit", ]
iv <- merge(iv, r, by = c("stat", "model"))
iv$check <- paste("ivprobit pr fixedasf", iv$stat, iv$model)
iv$est_diff <- abs(iv$stata_est - iv$estimate)
iv$se_rel <- abs(iv$stata_se/iv$se - 1)
print(iv[, c("check", "stata_est", "estimate", "stata_se", "se")], digits = 8, row.names = FALSE)
if (nrow(iv) != 4) stop("expected 4 ivprobit rows from round 22")
write_gaps(iv, c("est_diff", "se_rel"), check = "check", append = TRUE)
