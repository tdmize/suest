# R side of the zero-inflated rounds (November 8, 2026): the prediction types
# of zero-inflated Poisson and negative binomial models (pscl::zeroinfl;
# Stata zip and zinb with suest2, then margins). From
# zeroinfl_types_reference_v01.R, reading zeroinfl_biochem.dta (pscl's
# bioChemists data, Long 1990), with the count part added as in round 4.
# Writes zeroinfl_r_reference.csv (case, estimate, std.error).
suppressPackageStartupMessages({library(suest); library(marginaleffects)})
dat <- as.data.frame(haven::read_dta("zeroinfl_biochem.dta"))
dat[] <- lapply(dat, as.numeric)
zip <- pscl::zeroinfl(art ~ factor(female) + kid5 + ment | ment, data = dat)
zinb <- pscl::zeroinfl(art ~ factor(female) + kid5 + ment | ment, data = dat, dist = "negbin")
fit <- suest(zip, zinb, model_names = c("ZIP", "ZINB"))

rows <- list()
add <- function(case, x) {
  x <- as.data.frame(x)
  rows[[length(rows) + 1L]] <<- data.frame(
    case = paste0(case, "|", x$group), estimate = x$estimate, std.error = x$std.error)
}
stata_names <- c(response = "n", count = "count", zero = "pr")
for (type in names(stata_names)) {
  add(paste0("mean_", stata_names[[type]]), avg_predictions(fit, type = type))
  add(paste0("dydx_", stata_names[[type]]), avg_slopes(fit, variables = "ment", type = type))
}
keep <- function(x) x[as.character(x$group) %in% paste0(rep(c("ZIP", "ZINB"), each = 3), "::", 0:2), ]
add("mean_prk", keep(avg_predictions(fit, type = "prob")))
add("dydx_prk", keep(avg_slopes(fit, variables = "ment", type = "prob")))
add("diff_mean_pr", avg_predictions(fit, type = "zero", hypothesis = ~ pairwise))
add("diff_dydx_pr", avg_slopes(fit, variables = "ment", type = "zero", hypothesis = ~ pairwise))
ref <- do.call(rbind, rows)
ref$case <- sub("\\|$", "", ref$case)
write.csv(ref, "zeroinfl_r_reference.csv", row.names = FALSE)
