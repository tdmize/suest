# Compare test_suest_r_stata_v21.log with v20_r_reference.csv (from make_v20.R).
# Changes from compare_v20.R: mecompare's rows for models with several outcomes
# are matched by model and outcome (Stata lists them outcome by outcome); bare
# margins after suest2 gives the probability for xtlogit, re and xtprobit, re
# (R's response), and the linear prediction for xtpoisson, re and fe (R's link),
# where mecompare's change is the change in the linear prediction (R's
# coefficient, dydx_link). Stata gives only the linear prediction after ivprobit
# (suest2), so ivprobit is listed apart. Writes comparison_v21.csv.
# Usage: Rscript compare_v21.R [test_suest_r_stata_v21.log]
args <- commandArgs(TRUE)
logfile <- if (length(args)) args[1] else "test_suest_r_stata_v21.log"
lines <- readLines(logfile)
fits <- grep("^FIT\\|", lines, value = TRUE)
cat("Fits that did not converge or failed:\n")
writeLines(fits[!grepl("converged=1|converged=\\.", fits)])
failed <- grep("^(MG|MC|EV)\\|.*failed", lines, value = TRUE)
if (length(failed)) { cat("\nNot run:\n"); writeLines(failed) }

parse <- function(tag) {
  x <- grep(paste0("^", tag, "\\|"), lines, value = TRUE)
  x <- x[!grepl("failed", x)]
  p <- strsplit(x, "|", fixed = TRUE)
  data.frame(case = vapply(p, `[`, "", 2), stat = vapply(p, `[`, "", 3),
    name = vapply(p, `[`, "", 4), estimate = as.numeric(vapply(p, `[`, "", 5)),
    se = as.numeric(vapply(p, `[`, "", 6)))
}
stata <- rbind(parse("MG"), parse("MC"))
stata$order <- ave(seq_len(nrow(stata)), stata$case, stata$stat, FUN = seq_along)
ref <- read.csv("v20_r_reference.csv")
ref$outcome[is.na(ref$outcome)] <- ""
ref$model[ref$stat == "change_diff"] <- "A-B"
ref$stat[ref$stat == "change_diff"] <- "change"

linear <- c("xtreg_fe", "xtreg_re", "xtreg_mle", "xtpoisson_re", "xtpoisson_fe")
# R statistic for each Stata statistic
r_stat <- function(case, stat) {
  base <- sub("_.*$", "", stat)
  if (case %in% c("xtpoisson_re", "xtpoisson_fe") && stat %in% c("pred", "dydx", "change"))
    return(c(pred = "pred_link", dydx = "dydx_link", change = "dydx_link")[[stat]])
  base
}
out <- list()
for (key in unique(paste(stata$case, stata$stat, sep = "\r"))) {
  s <- stata[paste(stata$case, stata$stat, sep = "\r") == key, ]
  case <- s$case[1]; stat <- s$stat[1]
  r <- ref[ref$case == case & ref$stat == r_stat(case, stat), ]
  m <- regmatches(stat, regexpr("_[AB]$", stat))
  if (length(m)) r <- r[r$model == sub("_", "", m), ]
  if (stat == "change") {
    # Stata names: x:A, x:B, x:Difference, or x:A_1, x:B_1, x:Difference_1, ...
    nm <- sub("^x:", "", s$name)
    s$model <- sub("_.*$", "", nm)
    s$model[s$model == "Difference"] <- "A-B"
    s$outcome <- ifelse(grepl("_", nm), sub("^[^_]*_", "", nm), "")
    if (r_stat(case, stat) == "dydx_link") s <- s[s$model != "A-B", ]
    r <- r[match(paste(s$model, s$outcome), paste(r$model, r$outcome)), ]
  } else {
    r <- r[seq_len(min(nrow(r), nrow(s))), ]
    s <- s[seq_len(nrow(r)), ]
  }
  out[[key]] <- data.frame(case = case, stata_stat = stat, r_stat = r_stat(case, stat),
    stata_name = s$name, r_model = r$model, r_outcome = r$outcome,
    stata_est = s$estimate, r_est = r$estimate, stata_se = s$se, r_se = r$se)
}
x <- do.call(rbind, out)
rownames(x) <- NULL
x$est_diff <- abs(x$stata_est - x$r_est)
x$est_rel <- x$est_diff / pmax(abs(x$r_est), 1e-8)
x$se_rel <- abs(x$stata_se / x$r_se - 1)
options(width = 250)
cat("\nLargest gaps by case and statistic (ivprobit: Stata's linear prediction, R's probability):\n")
agg <- aggregate(cbind(est_diff, est_rel, se_rel) ~ case + stata_stat, x, max, na.action = na.pass)
print(format(agg[order(agg$case == "ivprobit", agg$case), ], digits = 3), row.names = FALSE)
write.csv(x, "comparison_v21.csv", row.names = FALSE)
