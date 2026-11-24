# Compare test_suest_r_stata_v19.log with survreg_median_r_reference_v01.csv.
# Stata's streg log likelihood should exceed survreg's by the sum of log times
# over the failures (397.040394847 in this fixture).
# Usage: Rscript compare_v19.R [test_suest_r_stata_v19.log]
args <- commandArgs(TRUE)
logfile <- if (length(args)) args[1] else "test_suest_r_stata_v19.log"
lines <- readLines(logfile)
cat("Fits:\n"); writeLines(grep("^FIT\\|", lines, value = TRUE))
failed <- grep("^M\\|.*failed", lines, value = TRUE)
if (length(failed)) { cat("\nNot run:\n"); writeLines(failed) }
m <- grep("^M\\|", lines, value = TRUE)
m <- m[!grepl("failed", m)]
part <- strsplit(m, "|", fixed = TRUE)
stata <- data.frame(case = vapply(part, `[`, "", 2),
  estimate = as.numeric(vapply(part, `[`, "", 3)),
  se = suppressWarnings(as.numeric(vapply(part, function(p) if (length(p) >= 4) p[4] else NA_character_, ""))))
ref <- read.csv("survreg_median_r_reference_v01.csv")
d <- haven::read_dta("batch4_fixture_v01.dta")
shift <- sum(d$ev*log(d$tm))
x <- merge(ref, stata, by = "case", suffixes = c("_r", "_stata"), all.x = TRUE, sort = FALSE)
ll <- grepl("_ll$", x$case)
x$estimate_r[ll] <- x$estimate_r[ll] + shift
x$abs_diff <- abs(x$estimate_stata - x$estimate_r)
x$rel_diff <- x$abs_diff/abs(x$estimate_r)
x$se_rel <- abs(x$se_stata/x$se_r - 1)
cat(sprintf("\nR log likelihoods shifted by the sum of log times over the failures (%.9f)\n\n", shift))
options(width = 200)
print(format(x, digits = 12), row.names = FALSE, right = FALSE)
write.csv(x, "comparison_v19.csv", row.names = FALSE)
