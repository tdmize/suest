# Compare the Stata logs of the zero-inflated rounds (ME| lines) with
# zeroinfl_r_reference.csv: round 1 for every prediction type except the
# count part, which round 4 checked (exp(xb) of the count equation). From
# compare_zeroinfl_v01.R. Writes comparison_zeroinfl.csv.
me <- function(file) {
  x <- grep("^ME\\|", readLines(file), value = TRUE)
  p <- strsplit(x, "|", fixed = TRUE)
  do.call(rbind, lapply(p, function(p) {
    n <- length(p)
    data.frame(case = paste(p[2:(n - 3)], collapse = "|"),
      stata_est = as.numeric(p[n - 1]), stata_se = as.numeric(p[n]))
  }))
}
v01 <- me("test_zeroinfl_types_stata_v01.log")
v04 <- me("test_zeroinfl_types_stata_v04.log")
st <- rbind(v01[!grepl("^(mean|dydx)_ir\\|", v01$case), ],
  v04[grepl("^(mean|dydx)_count\\|", v04$case), ])
ref <- read.csv("zeroinfl_r_reference.csv")
m <- merge(ref, st, by = "case", all.x = TRUE, sort = FALSE)
m$rel_est <- abs(m$estimate - m$stata_est) / abs(m$stata_est)
m$abs_est <- abs(m$estimate - m$stata_est)
m$rel_se <- abs(m$std.error - m$stata_se) / m$stata_se
print(m, digits = 6, row.names = FALSE)
write.csv(m, "comparison_zeroinfl.csv", row.names = FALSE)
