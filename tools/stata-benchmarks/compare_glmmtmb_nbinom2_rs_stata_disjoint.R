# Run here after parsing the returned disjoint log; optional comparison CSV.
args <- commandArgs(trailingOnly = TRUE)
f <- readRDS("../../tests/testthat/fixtures/glmmtmb-nbinom2-rs-stata-disjoint-v3.rds")
r <- readRDS("glmmtmb_nbinom2_rs_r_reference.rds")
compare <- function(s, r, label) {
  stopifnot(identical(names(s$coefficients), names(r$coefficients)))
  data.frame(model = label, max_abs_coefficient = max(abs(s$coefficients-r$coefficients)),
    max_abs_covariance = max(abs(s$vcov-r$vcov)),
    max_relative_se = max(abs(sqrt(diag(s$vcov)/diag(r$vcov))-1)),
    abs_logLik = abs(s$logLik-r$logLik))
}
comparison <- do.call(rbind, lapply(names(f$components), function(y)
  compare(f$components[[y]], r$components[[y]], y)))
rj <- r$systems$disjoint; rj$logLik <- sum(vapply(r$components[c("y1_left", "y2_right")], `[[`, numeric(1), "logLik"))
comparison <- rbind(comparison, compare(f$joint, rj, "disjoint_joint_raw"))
aligned <- f$joint; aligned$vcov <- aligned$vcov*99/98
comparison <- rbind(comparison, compare(aligned, rj, "disjoint_joint_normalized"))
cat("STATA_RAW_MAX_CROSS_BLOCK=", max(abs(f$joint$vcov[1:7, 8:14])), "\n")
cat("R_MAX_CROSS_BLOCK=", max(abs(rj$vcov[1:7, 8:14])), "\n")
print(comparison, digits = 10, row.names = FALSE)
for (y in names(f$components)) {
  s <- f$components[[y]]$native_margins; a <- r$components[[y]]$native_margins
  for (i in seq_along(s)) cat(y, names(s)[i], "estimate_difference=", s[[i]]["estimate"]-a$estimates[i],
    "relative_SE_difference=", sqrt(s[[i]]["variance"]/a$vcov[i, i])-1, "\n")
}
if (length(args)) write.csv(comparison, args[1L], row.names = FALSE)
