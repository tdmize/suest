# Run here after parsing the returned balanced log; optional comparison CSV.
args <- commandArgs(trailingOnly = TRUE)
f <- readRDS("../../tests/testthat/fixtures/glmmtmb-logit-rs-stata-balanced-v1.rds")
r <- readRDS("glmmtmb_logit_rs_r_reference.rds")
compare <- function(s, r, label) {
  stopifnot(identical(names(s$coefficients), names(r$coefficients)))
  data.frame(model = label, max_abs_coefficient = max(abs(s$coefficients-r$coefficients)),
    max_abs_covariance = max(abs(s$vcov-r$vcov)),
    max_relative_se = max(abs(sqrt(diag(s$vcov)/diag(r$vcov))-1)),
    abs_logLik = abs(s$logLik-r$logLik))
}
comparison <- do.call(rbind, lapply(names(f$components), function(y)
  compare(f$components[[y]], r$components[[y]], y)))
rj <- r$systems$balanced; rj$logLik <- sum(vapply(r$components[c("y1", "y2")], `[[`, numeric(1), "logLik"))
comparison <- rbind(comparison, compare(f$joint, rj, "balanced_joint_raw"))
print(comparison, digits = 10, row.names = FALSE)
for (y in names(f$components)) {
  s <- f$components[[y]]$native_margins; a <- r$components[[y]]$native_margins
  for (i in seq_along(s)) cat(y, names(s)[i], "estimate_difference=", s[[i]]["estimate"]-a$estimates[i],
    "relative_SE_difference=", sqrt(s[[i]]["variance"]/a$vcov[i, i])-1, "\n")
}
if (length(args)) write.csv(comparison, args[1L], row.names = FALSE)
