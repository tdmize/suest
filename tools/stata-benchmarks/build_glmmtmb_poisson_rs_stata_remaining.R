# Run here with the returned v5 log; an optional second argument sets output.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1L)
source("poisson_rs_joint_log_helpers.R")
lines <- readLines(args[1], warn = FALSE)
stopifnot("BENCHMARK_FILE_REVISION=5" %in% lines,
  "POISSON_RS_REMAINING_JOINT_MODELS_COMPLETE=1" %in% lines)
md5 <- unname(tools::md5sum(args[1]))
fixture <- list(
  partial_higher = read_poisson_rs_joint(lines, "partial", c("y1_partial", "y2_partial"), 20L, 5L, args[1], md5),
  disjoint = read_poisson_rs_joint(lines, "disjoint", c("y1_left", "y2_right"), 100L, 5L, args[1], md5))
stopifnot(all(vapply(fixture, function(x) !is.null(x$modelbased), logical(1))))
output <- if (length(args) >= 2L) args[2] else
  "../../tests/testthat/fixtures/glmmtmb-poisson-rs-stata-remaining-v5.rds"
saveRDS(fixture, output)
cat("POISSON_RS_REMAINING_STATA_FIXTURE_RECORDED=1\n")
