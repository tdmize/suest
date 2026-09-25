# Run here with a returned focused log; an optional second argument sets output.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1L)
source("poisson_rs_joint_log_helpers.R")
lines <- readLines(args[1], warn = FALSE)
revision <- as.integer(sub("BENCHMARK_FILE_REVISION=", "", grep("^BENCHMARK_FILE_REVISION=", lines, value = TRUE)))
stopifnot(length(revision) == 1L, revision %in% c(3L, 4L),
  any(c("POISSON_RS_BALANCED_WARMSTART_COMPLETE=1", "POISSON_RS_BALANCED_PRECISION_COMPLETE=1") %in% lines))
fixture <- read_poisson_rs_joint(lines, "balanced", c("y1", "y2"), 100L,
  revision, args[1], unname(tools::md5sum(args[1])))
output <- if (length(args) >= 2L) args[2] else
  paste0("../../tests/testthat/fixtures/glmmtmb-poisson-rs-stata-balanced-v", revision, ".rds")
saveRDS(fixture, output)
cat("POISSON_RS_BALANCED_STATA_FIXTURE_RECORDED=1\n")
