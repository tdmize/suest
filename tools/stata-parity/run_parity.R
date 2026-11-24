# Reruns the R side of each returned Stata check with the installed suest and
# compares it with the Stata log stored with the round. A check fails when a
# gap between R and Stata is more than twice what it was when the round was
# stored, plus 1e-6.
# Usage, from the package folder with suest installed:
#   Rscript tools/stata-parity/run_parity.R                  every round
#   Rscript tools/stata-parity/run_parity.R round09 round18  some rounds
#   Rscript tools/stata-parity/run_parity.R --update round09 store the current gaps
# The twopart round needs meps.dta: it is downloaded unless SUEST_DATA names a
# folder that has it.
args <- commandArgs(TRUE)
update <- "--update" %in% args
args <- setdiff(args, "--update")
here <- normalizePath(file.path("tools", "stata-parity"), mustWork = TRUE)
rounds <- if (length(args)) args else
  basename(list.dirs(here, recursive = FALSE))[file.exists(file.path(list.dirs(here, recursive = FALSE), "replay.R"))]

run_round <- function(round) {
  src <- file.path(here, round)
  work <- file.path(tempdir(), round)
  unlink(work, recursive = TRUE)
  dir.create(work)
  file.copy(c(list.files(src, full.names = TRUE), list.files(file.path(here, "data"), full.names = TRUE),
    file.path(here, "helpers.R")), work)
  for (gz in list.files(work, "\\.gz$", full.names = TRUE)) {
    writeLines(readLines(gzfile(gz)), sub("\\.gz$", "", gz))
    unlink(gz)
  }
  start <- Sys.time()
  owd <- setwd(work)
  status <- system2(file.path(R.home("bin"), "Rscript"), "replay.R", stdout = "replay.log",
    stderr = "replay.log")
  setwd(owd)
  minutes <- as.numeric(difftime(Sys.time(), start, units = "mins"))
  gaps <- file.path(work, "gaps.csv")
  if (status != 0 || !file.exists(gaps)) {
    cat(sprintf("\n%s: the replay failed (status %s). Last lines of its log:\n", round, status))
    writeLines(tail(readLines(file.path(work, "replay.log")), 30))
    return(data.frame(round = round, check = "(replay)", measure = "", stored = NA, now = NA,
      limit = NA, pass = FALSE))
  }
  now <- read.csv(gaps, stringsAsFactors = FALSE)
  baseline <- file.path(src, "baseline.csv")
  if (update || !file.exists(baseline)) {
    file.copy(gaps, baseline, overwrite = TRUE)
    cat(sprintf("%s: stored %d gaps (%.1f minutes)\n", round, nrow(now), minutes))
  }
  stored <- read.csv(baseline, stringsAsFactors = FALSE)
  x <- merge(stored, now, by = c("check", "measure"), all.x = TRUE, suffixes = c("_stored", "_now"))
  x$limit <- 2*x$gap_stored + 1e-6
  x$pass <- is.na(x$gap_stored) | (!is.na(x$gap_now) & x$gap_now <= x$limit)
  cat(sprintf("%s: %d of %d gaps within their limits (%.1f minutes)\n", round, sum(x$pass), nrow(x),
    minutes))
  data.frame(round = round, check = x$check, measure = x$measure, stored = x$gap_stored,
    now = x$gap_now, limit = x$limit, pass = x$pass)
}

results <- do.call(rbind, lapply(rounds, run_round))
bad <- results[!results$pass, ]
if (nrow(bad)) {
  cat("\nGaps larger than their limits:\n")
  print(bad, row.names = FALSE)
  quit(status = 1)
}
cat("\nAll", nrow(results), "gaps are within their limits.\n")
