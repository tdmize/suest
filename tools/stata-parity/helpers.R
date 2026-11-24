# Shared by the replay scripts (run_parity.R copies this file next to them).

# Puts the stored fixture's values (the data Stata used) into the data frame a
# reference script simulated, keeping the script's column types. Factor
# columns are left as simulated.
use_fixture <- function(d, file) {
  f <- as.data.frame(haven::read_dta(file))
  if (nrow(f) != nrow(d)) stop(file, " has ", nrow(f), " rows; the script made ", nrow(d))
  for (v in intersect(names(f), names(d))) if (is.numeric(d[[v]])) {
    x <- as.vector(f[[v]])
    attributes(x) <- NULL
    d[[v]] <- if (is.integer(d[[v]])) as.integer(round(x)) else x
  }
  d
}

# Writes gaps.csv (or adds to it): one row per check and measure, with the
# absolute gap between R and Stata.
write_gaps <- function(x, measures, check = "case", keep = NULL, append = FALSE) {
  x <- as.data.frame(x)
  if (!is.null(keep)) {
    missing <- setdiff(keep, x[[check]])
    if (length(missing)) stop("checks not found: ", paste(missing, collapse = ", "))
    x <- x[x[[check]] %in% keep, , drop = FALSE]
  }
  out <- do.call(rbind, lapply(measures, function(m)
    data.frame(check = as.character(x[[check]]), measure = m, gap = abs(as.numeric(x[[m]])))))
  if (append) out <- rbind(read.csv("gaps.csv", stringsAsFactors = FALSE), out)
  write.csv(out, "gaps.csv", row.names = FALSE)
  invisible(out)
}

# Path to a data file from tdmize.github.io/data: from the folder SUEST_DATA
# names, if it has the file, or downloaded.
data_file <- function(name) {
  local <- file.path(Sys.getenv("SUEST_DATA"), name)
  if (nzchar(Sys.getenv("SUEST_DATA")) && file.exists(local)) return(local)
  path <- file.path(tempdir(), name)
  utils::download.file(paste0("https://tdmize.github.io/data/data/", name), path, mode = "wb",
    quiet = TRUE)
  path
}
