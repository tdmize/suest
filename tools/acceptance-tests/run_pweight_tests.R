# Run focused pweight tests against the package source.
#
# From the repository root:
#   source("tools/acceptance-tests/run_pweight_tests.R")

get_script_root <- function() {
  normalizePath(getwd(), mustWork = TRUE)
}

run_suest_pweight_tests <- function() {
  conflicts <- intersect(
    c("suest", "suest_newdata"),
    ls(envir = .GlobalEnv, all.names = TRUE)
  )
  if (length(conflicts))
    rm(list = conflicts, envir = .GlobalEnv)

  root <- get_script_root()
  log_connection <- file(
    file.path(root, "suest_pweight_test_output.txt"),
    open = "wt"
  )
  sink(log_connection, split = TRUE)
  sink(log_connection, type = "message")

  on.exit({
    if (sink.number(type = "message") != 2)
      sink(type = "message")
    while (sink.number(type = "output") > 0)
      sink(type = "output")
    close(log_connection)
  }, add = TRUE)

  old_options <- options(warn = 1, width = 120, max.print = 200)
  on.exit(options(old_options), add = TRUE)

  pkgload::load_all(root, quiet = TRUE)

  source(file.path(
    root,
    "tools",
    "acceptance-tests",
    "tests",
    "helpers.R"
  ))
  test_initialize(root)
  require_test_packages()

  cat("SUEST PWEIGHT TEST LOG\n")
  cat("Started:", format(Sys.time()), "\n")
  cat("R version:", R.version.string, "\n")
  description <- read.dcf(file.path(root, "DESCRIPTION"))
  cat("suest version:", description[1L, "Version"], "\n")
  cat(
    "marginaleffects:",
    as.character(utils::packageVersion("marginaleffects")),
    "\n"
  )
  cat(
    "sandwich:",
    as.character(utils::packageVersion("sandwich")),
    "\n"
  )

  cat("\nTESTTHAT PWEIGHT TESTS\n")
  testthat::test_file(
    file.path(root, "tests", "testthat", "test-pweights.R"),
    reporter = "summary",
    stop_on_failure = TRUE,
    stop_on_warning = FALSE
  )

  source(file.path(
    root,
    "tools",
    "acceptance-tests",
    "tests",
    "test_pweights.R"
  ))

  results <- test_summary()
  utils::write.csv(
    results,
    file.path(root, "suest_pweight_test_results.csv"),
    row.names = FALSE
  )

  cat("\nCompleted:", format(Sys.time()), "\n")

  if (any(results$status == "FAIL"))
    stop(
      sum(results$status == "FAIL"),
      " pweight acceptance tests failed.",
      call. = FALSE
    )

  invisible(results)
}

run_suest_pweight_tests()
