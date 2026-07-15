# Run the complete numerical acceptance suite against the package source.
#
# From the repository root:
#   source("tools/acceptance-tests/run_all_tests.R")

get_script_root <- function() {
  normalizePath(getwd(), mustWork = TRUE)
}

run_suest_acceptance_tests <- function() {
  root <- get_script_root()

  log_connection <- file(
    file.path(root, "suest_full_test_output.txt"),
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

  devtools::load_all(root, quiet = TRUE)

  source(file.path(
    root,
    "tools",
    "acceptance-tests",
    "tests",
    "helpers.R"
  ))

  test_initialize(root)
  require_test_packages()

  cat("SUEST FULL ACCEPTANCE-TEST LOG\n")
  cat("Started:", format(Sys.time()), "\n")
  cat("R version:", R.version.string, "\n")
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
  cat("MASS:", as.character(utils::packageVersion("MASS")), "\n")
  cat("nnet:", as.character(utils::packageVersion("nnet")), "\n")

  tests <- file.path(root, "tools", "acceptance-tests", "tests")
  source(file.path(tests, "test_paper_examples.R"))
  source(file.path(tests, "test_invariants.R"))
  source(file.path(tests, "test_comparison_universe.R"))
  source(file.path(tests, "test_validation.R"))

  results <- test_summary()
  utils::write.csv(
    results,
    file.path(root, "suest_full_test_results.csv"),
    row.names = FALSE
  )

  cat("\nCompleted:", format(Sys.time()), "\n")

  if (any(results$status == "FAIL"))
    stop(
      sum(results$status == "FAIL"),
      " SUEST acceptance tests failed.",
      call. = FALSE
    )

  invisible(results)
}

run_suest_acceptance_tests()
