# Run the complete numerical acceptance suite against the package source.
#
# From the repository root:
#   source("tools/acceptance-tests/run_all_tests.R")

get_script_root <- function() {
  normalizePath(getwd(), mustWork = TRUE)
}

run_suest_acceptance_tests <- function() {
  conflicts <- intersect(
    c("suest", "suest_newdata"),
    ls(envir = .GlobalEnv, all.names = TRUE)
  )
  if (length(conflicts))
    rm(list = conflicts, envir = .GlobalEnv)

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

  pkgload::load_all(root, quiet = TRUE)

  namespace <- asNamespace("suest")
  if (!exists(".suest_model_adapter", envir = namespace, inherits = FALSE))
    stop(
      "The loaded suest package does not contain the model-adapter layer.",
      call. = FALSE
    )
  if (!identical(environment(suest), namespace))
    stop(
      paste0(
        "A global object named 'suest' is masking the package function. ",
        "Restart R or remove the global object before testing."
      ),
      call. = FALSE
    )

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
  description <- read.dcf(file.path(root, "DESCRIPTION"))
  cat("suest version:", description[1L, "Version"], "\n")
  cat(
    "suest function environment:",
    environmentName(environment(suest)),
    "\n"
  )
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
  source(file.path(tests, "test_pweights.R"))
  source(file.path(tests, "test_pweights_extended.R"))
  source(file.path(tests, "test_invariants.R"))
  source(file.path(tests, "test_comparison_universe.R"))
  source(file.path(tests, "test_multiple_models.R"))
  source(file.path(tests, "test_glm_offsets.R"))
  source(file.path(tests, "test_model_adapters.R"))
  source(file.path(tests, "test_extended_models.R"))
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
