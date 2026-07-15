# Run only the alternative-model-engine acceptance tests.
#
# From the repository root:
#   source("tools/acceptance-tests/run_model_adapter_tests.R")

run_suest_model_adapter_tests <- function() {
  conflicts <- intersect(
    c("suest", "suest_newdata"),
    ls(envir = .GlobalEnv, all.names = TRUE)
  )
  if (length(conflicts))
    rm(list = conflicts, envir = .GlobalEnv)

  root <- normalizePath(getwd(), mustWork = TRUE)

  log_connection <- file(
    file.path(root, "suest_model_adapter_test_output.txt"),
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

  cat("SUEST MODEL-ADAPTER TEST LOG\n")
  cat("Started:", format(Sys.time()), "\n")
  cat("R version:", R.version.string, "\n")
  description <- read.dcf(file.path(root, "DESCRIPTION"))
  cat("suest version:", description[1L, "Version"], "\n")
  cat(
    "suest function environment:",
    environmentName(environment(suest)),
    "\n"
  )
  cat("glm2:", as.character(utils::packageVersion("glm2")), "\n")
  cat("ordinal:", as.character(utils::packageVersion("ordinal")), "\n")
  cat("brglm2:", as.character(utils::packageVersion("brglm2")), "\n")
  cat("sandwich:", as.character(utils::packageVersion("sandwich")), "\n")
  cat("marginaleffects:",
      as.character(utils::packageVersion("marginaleffects")), "\n")

  source(file.path(
    root,
    "tools",
    "acceptance-tests",
    "tests",
    "test_model_adapters.R"
  ))

  results <- test_summary()
  utils::write.csv(
    results,
    file.path(root, "suest_model_adapter_test_results.csv"),
    row.names = FALSE
  )

  cat("\nCompleted:", format(Sys.time()), "\n")

  if (any(results$status == "FAIL"))
    stop(
      sum(results$status == "FAIL"),
      " SUEST model-adapter tests failed.",
      call. = FALSE
    )

  invisible(results)
}

run_suest_model_adapter_tests()
