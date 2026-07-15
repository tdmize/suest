# Run documentation, tests, checks, and the website build in clean processes.
#
# From the repository root:
#   source("tools/preflight.R")

root <- normalizePath(getwd(), mustWork = TRUE)
required <- c(
  "brglm2", "callr", "devtools", "glm2", "ordinal", "pkgdown",
  "rlang", "roxygen2", "testthat"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing)) {
  stop(
    "Install these packages first: ",
    paste(missing, collapse = ", "),
    call. = FALSE
  )
}

if (utils::packageVersion("rlang") < "1.3.0") {
  stop(
    "Restart R, update rlang with install.packages('rlang'), restart R again, ",
    "and rerun this script.",
    call. = FALSE
  )
}

run_clean <- function(code) {
  callr::r(
    function(root, code) {
      setwd(root)
      eval(parse(text = code), envir = .GlobalEnv)
    },
    args = list(root = root, code = code),
    show = TRUE
  )
}

run_clean("devtools::document()")
run_clean("devtools::test()")
run_clean("devtools::check()")
run_clean("pkgdown::clean_site(); pkgdown::build_site(new_process = FALSE, install = TRUE, quiet = FALSE)")

message("All local package and website checks completed.")
