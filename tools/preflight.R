# Run package documentation, tests, checks, and the pkgdown build locally.
#
# From the repository root:
#   source("tools/preflight.R")

required <- c("devtools", "pkgdown", "roxygen2", "testthat")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing)) {
  stop(
    "Install the following packages first: ",
    paste(missing, collapse = ", "),
    call. = FALSE
  )
}

devtools::document()
devtools::test()
devtools::check()
pkgdown::build_site()

message("All local package and website checks completed.")
