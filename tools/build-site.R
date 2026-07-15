# Build documentation and the pkgdown site in clean R subprocesses.
#
# From the repository root:
#   source("tools/build-site.R")

root <- normalizePath(getwd(), mustWork = TRUE)
required <- c("callr", "devtools", "pkgdown", "rlang", "roxygen2")
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
    "pkgdown requires a newer rlang installation. Restart R, run ",
    "install.packages('rlang'), restart R again, and rerun this script.",
    call. = FALSE
  )
}

callr::r(
  function(root) {
    setwd(root)
    devtools::document()
  },
  args = list(root = root),
  show = TRUE
)

callr::r(
  function(root) {
    setwd(root)
    pkgdown::clean_site()
    pkgdown::build_site(new_process = FALSE, install = TRUE, quiet = FALSE)
  },
  args = list(root = root),
  show = TRUE
)

message("Site built at: ", file.path(root, "docs", "articles", "suest.html"))
