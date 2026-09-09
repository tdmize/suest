# Compare the final returned Stata diagnostic with the R adapters.
args <- commandArgs(trailingOnly = TRUE)
workspace <- if (length(args)) normalizePath(args[1]) else normalizePath("../..")
pkgload::load_all(quiet = TRUE)

read_stata_matrices <- function(file) {
  lines <- readLines(file, warn = FALSE)
  starts <- grep("^e\\(b\\)\\[1,", lines)
  lapply(seq_along(starts), function(j) {
    end <- if (j < length(starts)) starts[j + 1L] - 1L else length(lines)
    block <- lines[starts[j]:end]
    vi <- grep("^symmetric e\\(V\\)\\[", block)[1L]
    b <- unlist(lapply(grep("^y1 +", block[seq_len(vi)], value = TRUE),
      function(x) as.numeric(strsplit(trimws(sub("^y1", "", x)), " +")[[1L]])))
    rows <- list()
    for (line in block[(vi + 1L):length(block)]) {
      token <- strsplit(trimws(line), " +")[[1L]]
      if (length(token) < 2L || !grepl(":.+", token[1L])) next
      values <- suppressWarnings(as.numeric(token[-1L]))
      if (anyNA(values)) next
      rows[[token[1L]]] <- c(rows[[token[1L]]], values)
    }
    stopifnot(length(rows) == length(b), identical(lengths(rows),
      setNames(seq_along(b), names(rows))))
    V <- matrix(0, length(b), length(b), dimnames = list(names(rows), names(rows)))
    for (i in seq_along(b)) V[i, seq_len(i)] <- rows[[i]]
    V <- V + t(V) - diag(diag(V))
    list(b = setNames(b, names(rows)), V = V)
  })
}

refs <- read_stata_matrices(file.path(workspace, "upload",
  "suest_r_followup_diagnostic.log"))
dat <- read.csv(file.path(workspace, "suest_r_panel_fe_crosslang_benchmark.csv"))
dat$time <- factor(dat$time)

compare <- function(label, b, V, reference) {
  keep <- diag(reference$V) != 0
  reference$b <- reference$b[keep]
  reference$V <- reference$V[keep, keep, drop = FALSE]
  result <- data.frame(case = label,
    max_b = max(abs(b - reference$b)),
    max_V = max(abs(V - reference$V)),
    max_relative_diagonal = max(abs(diag(V)/diag(reference$V) - 1)))
  print(result, row.names = FALSE, digits = 9)
  result
}

results <- list()
for (j in seq_along(c("y2", "y2_partial"))) {
  outcome <- c("y2", "y2_partial")[j]
  model <- nlme::lme(reformulate(c("x", "z", "time"), outcome),
    random = ~1 | id, data = dat, method = "ML", na.action = na.omit)
  component <- suest:::.suest_lme_ml_components(model)
  order <- c("x", "z", paste0("time", 2:5), "(Intercept)", "sigma_u", "sigma_e")
  results[[paste("ML native", outcome)]] <- compare(paste("ML native", outcome),
    component$parameters[order], component$bread[order, order]/nobs(model), refs[[j]])
}

models <- lapply(c("y1_partial", "y2_partial"), function(outcome)
  plm::plm(reformulate(c("x", "z"), outcome), data = dat,
    index = c("id", "time"), model = "random"))
fit <- suest(models[[1]], models[[2]], model_names = c("First", "Second"),
  cluster = "higher")
order <- unlist(lapply(c("First", "Second"), function(x)
  paste0(x, "::", c("x", "z", "(Intercept)"))))
results[["unbalanced RE simple"]] <- compare("unbalanced RE simple",
  coef(fit)[order], vcov(fit)[order, order], refs[[3L]])
saveRDS(results, file.path(workspace, "followup_comparison_results.rds"))
