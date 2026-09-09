# Run from the repository root. Argument 1 is the workspace containing uploads/CSVs.
args <- commandArgs(trailingOnly = TRUE)
workspace <- if (length(args)) normalizePath(args[1]) else normalizePath("..")
pkgload::load_all(quiet = TRUE)

read_stata_matrices <- function(file) {
  lines <- readLines(file, warn = FALSE)
  starts <- grep("^e\\(b\\)\\[1,", lines)
  lapply(seq_along(starts), function(j) {
    end <- if (j < length(starts)) starts[j + 1] - 1L else length(lines)
    block <- lines[starts[j]:end]
    vi <- grep("^symmetric e\\(V\\)\\[", block)[1]
    b <- unlist(lapply(grep("^y1 +", block[seq_len(vi)], value = TRUE),
      function(x) as.numeric(strsplit(trimws(sub("^y1", "", x)), " +")[[1]])))
    rows <- list()
    for (line in block[(vi + 1):length(block)]) {
      token <- strsplit(trimws(line), " +")[[1]]
      if (length(token) < 2 || !grepl(":.+", token[1])) next
      values <- suppressWarnings(as.numeric(token[-1]))
      if (anyNA(values)) next
      rows[[token[1]]] <- c(rows[[token[1]]], values)
    }
    stopifnot(length(rows) == length(b), identical(lengths(rows),
      setNames(seq_along(b), names(rows))))
    V <- matrix(0, length(b), length(b), dimnames = list(names(rows), names(rows)))
    for (i in seq_along(b)) V[i, seq_len(i)] <- rows[[i]]
    V <- V + t(V) - diag(diag(V))
    list(b = setNames(b, names(rows)), V = V)
  })
}

comparisons <- list()
compare <- function(label, fit, reference, order) {
  keep <- diag(reference$V) != 0
  b <- reference$b[keep]
  V <- reference$V[keep, keep, drop = FALSE]
  rb <- coef(fit)[order]
  rV <- vcov(fit)[order, order, drop = FALSE]
  stopifnot(length(rb) == length(b), !anyNA(rb), identical(dim(V), dim(rV)))
  result <- data.frame(case = label, max_b = max(abs(rb - b)),
    max_V = max(abs(rV - V)), max_relative_diagonal = max(abs(diag(rV)/diag(V) - 1)))
  print(result, row.names = FALSE, digits = 7)
  comparisons[[label]] <<- list(summary = result, R = list(b = rb, V = rV),
    Stata = list(b = b, V = V), fit = fit)
  saveRDS(comparisons, file.path(workspace, "returned_log_comparisons.rds"))
}
getlog <- function(stem) read_stata_matrices(file.path(workspace, "upload",
  paste0("suest_r_", stem, "_crosslang_benchmark.log")))
panel <- read.csv(file.path(workspace, "suest_r_panel_fe_crosslang_benchmark.csv"))
panel$time <- factor(panel$time)
panel_refs <- list(fe = getlog("panel_fe"), bere = getlog("panel_be_re"), ml = getlog("panel_ml"))
for (route in c("within", "between", "random", "ML")) {
  for (partial in c(FALSE, TRUE)) {
    ys <- if (partial) c("y1_partial", "y2_partial") else c("y1", "y2")
    rhs <- if (route == "between") c("x", "z") else c("x", "z", "time")
    models <- tryCatch(lapply(ys, function(y) {
      f <- reformulate(rhs, y)
      if (route == "ML") nlme::lme(f, random = ~1 | id, data = panel,
        method = "ML", na.action = na.omit) else
        plm::plm(f, data = panel, index = c("id", "time"), model = route)
    }), error = function(e) { message(route, " ", partial, ": ", conditionMessage(e)); NULL })
    if (is.null(models)) next
    fit <- suest(models[[1]], models[[2]], model_names = c("First", "Second"),
      cluster = if (partial && route != "within") "higher" else NULL)
    local <- c("x", "z", if (route != "between") paste0("time", 2:5),
      "(Intercept)", if (route == "ML") c("sigma_u", "sigma_e"))
    order <- unlist(lapply(c("First", "Second"), function(x) paste0(x, "::", local)))
    ref <- switch(route, within = panel_refs$fe[[1 + partial]],
      between = panel_refs$bere[[1 + partial]], random = panel_refs$bere[[3 + partial]],
      ML = panel_refs$ml[[1 + partial]])
    compare(paste(route, partial), fit, ref, order)
    if (route == "within" && partial) {
      fit <- suest(models[[1]], models[[2]], model_names = c("First", "Second"), cluster = "higher")
      compare("within higher", fit, panel_refs$fe[[3]], order)
    }
  }
}
dat <- read.csv(file.path(workspace, "suest_r_gee_crosslang_benchmark.csv"))
refs <- getlog("gee")
families <- list(gaussian(), binomial(), binomial("probit"), binomial("cloglog"), poisson())
outcomes <- c("y", "binary", "binary", "binary", "count")
i <- 0L
for (correlation in c("independence", "exchangeable")) {
  for (j in seq_along(families)) {
    models <- lapply(list("x", c("x", "z")), function(rhs)
      geepack::geeglm(reformulate(rhs, outcomes[j]), id = id, data = dat,
        family = families[[j]], corstr = correlation,
        control = geepack::geese.control(epsilon = 1e-10, maxit = 100)))
    fit <- suest(models[[1]], models[[2]], model_names = c("Base", "Adjusted"))
    i <- i + 1L
    compare(paste("GEE", correlation, outcomes[j], families[[j]]$link), fit, refs[[i]],
      c("Base::x", "Base::(Intercept)", "Adjusted::x", "Adjusted::z", "Adjusted::(Intercept)"))
  }
}
dat <- read.csv(file.path(workspace, "suest_r_endogenous_crosslang_benchmark.csv"))
refs <- getlog("endogenous")
cache <- file.path(workspace, "endogenous_comparison_fits.rds")
if (file.exists(cache)) {
  fits <- readRDS(cache)
} else {
  bp1 <- mvProbit::mvProbit(cbind(binary1, binary2) ~ x, data = dat,
    algorithm = mvtnorm::TVPACK(), method = "BFGS", finalHessian = TRUE,
    iterlim = 300, reltol = 1e-12, printLevel = 0)
  bp2 <- mvProbit::mvProbit(cbind(binary1, binary2) ~ x + z, data = dat,
    algorithm = mvtnorm::TVPACK(), method = "BFGS", finalHessian = TRUE,
    iterlim = 300, reltol = 1e-12, printLevel = 0)
  fits <- list(bp = suest(bp1, bp2))
  saveRDS(fits, cache)
}
iv1 <- Rchoice::ivpml(binary_iv ~ x + endogenous | x + z1, data = dat,
  messages = FALSE, printLevel = 0, reltol = 1e-12)
iv2 <- Rchoice::ivpml(binary_iv ~ x + endogenous | x + z1 + z2, data = dat,
  messages = FALSE, printLevel = 0, reltol = 1e-12)
fits$iv <- suest(iv1, iv2)
fits$bp <- do.call(suest, fits$bp$models)
compare("biprobit", fits$bp, refs[[1]], c(2, 1, 4, 3, 5, 7, 8, 6, 10, 11, 9, 12))
compare("ivprobit", fits$iv, refs[[2]], c(3, 2, 1, 5, 6, 4, 8, 7, 11, 10, 9, 13, 14, 15, 12, 17, 16))
saveRDS(comparisons, file.path(workspace, "returned_log_comparisons.rds"))
