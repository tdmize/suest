# Compare test_suest_r_stata_v18.log with batch4_r_reference_v01.csv.
# Every Stata parameter equals one R parameter:
#   streg gompertz:  _t:x = x, _t:_cons = rate, gamma = shape,
#                    gamma:z = shape(z)
#   streg ggamma:    _t:x = x, _t:_cons = mu, lnsigma = sigma,
#                    lnsigma:z = sigma(z), kappa = Q, kappa:b = Q(b)
#   biprobit:        y1:x = y1:x, y1:_cons = y1:(Intercept), athrho = athrho
#   probit:          yp:x = x, yp:_cons = (Intercept)
#   Table 1: checks that should agree to rounding (suest2 systems and each
#            model's own robust fit).
#   Table 2: average predictions (margins) and log likelihoods.
# Usage: Rscript compare_v18.R [test_suest_r_stata_v18.log]
args <- commandArgs(TRUE)
logfile <- if (length(args)) args[1] else "test_suest_r_stata_v18.log"
lines <- readLines(logfile)
field <- function(x, k) vapply(strsplit(x, "|", fixed = TRUE), function(p)
  if (length(p) >= k) p[k] else "", "")
ev <- grep("^EV\\|", lines, value = TRUE)
stata <- data.frame(case = field(ev, 2), kind = field(ev, 3), row = field(ev, 4),
  col = field(ev, 5), value = suppressWarnings(as.numeric(field(ev, 6))), stringsAsFactors = FALSE)
il <- grep("^INFO\\|", lines, value = TRUE)
stata_info <- data.frame(case = field(il, 2), key = field(il, 3),
  value = suppressWarnings(as.numeric(field(il, 4))), stringsAsFactors = FALSE)
cat("Fits:\n"); writeLines(grep("^FIT\\|", lines, value = TRUE))
failed <- grep("^(EV|FIT)\\|.*(failed|no results)", lines, value = TRUE)
if (length(failed)) { cat("\nNot run:\n"); writeLines(failed) }

ref <- read.csv("batch4_r_reference_v01.csv", stringsAsFactors = FALSE)
reference <- function(case) {
  r <- ref[ref$case == case, ]
  if (!nrow(r)) stop("no R reference ", case)
  b <- setNames(r$value[r$kind == "b"], r$row[r$kind == "b"])
  V <- matrix(NA, length(b), length(b), dimnames = list(names(b), names(b)))
  v <- r[r$kind == "V", ]
  V[cbind(v$row, v$col)] <- v$value
  list(b = b, V = V)
}

# Kind of each model in each system.
types <- list(
  gomp_pair = c(A = "gompertz", B = "gompertz"), entry_pair = c(A = "gompertz", B = "gompertz"),
  ggamma_pair = c(A = "gengamma", B = "gengamma"),
  biprobit_pair = c(A = "biprobit", B = "biprobit"),
  three_mixed = c(A = "gompertz", B = "biprobit", C = "glm"))
ancillary <- c("gamma", "lnsigma", "kappa", "athrho")

# Split a Stata name into model letter, equation, and variable. Joint suest2
# names carry the letter ("A__t:x", "A_gamma:_cons", "/A:gamma"); a single
# model's names do not ("_t:x", "/:gamma"), so letter is given.
split_name <- function(name, letter) {
  name <- sub("^/", "", name)
  eq <- if (grepl(":", name)) sub(":.*$", "", name) else ""
  var <- sub("^[^:]*:", "", name)
  if (is.null(letter)) {
    model <- sub("_.*$", "", eq)
    eq <- if (grepl("_", eq)) sub("^[^_]*_", "", eq) else ""
  } else model <- letter
  list(model = model, eq = eq, var = var)
}

# Stata names (in order) -> R names
map_names <- function(names, system, letter = NULL) {
  vapply(names, function(name) {
    p <- split_name(name, letter)
    type <- types[[system]][[p$model]]
    if (is.null(type)) stop("no model ", p$model, " in ", system)
    term <- if (p$var == "_cons") "(Intercept)" else p$var
    # an ancillary parameter without covariates ("/:gamma") or with them
    # ("gamma:_cons", "gamma:z")
    anc <- if (p$var %in% ancillary) p$var else if (p$eq %in% ancillary) p$eq else ""
    base <- p$var %in% c(ancillary, "_cons")
    r <- switch(type,
      gompertz = if (anc == "gamma") (if (base) "shape" else paste0("shape(", p$var, ")")) else
        if (p$var == "_cons") "rate" else p$var,
      gengamma = if (anc == "lnsigma") (if (base) "sigma" else paste0("sigma(", p$var, ")")) else
        if (anc == "kappa") (if (base) "Q" else paste0("Q(", p$var, ")")) else
        if (p$var == "_cons") "mu" else p$var,
      biprobit = if (anc == "athrho") "athrho" else paste0(p$eq, ":", term),
      glm = term)
    paste0(p$model, "::", r)
  }, "", USE.NAMES = FALSE)
}

compare <- function(scase, rcase, label = scase, letter = NULL, system = NULL) {
  if (is.null(system)) system <- sub("_native_[A-Z]$", "", scase)
  s <- stata[stata$case == scase, ]
  if (!nrow(s)) return(data.frame(case = label, params = 0, max_b = NA, max_se_rel = NA,
    max_V_rel = NA, note = "no Stata output"))
  r <- reference(rcase)
  sb <- s[s$kind == "b", ]; sV <- s[s$kind == "V", ]
  rn <- map_names(sb$row, system, letter)
  if (!all(rn %in% names(r$b))) stop(label, ": no R parameter ",
    paste(setdiff(rn, names(r$b)), collapse = ", "))
  Vr <- r$V[rn, rn, drop = FALSE]
  sn <- sb$row
  Vs <- matrix(sV$value[match(paste(rep(sn, each = length(sn)), rep(sn, length(sn))),
    paste(sV$row, sV$col))], length(sn), byrow = TRUE)
  left <- setdiff(names(r$b), rn)
  data.frame(case = label, params = length(sn),
    max_b = max(abs(r$b[rn] - sb$value)),
    max_se_rel = max(abs(sqrt(diag(Vs))/sqrt(diag(Vr)) - 1)),
    max_V_rel = max(abs(Vs - Vr)/sqrt(outer(diag(Vr), diag(Vr)))),
    note = if (length(left)) paste("R not compared:", paste(left, collapse = " ")) else "")
}
safe <- function(...) tryCatch(compare(...), error = function(e)
  data.frame(case = list(...)[[1]], params = NA, max_b = NA, max_se_rel = NA, max_V_rel = NA,
    note = paste("not comparable:", conditionMessage(e))))

# Table 1
exact <- list()
for (system in names(types))
  exact[[length(exact) + 1]] <- safe(system, system, paste(system, "suest2"))
for (system in c("gomp_pair", "ggamma_pair", "entry_pair", "biprobit_pair")) for (m in c("A", "B"))
  exact[[length(exact) + 1]] <- safe(paste0(system, "_native_", m), paste0(system, "_native_", m),
    paste(system, m, "alone, robust"), letter = m)
exact <- do.call(rbind, exact)

# Table 2: average predictions and log likelihoods
margins_row <- function(case, label) {
  s <- stata[stata$case == case, ]
  r <- reference(case)
  if (!nrow(s)) return(data.frame(check = label, stata = NA, r = r$b[[1]], abs_diff = NA,
    se_stata = NA, se_r = sqrt(r$V[1, 1]), se_rel = NA))
  b <- s$value[s$kind == "b"]; se <- sqrt(s$value[s$kind == "V"])
  data.frame(check = label, stata = b, r = r$b[[1]], abs_diff = abs(b - r$b[[1]]),
    se_stata = se, se_r = sqrt(r$V[1, 1]), se_rel = abs(se/sqrt(r$V[1, 1]) - 1))
}
ll_row <- function(scase, rcase, label) {
  s <- stata_info$value[stata_info$case == scase & stata_info$key == "ll"]
  r <- ref$value[ref$case == rcase & ref$row == "ll"]
  data.frame(check = label, stata = if (length(s)) s else NA, r = r,
    abs_diff = if (length(s)) abs(s - r) else NA, se_stata = NA, se_r = NA, se_rel = NA)
}
other <- rbind(
  margins_row("gomp_A_margins", "Gompertz A: average median survival time"),
  margins_row("ggamma_B_margins", "Generalized gamma B: average median survival time"),
  margins_row("biprobit_A_margins", "Bivariate probit A: average Pr(y1 = 1, y2 = 1)"),
  ll_row("gomp_pair_native_A", "gA_info", "log likelihood: Gompertz A"),
  ll_row("gomp_pair_native_B", "gB_info", "log likelihood: Gompertz B"),
  ll_row("ggamma_pair_native_A", "qA_info", "log likelihood: generalized gamma A"),
  ll_row("ggamma_pair_native_B", "qB_info", "log likelihood: generalized gamma B"),
  ll_row("entry_pair_native_A", "eA_info", "log likelihood: Gompertz with entry A"),
  ll_row("entry_pair_native_B", "eB_info", "log likelihood: Gompertz with entry B"),
  ll_row("biprobit_pair_native_A", "bA_info", "log likelihood: bivariate probit A (GJRM's)"),
  ll_row("biprobit_pair_native_B", "bB_info", "log likelihood: bivariate probit B (GJRM's)"))

num <- function(d) { k <- vapply(d, is.numeric, NA); d[k] <- lapply(d[k], signif, 3); d }
options(width = 220)
cat("\nTable 1. Checks that should agree to rounding\n")
print(num(exact), row.names = FALSE, right = FALSE)
cat("\nTable 2. Average predictions (margins) and log likelihoods\n")
o <- other; o[c("abs_diff", "se_rel")] <- lapply(o[c("abs_diff", "se_rel")], signif, 3)
print(format(o, digits = 12), row.names = FALSE, right = FALSE)
write.csv(num(exact), "comparison_v18.csv", row.names = FALSE)
write.csv(other, "comparison_v18_other.csv", row.names = FALSE)
