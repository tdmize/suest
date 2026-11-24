# Compare test_suest_r_stata_v14.log and test_suest_r_stata_v17.log with
# batch3_r_reference_v02.csv (suest2 refused the xtgee fits of round 14; round
# 17 fit the same models with xtreg, xtlogit, and xtpoisson, pa).
# Each Stata parameter is a function of one R parameter (same value, minus it
# for log-log models and Stata's lnalpha), and covariances follow from those
# functions. From compare_v14.R (November 8, 2026): suest2 systems and each
# model's own robust fit (GEE: through suest2 only). Writes comparison_v14.csv.
lines <- c(readLines("test_suest_r_stata_v14.log"), readLines("test_suest_r_stata_v17.log"))
field <- function(x, k) vapply(strsplit(x, "|", fixed = TRUE), function(p)
  if (length(p) >= k) p[k] else "", "")
ev <- grep("^EV\\|", lines, value = TRUE)
stata <- data.frame(case = field(ev, 2), kind = field(ev, 3), row = field(ev, 4),
  col = field(ev, 5), value = suppressWarnings(as.numeric(field(ev, 6))), stringsAsFactors = FALSE)

ref <- read.csv("batch3_r_reference_v02.csv", stringsAsFactors = FALSE)
reference <- function(case) {
  r <- ref[ref$case == case, ]
  if (!nrow(r)) stop("no R reference ", case)
  b <- setNames(r$value[r$kind == "b"], r$row[r$kind == "b"])
  V <- matrix(NA, length(b), length(b), dimnames = list(names(b), names(b)))
  v <- r[r$kind == "V", ]
  V[cbind(v$row, v$col)] <- v$value
  list(b = b, V = V)
}
transforms <- list(
  id = list(f = function(x) x, d = function(x) 1),
  neg = list(f = function(x) -x, d = function(x) -1))

# Kind of each model in each system.
types <- list(
  nb1_pair = c(A = "nb1", B = "nb1"), nb_three = c(A = "nb1", B = "nb2", C = "glm"),
  gee_ar1 = c(A = "glm", B = "glm"), gee_un = c(A = "glm", B = "glm"),
  gee_mixed = c(A = "glm", B = "glm"), gee_short = c(A = "glm", B = "glm"),
  cauchit_pair = c(A = "gologit", B = "gologit"), loglog_pair = c(A = "loglog", B = "loglog"),
  three_mixed = c(A = "gologit", B = "nb1", C = "glm"))

# Split a Stata name into model letter, equation, and variable. Joint suest2
# names carry the letter ("A_cnt1:x", "A_lndelta:_cons", "/A:lndelta");
# a single model's names do not ("cnt1:x", "/:lndelta"), so letter is given.
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

# Stata names (in order) -> R names and transformations; NULL drops a name.
map_names <- function(names, system, letter = NULL) {
  parts <- lapply(names, split_name, letter = letter)
  # gologit2 equations, numbered in order of appearance within each model
  eq_index <- list()
  for (p in parts) if (nzchar(p$eq) && !p$eq %in% eq_index[[p$model]])
    eq_index[[p$model]] <- c(eq_index[[p$model]], p$eq)
  lapply(parts, function(p) {
    type <- types[[system]][[p$model]]
    if (is.null(type)) stop("no model ", p$model, " in ", system)
    if (startsWith(p$var, "o.") || grepl("^[0-9]+b?\\.", p$var)) return(NULL)
    term <- if (p$var == "_cons") "(Intercept)" else p$var
    ancillary <- if (p$var %in% c("lndelta", "lnalpha")) p$var else
      if (p$eq %in% c("lndelta", "lnalpha")) p$eq else ""
    out <- function(r, tf = "id") list(r = paste0(p$model, "::", r), tf = tf)
    if (ancillary == "lndelta") return(out("log_phi"))
    if (ancillary == "lnalpha") return(out("ln_theta", "neg"))
    if (type %in% c("gologit", "loglog")) {
      tf <- if (type == "loglog") "neg" else "id"
      j <- match(p$eq, eq_index[[p$model]])
      # z is parallel (pl(z)) wherever it appears in these models
      return(if (term == "z") out("z", tf) else out(paste0(term, ":", j), tf))
    }
    out(term)
  })
}

compare <- function(scase, rcase, label = scase, letter = NULL) {
  system <- sub("_native_[A-Z]$", "", scase)
  s <- stata[stata$case == scase, ]
  if (!nrow(s)) stop(label, ": no Stata output")
  r <- reference(rcase)
  sb <- s[s$kind == "b", ]; sV <- s[s$kind == "V", ]
  maps <- map_names(sb$row, system, letter)
  keep <- !vapply(maps, is.null, NA)
  maps <- maps[keep]; sn <- sb$row[keep]
  rn <- vapply(maps, `[[`, "", "r"); tf <- vapply(maps, `[[`, "", "tf")
  if (!all(rn %in% names(r$b))) stop(label, ": no R parameter ",
    paste(setdiff(rn, names(r$b)), collapse = ", "))
  fb <- mapply(function(v, t) transforms[[t]]$f(v), r$b[rn], tf)
  J <- mapply(function(v, t) transforms[[t]]$d(v), r$b[rn], tf)
  Vr <- outer(J, J)*r$V[rn, rn, drop = FALSE]
  Vs <- matrix(sV$value[match(paste(rep(sn, each = length(sn)), rep(sn, length(sn))),
    paste(sV$row, sV$col))], length(sn), byrow = TRUE)
  data.frame(case = label, params = length(sn),
    max_b = max(abs(fb - sb$value[keep])),
    max_se_rel = max(abs(sqrt(diag(Vs))/sqrt(diag(Vr)) - 1)),
    max_V_rel = max(abs(Vs - Vr)/sqrt(outer(diag(Vr), diag(Vr)))))
}

exact <- list()
for (system in names(types))
  exact[[length(exact) + 1]] <- compare(system, system,
    paste(system, if (startsWith(system, "gee")) "suest2 (round 17)" else "suest2"))
for (system in c("nb1_pair", "cauchit_pair", "loglog_pair")) for (m in c("A", "B"))
  exact[[length(exact) + 1]] <- compare(paste0(system, "_native_", m), paste0(system, "_native_", m),
    paste(system, m, "alone, robust"), letter = m)
exact <- do.call(rbind, exact)
num <- c("max_b", "max_se_rel", "max_V_rel")
exact[num] <- lapply(exact[num], signif, 3)
options(width = 220)
print(exact, row.names = FALSE, right = FALSE)
write.csv(exact, "comparison_v14.csv", row.names = FALSE)
