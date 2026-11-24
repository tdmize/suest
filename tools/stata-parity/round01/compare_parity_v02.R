# Compare Stata suest2 joint e(b)/e(V) (parity_stata_v01.csv) with the R
# reference (parity_r_reference_v01.csv). Location coefficients are compared
# directly; each ancillary parameter is mapped to Stata's scale by the
# transformation (identity, sign change, exp, exp(2x), log) that matches its
# estimate, and its covariance rows are transformed with the same Jacobian.
args <- commandArgs(trailingOnly = TRUE)
stata_file <- if (length(args) >= 1) args[1] else "parity_stata_v01.csv"
r_file <- if (length(args) >= 2) args[2] else "parity_r_reference_v01.csv"
S <- read.csv(stata_file, stringsAsFactors = FALSE)
R <- read.csv(r_file, stringsAsFactors = FALSE)

transforms <- list(
  identity = list(f = function(x) x, d = function(x) 1),
  negate   = list(f = function(x) -x, d = function(x) -1),
  exp      = list(f = exp, d = exp),
  exp2     = list(f = function(x) exp(2 * x), d = function(x) 2 * exp(2 * x)),
  log      = list(f = log, d = function(x) 1 / x)
)

map_case <- function(case) {
  sb <- S[S$case == case & S$kind == "b", ]
  rb <- R[R$case == case & R$kind == "b", ]
  if (!nrow(sb)) return(NULL)
  rows <- list()
  for (m in c("A", "B")) {
    s_names <- sb$row[startsWith(sb$row, paste0(m, "_")) | startsWith(sb$row, paste0("/", m, ":"))]
    slash <- startsWith(s_names, "/")
    eq <- ifelse(slash, sub("^/[^:]*:", "", s_names), sub(":.*$", "", sub(paste0("^", m, "_"), "", s_names)))
    var <- ifelse(slash, "_cons", sub("^.*:", "", s_names))
    loc_eq <- unique(eq[var == "x1" & eq != "inflate"])[1]
    r_names <- rb$row[startsWith(rb$row, paste0(m, "::"))]
    r_short <- sub(paste0("^", m, "::"), "", r_names)
    for (k in seq_along(s_names)) {
      v <- if (var[k] == "_cons") "(Intercept)" else var[k]
      target <- if (eq[k] == loc_eq) {
        c(v, paste0("count_", v))
      } else if (eq[k] == "inflate") {
        paste0("zero_", v)
      } else if (var[k] == "x3" && eq[k] != loc_eq) {
        paste0("het.", v)
      } else NA_character_
      hit <- if (all(is.na(target))) character(0) else intersect(target, r_short)
      if (length(hit)) {
        rows[[length(rows) + 1]] <- data.frame(stata = s_names[k], r = paste0(m, "::", hit[1]),
                                               type = "location", transform = "identity")
      } else {
        anc <- setdiff(r_short, c("(Intercept)", "x1", "x2", "x3",
          paste0("count_", c("(Intercept)", "x1", "x2", "x3")),
          paste0("zero_", c("(Intercept)", "x2")), "het.x3"))
        if (length(anc) != 1) next
        r_val <- rb$value[rb$row == paste0(m, "::", anc)]
        s_val <- sb$value[sb$row == s_names[k]]
        gaps <- vapply(transforms, function(t) suppressWarnings(abs(t$f(r_val) - s_val)), numeric(1))
        best <- names(which.min(gaps))
        rows[[length(rows) + 1]] <- data.frame(stata = s_names[k], r = paste0(m, "::", anc),
                                               type = "ancillary", transform = best)
      }
    }
  }
  do.call(rbind, rows)
}

result <- list()
for (case in unique(R$case)) {
  if (any(S$case == case & S$kind == "rc")) {
    result[[case]] <- data.frame(case = case, status = "Stata failed", n_match = NA,
      max_b_rel = NA, loc_V_rel = NA, all_V_rel = NA, mapped = NA)
    next
  }
  map <- map_case(case)
  if (is.null(map)) { result[[case]] <- data.frame(case = case, status = "missing in Stata",
    n_match = NA, max_b_rel = NA, loc_V_rel = NA, all_V_rel = NA, mapped = NA); next }
  rb <- R[R$case == case & R$kind == "b", ]; rv <- R[R$case == case & R$kind == "V", ]
  sb <- S[S$case == case & S$kind == "b", ]; sv <- S[S$case == case & S$kind == "V", ]
  bR <- setNames(rb$value, rb$row); bS <- setNames(sb$value, sb$row)
  VR <- matrix(rv$value, length(bR), dimnames = list(unique(rv$row), unique(rv$col)))[names(bR), names(bR)]
  VS <- matrix(NA, length(bS), length(bS), dimnames = list(names(bS), names(bS)))
  VS[cbind(sv$row, sv$col)] <- sv$value
  J <- vapply(seq_len(nrow(map)), function(i)
    transforms[[map$transform[i]]]$d(bR[[map$r[i]]]), numeric(1))
  bRt <- vapply(seq_len(nrow(map)), function(i)
    transforms[[map$transform[i]]]$f(bR[[map$r[i]]]), numeric(1))
  VRt <- diag(J, length(J)) %*% VR[map$r, map$r] %*% diag(J, length(J))
  VSm <- VS[map$stata, map$stata]
  scale <- sqrt(outer(diag(VSm), diag(VSm)))
  rel <- abs(VRt - VSm) / scale
  loc <- map$type == "location"
  nR <- R$value[R$case == case & R$kind == "N"]; nS <- S$value[S$case == case & S$kind == "N"]
  b_rel <- max(abs(bRt - bS[map$stata]) / sqrt(diag(VSm)))
  result[[case]] <- data.frame(case = case,
    status = if (nR != nS) sprintf("N differs (R %d, Stata %d)", nR, nS) else "compared",
    n_match = paste0(nrow(map), "/", length(bS)),
    max_b_rel = signif(b_rel, 3), loc_V_rel = signif(max(rel[loc, loc]), 3),
    all_V_rel = signif(max(rel), 3),
    mapped = paste(unique(paste(map$transform[!loc], sub("^/?A[_:]", "", map$stata[!loc]))), collapse = "; "))
}
out <- do.call(rbind, result)
options(width = 200)
cat("max_b_rel: largest coefficient gap in Stata SE units\n")
cat("loc_V_rel / all_V_rel: largest covariance gap relative to sqrt(V_ii V_jj), location block / all mapped parameters\n\n")
print(out, row.names = FALSE)
write.csv(out, "parity_comparison_v02.csv", row.names = FALSE)
