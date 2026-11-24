# Compare the returned Stata logs v05-v07 (lines starting EV|) with the R
# references. Each Stata parameter is mapped to an R parameter and a
# transformation (Stata value = f(R value)); covariances use the Jacobian.
# From compare_v05_v07_v02.R (October 5, 2026), keeping the cases where R and
# Stata fit the same models. Writes comparison_v05_v07.csv.
suppressPackageStartupMessages({library(suest); library(haven)})

parse_ev <- function(file) {
  lines <- grep("^EV\\|", readLines(file), value = TRUE)
  parts <- strsplit(lines, "|", fixed = TRUE)
  data.frame(case = sapply(parts, `[`, 2), kind = sapply(parts, `[`, 3),
    row = sapply(parts, `[`, 4), col = sapply(parts, `[`, 5),
    value = as.numeric(sapply(parts, `[`, 6)), stringsAsFactors = FALSE)
}
transforms <- list(
  id = list(f = function(x) x, d = function(x) 1),
  exp2 = list(f = function(x) exp(2*x), d = function(x) 2*exp(2*x)),
  log = list(f = log, d = function(x) 1/x))

# Stata name -> R name and transformation.
map_one <- function(name, rnames, outcomes = character(), heckman = FALSE) {
  name <- sub("^/", "", name)
  eq <- sub(":.*$", "", name); var <- sub("^[^:]*:", "", name)
  if (eq == "" || !grepl(":", name)) { eq <- ""; var <- name }
  if (eq %in% outcomes) {
    model <- c("A", "B", "C")[match(eq, outcomes)]; suffix <- ""
  } else {
    model <- sub("_.*$", "", eq); suffix <- if (grepl("_", eq)) sub("^[^_]*_", "", eq) else ""
  }
  pick <- function(candidates, tf) {
    hit <- candidates[paste0(model, "::", candidates) %in% rnames]
    if (length(hit)) list(r = paste0(model, "::", hit[1]), tf = tf) else NULL
  }
  base <- if (var == "_cons") "(Intercept)" else var
  if (suffix == "lns1_1_1") return(c(pick("log_sigma", "id"), pick("sigma_u", "log"))[1:2])
  if (suffix == "lnsig_e") return(c(pick("log_sigma_e", "id"), pick("sigma_e", "log"))[1:2])
  if (var == "var(_cons[id])") return(pick("log_sigma", "exp2"))
  if (suffix == "lnvar") return(pick("lnvar", "id"))
  if (var %in% c("athrho", "lnsigma")) return(pick(var, "id"))
  if (heckman && model == "A" && suffix == "s") return(pick(paste0("selection:", base), "id"))
  if (heckman && model == "A" && suffix != "") return(pick(paste0("outcome:", base), "id"))
  if (grepl("^[0-9]+$", suffix)) return(pick(c(paste0(base, ":", suffix), base), "id"))
  pick(base, "id")
}

compare_case <- function(stata, case, rb, rV, outcomes = character(), heckman = FALSE,
                         label = case) {
  s <- stata[stata$case == case, ]
  if (!nrow(s)) stop(label, ": no Stata output")
  sb <- s[s$kind == "b", ]
  maps <- lapply(sb$row, map_one, rnames = names(rb), outcomes = outcomes, heckman = heckman)
  keep <- !vapply(maps, is.null, logical(1))
  sV <- s[s$kind == "V", ]
  diag_s <- setNames(sV$value[sV$row == sV$col], sV$row[sV$row == sV$col])
  keep <- keep & diag_s[sb$row] > 0
  maps <- maps[keep]; sn <- sb$row[keep]
  rn <- vapply(maps, `[[`, "", "r"); tf <- vapply(maps, `[[`, "", "tf")
  if (!all(rn %in% rownames(rV))) stop(label, ": missing R names ", paste(setdiff(rn, rownames(rV)), collapse = ", "))
  rval <- rb[rn]
  fb <- mapply(function(v, t) transforms[[t]]$f(v), rval, tf)
  J <- mapply(function(v, t) transforms[[t]]$d(v), rval, tf)
  Vr <- outer(J, J)*rV[rn, rn]
  Vs <- matrix(NA, length(sn), length(sn))
  for (i in seq_along(sn)) for (j in seq_along(sn))
    Vs[i, j] <- sV$value[sV$row == sn[i] & sV$col == sn[j]]
  scale <- sqrt(outer(diag(Vr), diag(Vr)))
  data.frame(case = label, matched = length(sn), unmatched = sum(!keep),
    max_b = max(abs(fb - sb$value[keep])),
    max_se_rel = max(abs(sqrt(diag(Vs))/sqrt(diag(Vr)) - 1)),
    max_V_rel = max(abs(Vs - Vr)/scale))
}

ref <- function(case) {
  r <- read.csv("v05_v07_r_reference.csv", stringsAsFactors = FALSE)
  r <- r[r$case == case, ]
  b <- setNames(r$value[r$kind == "b"], r$row[r$kind == "b"])
  V <- matrix(NA, length(b), length(b), dimnames = list(names(b), names(b)))
  v <- r[r$kind == "V", ]
  V[cbind(v$row, v$col)] <- v$value
  list(b = b, V = V)
}

stata <- rbind(parse_ev("test_suest_r_stata_v05.log"), parse_ev("test_suest_r_stata_v06.log"),
  parse_ev("test_suest_r_stata_v07.log"))
out <- list()
run <- function(case, rcase, ..., label = case) {
  r <- ref(rcase)
  out[[length(out) + 1]] <<- compare_case(stata, case, r$b, r$V, ..., label = label)
}

# Round 5, multiple imputation: the R pooled system from the fixture.
mi <- as.data.frame(read_dta("mi_fixture_v01.dta"))
fits <- lapply(1:5, function(m) {
  d <- mi[mi$imp_m == m, ]
  suest(glm(y ~ x, binomial, d), glm(y ~ x + z, binomial, d), model_names = c("mi1", "mi2"))
})
pooled <- suest_mi(fits)
out[[1]] <- compare_case(stata, "mi_suest2", coef(pooled), vcov(pooled), label = "v05 mi_suest2")

run("gaussian_suest2", "lmer", label = "v06 gaussian pair (lmer) mixed")
run("mixed_logit_suest2", "lmer_logit", label = "v06 lmer+logit vs mixed+logit")
run("regress_logit_poisson", "regress_logit_poisson", label = "v07 svy regress+logit+poisson")
r <- ref("logit_sub_poisson")
out[[length(out) + 1]] <- compare_case(stata, "logit_sub_native", r$b, r$V,
  outcomes = "yb", label = "v07 svy logit domain, native block")
run("gologit_full_partial", "gologit_full_partial", label = "v07 gologit2 full+partial")
run("gologit_logit", "gologit_logit", label = "v07 gologit2 partial+logit")
run("heckman_regress", "heckman_regress", heckman = TRUE, label = "v07 heckman+regress")

res <- do.call(rbind, out)
num <- c("max_b", "max_se_rel", "max_V_rel")
res[num] <- lapply(res[num], signif, 3)
print(res, right = FALSE, row.names = FALSE)
write.csv(res, "comparison_v05_v07.csv", row.names = FALSE)
