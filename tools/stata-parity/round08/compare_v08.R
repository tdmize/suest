# Compare test_suest_r_stata_v08.log (EV| lines) with v08_r_reference_v01.csv.
# From compare_v08.R (October 5, 2026), without the Gamma and Laplace
# diagnostics. Writes comparison_v08.csv.
lines <- readLines("test_suest_r_stata_v08.log")
x <- grep("^EV\\|", lines, value = TRUE)
p <- strsplit(x, "|", fixed = TRUE)
stata <- data.frame(case = sapply(p, `[`, 2), kind = sapply(p, `[`, 3), row = sapply(p, `[`, 4),
  col = sapply(p, `[`, 5), value = suppressWarnings(as.numeric(sapply(p, `[`, 6))),
  stringsAsFactors = FALSE)

tf <- list(id = c(function(x) x, function(x) 1),
  exp2 = c(function(x) exp(2*x), function(x) 2*exp(2*x)))

# Stata name -> list(R name, transform)
map_name <- function(name, rnames) {
  name <- sub("^/", "", name)
  eq <- if (grepl(":", name)) sub(":.*$", "", name) else ""
  var <- sub("^[^:]*:", "", name)
  hit <- function(model, candidates, f = "id") {
    r <- paste0(model, "::", candidates)
    r <- r[r %in% rnames]
    if (length(r)) list(r = r[1], tf = f) else NULL
  }
  model <- sub("_.*$", "", eq); suffix <- sub("^[^_]*_", "", eq)
  base <- if (var == "_cons") "(Intercept)" else var
  if (var == "var(_cons[id])") return(hit(model, "log_sigma", "exp2"))
  if (grepl("^cut[0-9]+$", var)) {
    k <- as.integer(sub("cut", "", var))
    cuts <- grep("\\|", sub("^[^:]*::", "", rnames[startsWith(rnames, paste0(model, "::"))]), value = TRUE)
    return(hit(model, cuts[k]))
  }
  if (suffix == "lnsig_e") return(hit(model, "log_sigma_e"))
  if (suffix == "lns1_1_1") return(hit(model, c("log_sd_slope", "log_sigma")))
  if (suffix == "lns1_1_2") return(hit(model, "log_sd_intercept"))
  if (suffix == "atr1_1_1_2") return(hit(model, "atanh_rho"))
  hit(model, base)
}

reference <- function(case) {
  r <- read.csv("v08_r_reference_v01.csv", stringsAsFactors = FALSE)
  r <- r[r$case == case, ]
  b <- setNames(r$value[r$kind == "b"], r$row[r$kind == "b"])
  V <- matrix(NA, length(b), length(b), dimnames = list(names(b), names(b)))
  v <- r[r$kind == "V", ]
  V[cbind(v$row, v$col)] <- v$value
  list(b = b, V = V)
}

compare <- function(scase, rcase = scase, label = rcase) {
  s <- stata[stata$case == scase, ]
  if (!nrow(s)) stop(label, ": no Stata output")
  r <- reference(rcase)
  sb <- s[s$kind == "b", ]
  sV <- s[s$kind == "V", ]
  dv <- setNames(sV$value[sV$row == sV$col], sV$row[sV$row == sV$col])
  maps <- lapply(sb$row, map_name, rnames = names(r$b))
  keep <- !vapply(maps, is.null, NA) & dv[sb$row] > 0
  maps <- maps[keep]; sn <- sb$row[keep]
  rn <- vapply(maps, `[[`, "", "r"); f <- vapply(maps, `[[`, "", "tf")
  fb <- mapply(function(v, t) tf[[t]][[1]](v), r$b[rn], f)
  J <- mapply(function(v, t) tf[[t]][[2]](v), r$b[rn], f)
  Vr <- outer(J, J)*r$V[rn, rn]
  Vs <- outer(sn, sn, Vectorize(function(a, b) sV$value[sV$row == a & sV$col == b]))
  data.frame(case = label, matched = length(sn), max_b = max(abs(fb - sb$value[keep])),
    max_se_rel = max(abs(sqrt(diag(Vs))/sqrt(diag(Vr)) - 1)),
    max_V_rel = max(abs(Vs - Vr)/sqrt(outer(diag(Vr), diag(Vr)))))
}

out <- rbind(
  compare("melogit_aghq25", "logit_aghq25"),
  compare("meprobit_aghq25", "probit_aghq25"),
  compare("mecloglog_aghq25", "cloglog_aghq25"),
  compare("probit_ordinary_aghq25"),
  compare("logit_melogit_aghq25"),
  compare("regress_logit_melogit_aghq25"),
  compare("poisson_aghq25"),
  compare("ologit_oprobit_aghq25"),
  compare("ologit_logit_aghq25"),
  compare("mixed_slopes", "lmer_slopes", "mixed slopes vs lmer"),
  compare("mixed_slopes", "glmmtmb_slopes", "mixed slopes vs glmmTMB"))
num <- c("max_b", "max_se_rel", "max_V_rel")
out[num] <- lapply(out[num], signif, 3)
print(out, row.names = FALSE)
write.csv(out, "comparison_v08.csv", row.names = FALSE)
