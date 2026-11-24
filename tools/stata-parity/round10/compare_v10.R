# Compare test_suest_r_stata_v10.log (EV| lines) with rs_r_reference_v01.csv.
# Each Stata parameter is written as a function of the R parameters (Stata
# reports variances and covariances; R uses log SDs and Fisher z
# correlations), and covariances use the Jacobian of those functions.
# From compare_v10.R (November 7, 2026), keeping the linear multilevel and
# survey cases. Writes comparison_v10.csv.
lines <- readLines("test_suest_r_stata_v10.log")
x <- grep("^EV\\|", lines, value = TRUE)
p <- strsplit(x, "|", fixed = TRUE)
stata <- data.frame(case = sapply(p, `[`, 2), kind = sapply(p, `[`, 3), row = sapply(p, `[`, 4),
  col = sapply(p, `[`, 5), value = suppressWarnings(as.numeric(sapply(p, `[`, 6))),
  stringsAsFactors = FALSE)

ref <- read.csv("rs_r_reference_v01.csv", stringsAsFactors = FALSE)
reference <- function(case) {
  r <- ref[ref$case == case, ]
  b <- setNames(r$value[r$kind == "b"], r$row[r$kind == "b"])
  V <- matrix(NA, length(b), length(b), dimnames = list(names(b), names(b)))
  v <- r[r$kind == "V", ]
  V[cbind(v$row, v$col)] <- v$value
  list(b = b, V = V)
}

column <- function(x) if (x == "_cons") "intercept" else x
group <- function(x) gsub(">", ".", x, fixed = TRUE)   # school>class -> school.class
log_sd <- function(m, g, c) function(b) b[[paste0(m, "::log_sd_", g, "_", c)]]
rho <- function(m, g, c1, c2) {
  n1 <- paste0(m, "::atanh_rho_", g, "_", c1, "_", c2)
  n2 <- paste0(m, "::atanh_rho_", g, "_", c2, "_", c1)
  function(b) tanh(if (n1 %in% names(b)) b[[n1]] else b[[n2]])
}
variance <- function(m, g, c) { s <- log_sd(m, g, c); function(b) exp(2*s(b)) }
covariance <- function(m, g, c1, c2) {
  s1 <- log_sd(m, g, c1); s2 <- log_sd(m, g, c2); r <- rho(m, g, c1, c2)
  function(b) r(b)*exp(s1(b) + s2(b))
}
plain <- function(r) function(b) b[[r]]

# Stata name -> function of the R coefficient vector (NULL: not compared).
#   letter: model letter for one-model results; levels: mixed's random-effect
#   columns at each level, in Stata's order.
map_name <- function(name, rnames, letter = NULL, levels = NULL) {
  name <- sub("^/", "", name)
  eq <- if (grepl(":", name)) sub(":.*$", "", name) else ""
  var <- sub("^[^:]*:", "", name)
  if (!is.null(letter)) { model <- letter; suffix <- eq }
  else { model <- sub("_.*$", "", eq); suffix <- if (grepl("_", eq)) sub("^[^_]*_", "", eq) else "" }
  lv <- levels[[model]]
  m <- regmatches(var, regexec("^var\\((.+)\\[(.+)\\]\\)$", var))[[1]]
  if (length(m)) return(variance(model, group(m[3]), column(m[2])))
  m <- regmatches(var, regexec("^cov\\((.+)\\[(.+)\\],(.+)\\[(.+)\\]\\)$", var))[[1]]
  if (length(m)) return(covariance(model, group(m[3]), column(m[2]), column(m[4])))
  if (grepl("^lns[0-9]+_1_[0-9]+$", suffix)) {
    k <- as.integer(strsplit(sub("^lns", "", suffix), "_")[[1]])
    return(log_sd(model, names(lv)[k[1]], lv[[k[1]]][k[3]]))
  }
  if (grepl("^atr[0-9]+_1_[0-9]+_[0-9]+$", suffix)) {
    k <- as.integer(strsplit(sub("^atr", "", suffix), "_")[[1]])
    cols <- lv[[k[1]]]
    r <- rho(model, names(lv)[k[1]], cols[k[3]], cols[k[4]])
    return(function(b) atanh(r(b)))
  }
  if (suffix == "lnsig_e") return(plain(paste0(model, "::log_sigma_e")))
  if (startsWith(var, "o.") || grepl("^[0-9]+b?\\.", var)) return(NULL)
  r <- paste0(model, "::", if (var == "_cons") "(Intercept)" else var)
  if (!r %in% rnames) stop("no R parameter for Stata ", name)
  plain(r)
}

compare <- function(scase, rcase, label = scase, ...) {
  s <- stata[stata$case == scase, ]
  if (!nrow(s)) stop(label, ": no Stata output")
  r <- reference(rcase)
  sb <- s[s$kind == "b", ]; sV <- s[s$kind == "V", ]
  dv <- setNames(sV$value[sV$row == sV$col], sV$row[sV$row == sV$col])
  keep <- !is.na(dv[sb$row]) & dv[sb$row] > 0
  sn <- sb$row[keep]
  maps <- lapply(sn, map_name, rnames = names(r$b), ...)
  use <- !vapply(maps, is.null, NA)
  maps <- maps[use]; sn <- sn[use]
  f <- function(b) vapply(maps, function(g) g(b), 0)
  fb <- f(r$b)
  J <- numDeriv::jacobian(f, r$b)
  used <- colSums(abs(J)) > 0
  Vr <- J[, used, drop = FALSE] %*% r$V[used, used] %*% t(J[, used, drop = FALSE])
  Vs <- matrix(sV$value[match(paste(rep(sn, each = length(sn)), rep(sn, length(sn))),
    paste(sV$row, sV$col))], length(sn), byrow = TRUE)
  free <- diag(Vr) > 0
  data.frame(case = label, matched = length(sn),
    max_b = max(abs(fb - sb$value[keep][use])),
    max_se_rel = max(abs(sqrt(diag(Vs)[free])/sqrt(diag(Vr)[free]) - 1)),
    max_V_rel = max((abs(Vs - Vr)/sqrt(outer(diag(Vr), diag(Vr))))[free, free]))
}

mixed_levels <- list(
  A = list(school = c("x1", "x2", "intercept")),
  B = list(school = c("x1", "intercept"), school.class = "intercept"))
out <- rbind(
  compare("mixed_pair", "mixed_pair", "mixed pair, suest2", levels = mixed_levels),
  compare("mixed_pair_cluster_A", "mixed_pair", "mixed A alone, clustered", letter = "A",
    levels = mixed_levels),
  compare("mixed_pair_cluster_B", "mixed_pair", "mixed B alone, clustered", letter = "B",
    levels = mixed_levels),
  compare("mixed_logit", "mixed_logit", "mixed + logit, suest2",
    levels = list(A = list(school = c("x1", "intercept"), school.class = "intercept"))),
  compare("mixed_logit_cluster_B", "mixed_logit", "logit alone, clustered", letter = "B"),
  compare("svy_cloglog_pair", "svy_cloglog_pair", "svy cloglog pair, suest2"),
  compare("svy_cloglog_pair_cluster_A", "svy_cloglog_pair", "svy cloglog A alone", letter = "A"),
  compare("svy_cloglog_pair_cluster_B", "svy_cloglog_pair", "svy cloglog B alone", letter = "B"))
num <- c("max_b", "max_se_rel", "max_V_rel")
out[num] <- lapply(out[num], signif, 3)
print(out, row.names = FALSE, right = FALSE)
write.csv(out, "comparison_v10.csv", row.names = FALSE)
