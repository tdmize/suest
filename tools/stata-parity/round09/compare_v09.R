# Compare the returned Stata log v09 (lines starting EV|) with the R reference
# survey_models_r_reference_v01.csv. Each Stata parameter is mapped to an R
# parameter and a transformation (Stata value = f(R value)); covariances use
# the Jacobian. Native cases (svy: e(V) of one model) are compared with that
# model's block of the R system.
parse_ev <- function(file) {
  lines <- grep("^EV\\|", readLines(file), value = TRUE)
  parts <- strsplit(lines, "|", fixed = TRUE)
  data.frame(case = sapply(parts, `[`, 2), kind = sapply(parts, `[`, 3),
    row = sapply(parts, `[`, 4), col = sapply(parts, `[`, 5),
    value = as.numeric(sapply(parts, `[`, 6)), stringsAsFactors = FALSE)
}
transforms <- list(
  id = list(f = function(x) x, d = function(x) 1),
  neg = list(f = function(x) -x, d = function(x) -1),
  exp2 = list(f = function(x) exp(2 * x), d = function(x) 2 * exp(2 * x)))

# Model type of each model letter in each system.
types <- list(
  ologit_pair = c(A = "ordered", B = "ordered"), oprobit_pair = c(A = "ordered", B = "ordered"),
  mlogit_pair = c(A = "mlogit", B = "mlogit"), nbreg_pair = c(A = "nbreg", B = "nbreg"),
  gologit_pair = c(A = "gologit", B = "gologit_plz"), tobit_pair = c(A = "tobit", B = "tobit"),
  intreg_pair = c(A = "intreg", B = "intreg"), weibull_pair = c(A = "weibull", B = "weibull"),
  ivreg_pair = c(A = "linear", B = "linear"),
  mixed_six = c(A = "linear", B = "linear", C = "ordered", D = "mlogit", E = "nbreg", F = "tobit"),
  domain_pair = c(A = "ordered", B = "mlogit"))

# Native names ("o:x", "/cut1", "endo") to the joint form ("A_o:x", "/A:cut1", "A:endo").
joint_name <- function(name, model) {
  if (startsWith(name, "/")) return(paste0("/", model, ":", sub("^/", "", name)))
  if (grepl(":", name)) return(paste0(model, "_", name))
  paste0(model, ":", name)
}

map_one <- function(name, system) {
  if (startsWith(name, "/")) {
    model <- sub("^/([A-Z]):.*$", "\\1", name); eq <- ""; var <- sub("^/[A-Z]:", "", name)
  } else {
    model <- substr(name, 1, 1)
    eq <- sub("^[A-Z]_?", "", sub(":.*$", "", name)); var <- sub("^[^:]*:", "", name)
  }
  if (startsWith(var, "o.")) return(NULL)
  type <- types[[system]][[model]]
  term <- if (var == "_cons") "(Intercept)" else var
  out <- function(r, tf = "id") list(r = paste0(model, "::", r), tf = tf)
  switch(type,
    ordered = if (grepl("^cut", var)) {
      k <- as.integer(sub("cut", "", var)); out(paste0(k, "|", k + 1))
    } else out(term),
    mlogit = out(paste0(term, ":", as.integer(eq) - 1L)),
    nbreg = if (var == "lnalpha") out("(Intercept):2", "neg") else
      out(if (var == "_cons") "(Intercept):1" else term),
    gologit = out(paste0(term, ":", eq)),
    gologit_plz = if (var == "z") out("z") else out(paste0(term, ":", eq)),
    tobit = if (var == "var(e.yc)") out("Log(scale)", "exp2") else out(term),
    intreg = if (eq == "lnsigma") out("Log(scale)") else out(term),
    weibull = if (var == "ln_p") out("Log(scale)", "neg") else out(term),
    linear = out(term))
}

compare_case <- function(s, system, rb, rV, label) {
  sb <- s[s$kind == "b", ]
  maps <- lapply(sb$row, map_one, system = system)
  keep <- !vapply(maps, is.null, logical(1))
  maps <- maps[keep]; sn <- sb$row[keep]
  rn <- vapply(maps, `[[`, "", "r"); tf <- vapply(maps, `[[`, "", "tf")
  if (!all(rn %in% rownames(rV))) stop(label, ": missing R names ", paste(setdiff(rn, rownames(rV)), collapse = ", "))
  rval <- rb[rn]
  fb <- mapply(function(v, t) transforms[[t]]$f(v), rval, tf)
  J <- mapply(function(v, t) transforms[[t]]$d(v), rval, tf)
  Vr <- outer(J, J) * rV[rn, rn]
  sV <- s[s$kind == "V", ]
  Vs <- matrix(sV$value[match(paste(rep(sn, each = length(sn)), rep(sn, length(sn))),
    paste(sV$row, sV$col))], length(sn), byrow = TRUE)
  scale <- sqrt(outer(diag(Vr), diag(Vr)))
  data.frame(case = label, params = length(sn),
    max_b = max(abs(fb - sb$value[keep])),
    max_se_rel = max(abs(sqrt(diag(Vs)) / sqrt(diag(Vr)) - 1)),
    max_V_rel = max(abs(Vs - Vr) / scale))
}

stata <- parse_ev("test_suest_r_stata_v09.log")
ref <- read.csv("survey_models_r_reference_v01.csv", stringsAsFactors = FALSE)
out <- list()
for (system in names(types)) {
  r <- ref[ref$case == system, ]
  rb <- setNames(r$value[r$kind == "b"], r$row[r$kind == "b"])
  v <- r[r$kind == "V", ]
  rV <- matrix(NA, length(rb), length(rb), dimnames = list(names(rb), names(rb)))
  rV[cbind(v$row, v$col)] <- v$value
  out[[length(out) + 1]] <- compare_case(stata[stata$case == system, ], system, rb, rV,
    paste(system, "suest2"))
  for (model in names(types[[system]])) {
    s <- stata[stata$case == paste0(system, "_native_", model), ]
    s$row <- vapply(s$row, joint_name, "", model = model)
    s$col[s$kind == "V"] <- vapply(s$col[s$kind == "V"], joint_name, "", model = model)
    if (types[[system]][[model]] == "linear" && system == "ivreg_pair")
      s$row <- sub("^([A-Z]):", "\\1_mean:", s$row)
    if (types[[system]][[model]] == "linear" && system == "ivreg_pair")
      s$col <- sub("^([A-Z]):", "\\1_mean:", s$col)
    out[[length(out) + 1]] <- compare_case(s, system, rb, rV,
      paste0(system, " native ", model))
  }
}
res <- do.call(rbind, out)
res[c("max_b", "max_se_rel", "max_V_rel")] <- lapply(res[c("max_b", "max_se_rel", "max_V_rel")], signif, 3)
print(res, right = FALSE, row.names = FALSE)
write.csv(res, "comparison_v09.csv", row.names = FALSE)
