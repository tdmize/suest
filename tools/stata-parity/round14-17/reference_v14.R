# R side of round 14 (Batch 3): constant-dispersion negative binomial, GEE with
# ar1 and unstructured working correlations, and generalized ordered models
# with cauchit and log-log links. From make_v14.R (November 8, 2026), without
# the do-file template and the diagnostics, reading batch3_fixture_v01.dta,
# the data Stata used. Writes batch3_r_reference_v02.csv. glmmTMB fits use
# BFGS with a tight tolerance, because glmmTMB's default optimizer stops just
# short of the maximum.
suppressPackageStartupMessages({
  library(suest); library(haven); library(geepack); library(glmmTMB); library(VGAM); library(MASS)
})
d <- as.data.frame(read_dta("batch3_fixture_v01.dta"))
d$ocf <- ordered(d$oc); d$ogf <- ordered(d$og)
short <- d[d$short == 1, ]
gctl <- geese.control(epsilon = 1e-12, maxit = 200)
vctl <- vglm.control(criterion = "coefficients", epsilon = 1e-10, maxit = 500)
tctl <- glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"),
  optCtrl = list(reltol = 1e-15, maxit = 2000))
# Fit in the global environment with the data named in the call, so that
# suest can find each model's data for observation_id.
gee <- function(f, data, corstr, family = gaussian(), waves = FALSE) {
  call <- bquote(geeglm(.(f), family = .(as.name(family$family))(), data = .(substitute(data)),
    id = id, corstr = .(corstr), control = gctl))
  if (waves) call$waves <- quote(time)
  eval(call, globalenv())
}

out <- list()
long <- function(case, b, V) rbind(
  data.frame(case = case, kind = "b", row = names(b), col = "", value = unname(b)),
  data.frame(case = case, kind = "V", row = rep(rownames(V), ncol(V)),
    col = rep(colnames(V), each = nrow(V)), value = as.vector(V)))
add <- function(case, ...) {
  models <- list(...)
  fit <- suest(..., model_names = LETTERS[seq_along(models)], observation_id = "rid")
  b <- coef(fit); V <- vcov(fit)
  out[[case]] <<- long(case, b, V)
  # Each model's own block (Stata's native vce(robust) fits).
  for (k in seq_along(models)) {
    rows <- startsWith(names(b), paste0(LETTERS[k], "::"))
    out[[paste0(case, "_native_", LETTERS[k])]] <<- long(paste0(case, "_native_", LETTERS[k]),
      b[rows], V[rows, rows, drop = FALSE])
  }
  invisible(fit)
}

add("nb1_pair", glmmTMB(cnt1 ~ x + z, d, family = nbinom1(), control = tctl),
  glmmTMB(cnt1 ~ x, d, family = nbinom1(), control = tctl))
add("nb_three", glmmTMB(cnt1 ~ x + z, d, family = nbinom1(), control = tctl),
  glm.nb(cnt1 ~ x + z, d, control = glm.control(epsilon = 1e-12, maxit = 100)),
  glm(cnt1 ~ x + z, poisson, d, control = glm.control(epsilon = 1e-12)))
g1 <- gee(yg ~ x + z, d, "ar1"); g2 <- gee(yg ~ x, d, "ar1")
add("gee_ar1", g1, g2)
g1 <- gee(yg ~ x + z, d, "unstructured"); g2 <- gee(yg ~ x, d, "unstructured")
add("gee_un", g1, g2)
g1 <- gee(yb ~ x + z, d, "ar1", binomial()); g2 <- gee(yc ~ x + z, d, "unstructured", poisson())
add("gee_mixed", g1, g2)
g1 <- gee(yg ~ x + z, short, "unstructured", waves = TRUE)
g2 <- gee(yb ~ x + z, short, "ar1", binomial(), waves = TRUE)
add("gee_short", g1, g2)
add("cauchit_pair",
  vglm(ocf ~ x + z, cumulative(link = "cauchitlink", parallel = FALSE ~ x, reverse = TRUE),
    data = d, control = vctl),
  vglm(ocf ~ x, cumulative(link = "cauchitlink", parallel = FALSE, reverse = TRUE),
    data = d, control = vctl))
# Log-log (gologit2, link(loglog)): VGAM cloglog with reverse = FALSE; every
# coefficient is minus Stata's.
add("loglog_pair",
  vglm(ogf ~ x + z, cumulative(link = "clogloglink", parallel = FALSE ~ x, reverse = FALSE),
    data = d, control = vctl),
  vglm(ogf ~ x, cumulative(link = "clogloglink", parallel = FALSE, reverse = FALSE),
    data = d, control = vctl))
add("three_mixed",
  vglm(ocf ~ x + z, cumulative(link = "cauchitlink", parallel = FALSE ~ x, reverse = TRUE),
    data = d, control = vctl),
  glmmTMB(cnt1 ~ x + z, d, family = nbinom1(), control = tctl),
  glm(yb ~ x + z, binomial, d, control = glm.control(epsilon = 1e-12)))

write.csv(do.call(rbind, out), "batch3_r_reference_v02.csv", row.names = FALSE)
