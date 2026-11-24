# R side of Stata round 18 (Batch 4): Gompertz and generalized gamma survival
# models (flexsurv::flexsurvreg, Stata's streg with distribution(gompertz) and
# distribution(ggamma)), with ancillary covariates and delayed entry, and
# bivariate probit with different predictors in each equation (GJRM::gjrm,
# Stata's biprobit).
# Writes batch4_r_reference_v01.csv. Replay copy: the data come from
# batch4_fixture_v01.dta, the file Stata used.
# The flexsurv fits are polished: flexsurv's optimizer (BFGS with numerical
# gradients for the generalized gamma) stops short of the maximum (by up to
# 1.5e-4 in a coefficient here with flexsurv's defaults, and 4e-5 with
# reltol = 1e-14), so each fit takes Newton steps with suest's scores and information
# and is refit at the result (inits, maxit = 0). This keeps optimizer
# differences out of the comparison. GJRM fits are used as they are: GJRM
# approximates the bivariate normal (relative error about 1e-8), so its
# estimates can differ from Stata's by about 1e-7 (1e-7 and 5e-8 here).
# Usage: Rscript make_v18.R   (SUEST_DEV=<path> loads a development copy)
suppressPackageStartupMessages({
  if (nzchar(Sys.getenv("SUEST_DEV"))) pkgload::load_all(Sys.getenv("SUEST_DEV"), quiet = TRUE) else
    library(suest)
  library(haven); library(flexsurv); library(GJRM); library(marginaleffects)
})
set.seed(20261108)
n <- 800L
d <- data.frame(rid = seq_len(n))
d$x <- round(rnorm(n), 6); d$z <- round(rnorm(n), 6); d$w <- round(rnorm(n), 6)
d$b <- rbinom(n, 1, 0.4)
# Gompertz times (shape depends on z), censored at a uniform time
Tg <- rgompertz(n, shape = 0.2 + 0.08*d$z, rate = exp(-2 + 0.5*d$x - 0.3*d$z))
C <- runif(n, 1, 12)
d$tm <- round(pmin(Tg, C), 6); d$ev <- as.integer(Tg <= C)
# Generalized gamma times (sigma depends on z, Q on b; Q stays away from 0)
Tq <- rgengamma(n, mu = 1 + 0.4*d$x - 0.2*d$z, sigma = exp(-0.3 + 0.2*d$z), Q = 0.5 + 0.4*d$b)
C2 <- runif(n, 0.5, 15)
d$tg <- round(pmin(Tq, C2), 6); d$eg <- as.integer(Tq <= C2)
# Delayed entry for 40% of the rows (always before the exit time)
d$entry <- ifelse(runif(n) < 0.4, round(runif(n)*0.6*d$tm, 6), 0)
# Two binary outcomes with correlated errors (0.4), and a third probit outcome
e <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, 0.4, 0.4, 1), 2))
d$y1 <- as.integer(0.2 + 0.5*d$x - 0.4*d$z + e[, 1] > 0)
d$y2 <- as.integer(-0.1 + 0.3*d$x + 0.6*d$w + e[, 2] > 0)
d$yp <- as.integer(0.1 + 0.4*d$x + 0.3*d$w + 0.5*e[, 1] + rnorm(n, sd = sqrt(0.75)) > 0)
stopifnot(all(d$entry < d$tm), all(d$tm > 0), all(d$tg > 0))
d <- use_fixture(d, "batch4_fixture_v01.dta")

polish <- function(model) {
  for (i in 1:8) {
    comp <- suest:::.suest_flexsurv_components(model)
    s <- colSums(comp$score)
    if (max(abs(s)) < 1e-8) break
    theta <- comp$parameters + drop(comp$bread %*% s)/nrow(comp$score)
    inits <- unname(theta)
    for (j in seq_along(model$dlist$pars)) inits[j] <- model$dlist$inv.transforms[[j]](inits[j])
    call <- model$call
    call$inits <- inits
    call$control <- list(maxit = 0)
    model <- eval(call, globalenv())
  }
  cat(sprintf("  polished: %d Newton steps, largest score sum %.2g\n", i - 1, max(abs(s))))
  model
}
ctl <- list(reltol = 1e-14, maxit = 5000)
fs <- function(f, dist, anc = NULL) {
  call <- bquote(flexsurvreg(.(f), data = d, dist = .(dist), control = ctl))
  if (!is.null(anc)) call$anc <- anc
  default <- call; default$control <- NULL
  m0 <- eval(default, globalenv())
  m <- polish(eval(call, globalenv()))
  cat(sprintf("  flexsurv's default convergence: largest coefficient gap %.2g\n",
    max(abs(coef(m0) - coef(m)))))
  m
}
cat("Fits\n")
gA <- fs(Surv(tm, ev) ~ x + z, "gompertz")
gB <- fs(Surv(tm, ev) ~ x + z, "gompertz", quote(list(shape = ~ z)))
qA <- fs(Surv(tg, eg) ~ x + z, "gengamma")
qB <- fs(Surv(tg, eg) ~ x + z, "gengamma", quote(list(sigma = ~ z, Q = ~ b)))
eA <- fs(Surv(entry, tm, ev) ~ x + z, "gompertz")
eB <- fs(Surv(entry, tm, ev) ~ x + z, "gompertz", quote(list(shape = ~ z)))
bA <- gjrm(list(y1 ~ x + z, y2 ~ x + w), data = d, margins = c("probit", "probit"), model = "B")
bB <- gjrm(list(y1 ~ x, y2 ~ x + w + z), data = d, margins = c("probit", "probit"), model = "B")
pC <- glm(yp ~ x + w, family = binomial(link = "probit"), data = d,
  control = glm.control(epsilon = 1e-14, maxit = 100))

out <- list()
long <- function(case, b, V) rbind(
  data.frame(case = case, kind = "b", row = names(b), col = "", value = unname(b)),
  data.frame(case = case, kind = "V", row = rep(rownames(V), ncol(V)),
    col = rep(colnames(V), each = nrow(V)), value = as.vector(V)))
info <- function(case, ...) {
  x <- c(...)
  out[[case]] <<- data.frame(case = case, kind = "info", row = names(x), col = "", value = unname(x))
}
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
  cat("\n==", case, "\n"); print(signif(cbind(b = b, se = sqrt(diag(V))), 10))
  invisible(fit)
}
# Average prediction of one model in a system (Stata: that model's native
# vce(robust) fit, then margins): median survival time, or Pr(y1 = 1, y2 = 1).
avg <- function(case, fit, model) {
  p <- as.data.frame(avg_predictions(fit, newdata = d))
  p <- p[p$group == model, ]
  stopifnot(nrow(p) == 1)
  b <- c("_cons" = p$estimate); V <- matrix(p$std.error^2, 1, 1, dimnames = list("_cons", "_cons"))
  out[[case]] <<- long(case, b, V)
  cat("\n==", case, ": ", sprintf("%.10g (%.10g)", p$estimate, p$std.error), "\n", sep = "")
}

fit_g <- add("gomp_pair", gA, gB)
fit_q <- add("ggamma_pair", qA, qB)
add("entry_pair", eA, eB)
fit_b <- add("biprobit_pair", bA, bB)
add("three_mixed", gA, bA, pC)
avg("gomp_A_margins", fit_g, "A")
avg("ggamma_B_margins", fit_q, "B")
avg("biprobit_A_margins", fit_b, "A")
for (k in c("gA", "gB", "qA", "qB", "eA", "eB")) {
  m <- get(k); info(paste0(k, "_info"), ll = m$loglik, N = m$N, events = m$events)
}
for (k in c("bA", "bB")) info(paste0(k, "_info"), ll = as.numeric(logLik(get(k))), N = n)
info("pC_info", ll = as.numeric(logLik(pC)), N = n)
ref <- do.call(rbind, out)
write.csv(ref, "batch4_r_reference_v01.csv", row.names = FALSE)
