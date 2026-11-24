# R side of the R-vs-Stata parity battery for model families that have no
# returned Stata comparison yet, plus ordinary different-sample systems for the
# October 2 combined-sample N/(N-1) correction.
# Writes parity_r_reference_v01.csv (R joint coefficients and covariance, long
# form). Replay copy: the data come from parity_fixture_v01.dta, the file
# Stata used.
suppressPackageStartupMessages({
  library(suest); library(haven); library(survival)
})
set.seed(20261003)
N <- 800
d <- data.frame(id = seq_len(N), x1 = rnorm(N), x2 = rbinom(N, 1, .5), x3 = rnorm(N))
lp <- function(a, b1, b2, b3) a + b1 * d$x1 + b2 * d$x2 + b3 * d$x3
d$yb <- rbinom(N, 1, 1 - exp(-exp(lp(-.5, .6, .4, .3))))
d$yh <- as.integer(lp(-.2, .7, .5, .4) + rnorm(N) * exp(.4 * d$x3) > 0)
d$yf <- rbinom(N, 10, plogis(lp(-.3, .5, .4, .3))) / 10
mu <- plogis(lp(-.2, .5, .3, .3)); d$ybeta <- rbeta(N, mu * 20, (1 - mu) * 20)
lam <- exp(lp(.5, .4, .3, .2)); pz <- plogis(-1 + .8 * d$x2)
d$yz <- ifelse(runif(N) < pz, 0L, rnbinom(N, size = 2, mu = lam))
d$ygam <- rgamma(N, shape = 2, rate = 2 / exp(lp(1, .3, .2, .2)))
d$ytob <- pmax(0, lp(.5, 1, .5, .5) + rnorm(N))
ys <- lp(2, 1, .5, .5) + rnorm(N)
d$lo <- floor(ys); d$hi <- d$lo + 1
d$hi[ys < 0] <- 0; d$lo[ys < 0] <- NA
d$lo[ys >= 5] <- 5; d$hi[ys >= 5] <- NA
tt <- exp(lp(1, .4, .3, .2) + .7 * log(rexp(N)))
cc <- rexp(N, 1 / (3 * median(tt)))
d$time <- pmin(tt, cc); d$event <- as.integer(tt <= cc)
ytr <- lp(1, .8, .5, .5) + rnorm(N); d$ytr <- ifelse(ytr > 0, ytr, NA)
d$ylin <- lp(1, .5, .5, .5) + rnorm(N)
d$ylogit <- rbinom(N, 1, plogis(lp(-.2, .6, .4, .3)))
d <- use_fixture(d, "parity_fixture_v01.dta")

out <- list()
add <- function(case, fit) {
  b <- coef(fit); V <- vcov(fit)
  out[[case]] <<- rbind(
    data.frame(case = case, kind = "b", row = names(b), col = "", value = unname(b)),
    data.frame(case = case, kind = "V", row = rep(rownames(V), ncol(V)),
               col = rep(colnames(V), each = nrow(V)), value = as.vector(V)),
    data.frame(case = case, kind = "N", row = "union", col = "", value = fit$nobs_union))
}
two <- function(case, fA, fB, ...) add(case, suest(fA, fB, model_names = c("A", "B"), ...))

# Same-sample pairs: A = x1 + x2, B = x1 + x2 + x3
two("cloglog",
    glm(yb ~ x1 + x2, binomial("cloglog"), d), glm(yb ~ x1 + x2 + x3, binomial("cloglog"), d))
two("gamma_log",
    glm(ygam ~ x1 + x2, Gamma("log"), d), glm(ygam ~ x1 + x2 + x3, Gamma("log"), d))
two("fracreg_logit",
    glm(yf ~ x1 + x2, quasibinomial("logit"), d), glm(yf ~ x1 + x2 + x3, quasibinomial("logit"), d))
two("fracreg_probit",
    glm(yf ~ x1 + x2, quasibinomial("probit"), d), glm(yf ~ x1 + x2 + x3, quasibinomial("probit"), d))
two("zip",
    pscl::zeroinfl(yz ~ x1 + x2 | x2, d, dist = "poisson"),
    pscl::zeroinfl(yz ~ x1 + x2 + x3 | x2, d, dist = "poisson"))
two("zinb",
    pscl::zeroinfl(yz ~ x1 + x2 | x2, d, dist = "negbin"),
    pscl::zeroinfl(yz ~ x1 + x2 + x3 | x2, d, dist = "negbin"))
two("tobit",
    censReg::censReg(ytob ~ x1 + x2, left = 0, data = d),
    censReg::censReg(ytob ~ x1 + x2 + x3, left = 0, data = d))
two("intreg",
    survreg(Surv(lo, hi, type = "interval2") ~ x1 + x2, d, dist = "gaussian"),
    survreg(Surv(lo, hi, type = "interval2") ~ x1 + x2 + x3, d, dist = "gaussian"))
for (dist in c("weibull", "exponential", "lognormal", "loglogistic"))
  two(paste0("streg_", dist),
      survreg(Surv(time, event) ~ x1 + x2, d, dist = dist),
      survreg(Surv(time, event) ~ x1 + x2 + x3, d, dist = dist))
dtr <- d[!is.na(d$ytr), ]
two("truncreg",
    truncreg::truncreg(ytr ~ x1 + x2, dtr, point = 0, direction = "left"),
    truncreg::truncreg(ytr ~ x1 + x2 + x3, dtr, point = 0, direction = "left"))
two("betareg",
    betareg::betareg(ybeta ~ x1 + x2, d), betareg::betareg(ybeta ~ x1 + x2 + x3, d))
two("hetprobit",
    Rchoice::hetprob(yh ~ x1 + x2 | x3, d, link = "probit"),
    Rchoice::hetprob(yh ~ x1 + x2 + x3 | x3, d, link = "probit"))

# Different samples, ordinary unclustered systems (combined-sample N/(N-1))
pA <- d$id <= 500; pB <- d$id > 300; jA <- d$id <= 400; jB <- d$id > 400
two("regress_partial", lm(ylin ~ x1 + x2, d, subset = pA), lm(ylin ~ x1 + x2 + x3, d, subset = pB))
two("regress_disjoint", lm(ylin ~ x1 + x2, d, subset = jA), lm(ylin ~ x1 + x2 + x3, d, subset = jB))
two("logit_partial", glm(ylogit ~ x1 + x2, binomial, d, subset = pA), glm(ylogit ~ x1 + x2 + x3, binomial, d, subset = pB))
two("logit_disjoint", glm(ylogit ~ x1 + x2, binomial, d, subset = jA), glm(ylogit ~ x1 + x2 + x3, binomial, d, subset = jB))
two("cloglog_partial", glm(yb ~ x1 + x2, binomial("cloglog"), d, subset = pA), glm(yb ~ x1 + x2 + x3, binomial("cloglog"), d, subset = pB))
two("tobit_partial",
    censReg::censReg(ytob ~ x1 + x2, left = 0, data = d[pA, ]),
    censReg::censReg(ytob ~ x1 + x2 + x3, left = 0, data = d[pB, ]), observation_id = "id")

res <- do.call(rbind, out)
write.csv(res, "parity_r_reference_v01.csv", row.names = FALSE)
