# R side of test_suest_r_stata_v08.do. Writes v08_r_reference_v01.csv (joint b
# and V in long form, R names). Replay copy: the data come from the fixtures.
#  - Random-effects models with 25-point adaptive quadrature in both programs
#    (R: glmer/clmm nAGQ = 25, mode-curvature; Stata: mvaghermite intpoints(25)).
#    At 25 points both are effectively exact integrals, so the two agree
#    closely even though their adaptive rules differ.
#  - Gaussian random slopes: lmer and glmmTMB vs mixed, covariance(unstructured).
#  (The Gamma and Laplace diagnostics of the original script are left out.)
# Reuses families_fixture_v01.dta, combined_fixture_v01.dta, clmm_fixture_v01.dta.
suppressPackageStartupMessages({
  library(suest); library(haven); library(lme4); library(glmmTMB); library(ordinal)
})
set.seed(20261013)
G <- 100; m <- 8
v8 <- data.frame(id = rep(seq_len(G), each = m), x = rnorm(G * m), z = rnorm(G * m))
u0 <- rnorm(G, sd = 0.8); u1 <- 0.3 * u0 + rnorm(G, sd = 0.4)
v8$ys1 <- 1 + 0.5 * v8$x - 0.3 * v8$z + rep(u0, each = m) + rep(u1, each = m) * v8$x + rnorm(G * m)
v8$ys2 <- 0.5 + 0.4 * v8$x + 0.7 * rep(u0, each = m) + 0.8 * rep(u1, each = m) * v8$x +
  rnorm(G * m, sd = 1.2)
v8$yp1 <- rpois(G * m, exp(0.2 + 0.3 * v8$x - 0.2 * v8$z + 0.5 * rep(u0, each = m)))
v8$yp2 <- rpois(G * m, exp(-0.1 + 0.4 * v8$x + 0.6 * rep(u1, each = m)))
v8 <- use_fixture(v8, "v08_fixture_v01.dta")
v8$id <- factor(v8$id)

fam <- as.data.frame(read_dta("families_fixture_v01.dta")); fam$id <- factor(fam$id)
comb <- as.data.frame(read_dta("combined_fixture_v01.dta")); comb$id <- factor(comb$id)
ord <- as.data.frame(read_dta("clmm_fixture_v01.dta")); ord$id <- factor(ord$id)
ord$yo <- ordered(ord$y)

out <- list()
add_bv <- function(case, b, V) {
  out[[case]] <<- rbind(
    data.frame(case = case, kind = "b", row = names(b), col = "", value = unname(b)),
    data.frame(case = case, kind = "V", row = rep(rownames(V), ncol(V)),
               col = rep(colnames(V), each = nrow(V)), value = as.vector(V)))
  cat("\n==", case, "\n"); print(signif(cbind(b = b, se = sqrt(diag(V))), 10))
}
add <- function(case, fit) add_bv(case, coef(fit), vcov(fit))

ctl <- glmerControl(tolPwrss = 1e-12, optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
q <- function(f, family, data) glmer(f, family = family, data = data, nAGQ = 25, control = ctl)
AB <- c("A", "B")

# Part 1. Binary random-intercept pairs, 25-point adaptive quadrature
add("logit_aghq25", suest(q(yb1 ~ x + z + (1 | id), binomial, fam),
                          q(yb2 ~ x + (1 | id), binomial, fam), model_names = AB))
add("probit_aghq25", suest(q(yb1 ~ x + z + (1 | id), binomial("probit"), fam),
                           q(yb2 ~ x + (1 | id), binomial("probit"), fam), model_names = AB))
add("cloglog_aghq25", suest(q(yb1 ~ x + z + (1 | id), binomial("cloglog"), fam),
                            q(yb2 ~ x + (1 | id), binomial("cloglog"), fam), model_names = AB))
add("probit_ordinary_aghq25", suest(glm(yb1 ~ x + z, binomial("probit"), fam),
                                    q(yb2 ~ x + (1 | id), binomial("probit"), fam),
                                    model_names = AB))
# Ordinary + random-intercept logit (combined fixture)
add("logit_melogit_aghq25", suest(glm(y1 ~ x + z, binomial, comb),
                                  q(y2 ~ x + (1 | id), binomial, comb), model_names = AB))
add("regress_logit_melogit_aghq25", suest(lm(yc ~ x + z, comb), glm(y1 ~ x + z, binomial, comb),
                                          q(y2 ~ x + (1 | id), binomial, comb),
                                          model_names = c("A", "B", "C")))
# Poisson pair
add("poisson_aghq25", suest(q(yp1 ~ x + z + (1 | id), poisson, v8),
                            q(yp2 ~ x + (1 | id), poisson, v8), model_names = AB))

# Part 2. Random-intercept ordered models, 25-point adaptive quadrature
cctl <- clmm.control(gradTol = 1e-8)
ol <- clmm(yo ~ x + z + (1 | id), data = ord, nAGQ = 25, control = cctl)
op <- clmm(yo ~ x + (1 | id), data = ord, nAGQ = 25, link = "probit", control = cctl)
add("ologit_oprobit_aghq25", suest(ol, op, model_names = AB))
add("ologit_logit_aghq25", suest(ol, glm(yb ~ x, binomial, ord), model_names = AB))

# Part 3. Gaussian random slopes (exact in both programs)
lm1 <- lmer(ys1 ~ x + z + (1 + x | id), data = v8, REML = FALSE)
lm2 <- lmer(ys2 ~ x + (1 + x | id), data = v8, REML = FALSE)
add("lmer_slopes", suest(lm1, lm2, model_names = AB))
add("glmmtmb_slopes", suest(glmmTMB(ys1 ~ x + z + (1 + x | id), data = v8),
                            glmmTMB(ys2 ~ x + (1 + x | id), data = v8), model_names = AB))

write.csv(do.call(rbind, out), "v08_r_reference_v01.csv", row.names = FALSE)
