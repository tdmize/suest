# R side of the two-part round (November 23, 2026; suest(twopart = TRUE)
# against Stata's suest2 two-part systems, test_twopart_stata_v03.do) and the
# comparison with its log (ME| lines). From compare_twopart_v02.R, without
# the mecompare printouts. Writes comparison_twopart.csv.
suppressPackageStartupMessages({
  library(suest); library(marginaleffects); library(haven); library(survey)
})
meps <- zap_labels(read_dta(data_file("meps.dta")))
for (v in c("female", "ins_unins", "anylim")) meps[[v]] <- factor(meps[[v]])
meps$any <- as.numeric(meps$exp_tot > 0)
meps$lnexp <- ifelse(meps$any == 1, log(pmax(meps$exp_tot, 1)), NA)
meps$w <- meps$wtdper / mean(meps$wtdper)
tight <- glm.control(epsilon = 1e-10, maxit = 100)
f7 <- ~ age + female + lninc + pcs12 + mcs12 + ins_unins + anylim
pos <- meps$any == 1

any <- glm(update(f7, any ~ .), family = binomial, data = meps, control = tight)
amount <- glm(update(f7, exp_tot ~ .), family = quasipoisson, data = meps, subset = pos, control = tight)
gam <- glm(update(f7, exp_tot ~ .), family = Gamma(link = "log"), data = meps, subset = pos, control = tight)
nb <- MASS::glm.nb(update(f7, exp_tot ~ .), data = meps, subset = pos, control = tight)
lin <- lm(update(f7, exp_tot ~ .), data = meps, subset = pos)
lnamt <- lm(update(f7, lnexp ~ .), data = meps, subset = pos)
probit <- glm(update(f7, any ~ .), family = binomial("probit"), data = meps, control = tight)

rows <- list()
add <- function(case, term, x) rows[[length(rows) + 1L]] <<-
  data.frame(case = case, term = term, r_est = x$estimate, r_se = x$std.error)
two <- function(a, b, names, ...) suest(a, b, model_names = names, twopart = TRUE, ...)
mean_dydx <- function(fit, lab, wts = FALSE) {
  add(paste0(lab, "_mean"), "_cons", avg_predictions(fit, wts = wts))
  add(paste0(lab, "_dydx"), c("age", "1.ins_unins"),
    avg_slopes(fit, variables = c("age", "ins_unins"), wts = wts))
}

# 1. suest2 page
tp <- two(any, amount, c("Any", "Amount"), cluster = "duid")
add("tp_mean", "_cons", avg_predictions(tp))
add("tp_dydx", c("age", "1.anylim", "1.female", "1.ins_unins"),
  avg_slopes(tp, variables = c("age", "female", "ins_unins", "anylim")))
add("tp_any", "1.ins_unins", avg_comparisons(tp, variables = "ins_unins", type = "positive"))
add("tp_amount", "1.ins_unins", avg_comparisons(tp, variables = "ins_unins", type = "amount"))
add("tp_amount_ifany", "1.ins_unins",
  avg_comparisons(tp, variables = "ins_unins", type = "amount", newdata = subset(exp_tot > 0)))
by_sex <- avg_comparisons(tp, variables = "ins_unins", by = "female")
add("tp_over", c("female=0", "female=1"), by_sex)
add("tp_overr", "1 vs 0", hypotheses(by_sex, hypothesis = difference ~ pairwise))

# 2. Other parts
mean_dydx(two(any, gam, c("Any", "Gamma"), cluster = "duid"), "Gamma")
mean_dydx(two(any, nb, c("Any", "NegBin"), cluster = "duid"), "NegBin")
mean_dydx(two(any, lin, c("Any", "Linear"), cluster = "duid"), "Linear")
mean_dydx(two(probit, amount, c("Probit", "Amount"), cluster = "duid"), "Probit")

# 3. Log outcome
for (r in c("duan", "normal"))
  mean_dydx(two(any, lnamt, c("Any", "LnAmount"), cluster = "duid", retransform = r), r)

# 4. Weights (rescaled to mean 1: same estimates; glm() diverges with the raw weights)
anyW <- glm(update(f7, any ~ .), family = quasibinomial, data = meps, weights = w, control = tight)
amtW <- glm(update(f7, exp_tot ~ .), family = quasipoisson, data = meps, weights = w, subset = pos,
  control = tight)
gamW <- glm(update(f7, exp_tot ~ .), family = Gamma(link = "log"), data = meps, weights = w,
  subset = pos, control = tight)
for (m in list(list(amtW, "AmountW"), list(gamW, "GammaW"))) {
  fw <- two(anyW, m[[1]], c("AnyW", m[[2]]), cluster = "duid", weight_type = "pweight")
  mean_dydx(fw, m[[2]], wts = TRUE)
  mean_dydx(fw, paste0(m[[2]], "_now"))
  rows[[length(rows) - 1L]]$case <- paste0(m[[2]], "_mean_now")
  rows[[length(rows)]]$case <- paste0(m[[2]], "_dydx_now")
}

# 5. Survey design
des <- svydesign(ids = ~duid, weights = ~wtdper, data = meps)
anyS <- svyglm(update(f7, any ~ .), design = des, family = quasibinomial)
amtS <- svyglm(update(f7, exp_tot ~ .), design = subset(des, any == 1), family = quasipoisson)
mean_dydx(suest(anyS, amtS, model_names = c("AnyS", "AmountS"), survey_design = des,
  observation_id = "dupersid", twopart = TRUE), "svy", wts = TRUE)
r <- do.call(rbind, rows)

# Comparison: rows matched in order within each case (base levels dropped)
lines <- grep("^ME\\|", readLines("test_twopart_stata_v03.log"), value = TRUE)
s <- do.call(rbind, lapply(strsplit(lines, "|", fixed = TRUE), function(x)
  data.frame(case = x[2], stata_term = x[3], s_est = as.numeric(x[4]), s_se = as.numeric(x[5]))))
s <- s[!is.na(s$s_se) & !(s$s_est == 0 & s$s_se == 0), ]
s$order <- ave(seq_len(nrow(s)), s$case, FUN = seq_along)
r$order <- ave(seq_len(nrow(r)), r$case, FUN = seq_along)
# Stata lists dydx() terms in the order given; R sorts them
tpd <- s$case == "tp_dydx"
s$order[tpd] <- match(sub("^1\\.", "", s$stata_term[tpd]), c("age", "anylim", "female", "ins_unins"))
cmp <- merge(s, r, by = c("case", "order"))
cmp$check <- paste(cmp$case, cmp$term)
cmp$rel_est <- abs(cmp$r_est - cmp$s_est) / abs(cmp$s_est)
cmp$rel_se <- abs(cmp$r_se - cmp$s_se) / cmp$s_se
print(cmp[, c("check", "s_est", "r_est", "rel_est", "s_se", "r_se", "rel_se")], row.names = FALSE,
  digits = 8)
write.csv(cmp, "comparison_twopart.csv", row.names = FALSE)
