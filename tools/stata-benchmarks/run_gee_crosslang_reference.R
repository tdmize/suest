# Run from the repository root; supply the benchmark CSV path as argument 1.
args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args)) args[1] else "suest_r_gee_crosslang_benchmark.csv"
pkgload::load_all(quiet = TRUE)
dat <- read.csv(input)
families <- list(gaussian(), binomial(), binomial("probit"),
  binomial("cloglog"), poisson())
outcomes <- c("y", "binary", "binary", "binary", "count")
for (correlation in c("independence", "exchangeable")) {
  for (j in seq_along(families)) {
    cat("\nGEE CASE:", outcomes[j], families[[j]]$link, correlation, "\n")
    first <- geepack::geeglm(reformulate("x", outcomes[j]), id = id,
      data = dat, family = families[[j]], corstr = correlation,
      control = geepack::geese.control(epsilon = 1e-10, maxit = 100))
    second <- geepack::geeglm(reformulate(c("x", "z"), outcomes[j]), id = id,
      data = dat, family = families[[j]], corstr = correlation,
      control = geepack::geese.control(epsilon = 1e-10, maxit = 100))
    fit <- suest(first, second, model_names = c("Base", "Adjusted"))
    print(coef(fit), digits = 16)
    print(vcov(fit), digits = 16)
  }
}
cat("\nGEE CASE: MIXED FAMILY, UNEQUAL SAMPLES, HIGHER CLUSTERS\n")
first <- geepack::geeglm(y ~ x + z, id = id, data = dat, subset = id <= 100,
  corstr = "exchangeable")
second <- geepack::geeglm(binary ~ x + z, id = id, data = dat,
  subset = id >= 21 & time != 4, family = binomial(), corstr = "exchangeable")
fit <- suest(first, second, model_names = c("Gaussian", "Logit"), cluster = "higher")
print(coef(fit), digits = 16)
print(vcov(fit), digits = 16)
