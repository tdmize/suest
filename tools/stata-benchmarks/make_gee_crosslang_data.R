args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[1] else "suest_r_gee_crosslang_benchmark.csv"
set.seed(7630)
id <- rep(1:120, each = 4)
time <- rep(1:4, 120)
x <- rnorm(length(id))
z <- rnorm(length(id))
u <- rnorm(120, sd = 0.5)[id]
dat <- data.frame(id, time, higher = ceiling(id/4), x, z,
  y = 0.4 + 0.6*x - 0.3*z + u + rnorm(length(id)),
  binary = rbinom(length(id), 1, plogis(-0.2 + 0.4*x - 0.2*z + u)),
  count = rpois(length(id), exp(0.2 + 0.3*x - 0.15*z + 0.3*u)))
write.csv(dat, output, row.names = FALSE, quote = FALSE)
