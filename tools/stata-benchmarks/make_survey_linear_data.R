# Run from the package root: source("tools/stata-benchmarks/make_survey_linear_data.R")
# No random draws: the CSV is identical across R platforms.
i <- 1:240
psu <- rep(1:24, each = 10)
strata <- rep(1:4, each = 60)
unit <- rep(1:10, 24)
x <- sin(i * .71) + cos(psu * .37)
z <- cos(i * .43) + sin(psu * .29)
y <- 2 + .7 * x - .4 * z + sin(psu * 1.13) + cos(i * 1.77)
y2 <- 1 - .2 * x + .8 * z + .6 * y + sin(i * 1.31)
d <- data.frame(id = i, strata, psu, unit, x, z, y, y2,
  w = 1 + (i %% 9) / 4, pop = rep(c(12, 18, 24, 30), each = 60))
write.csv(d, "tools/stata-benchmarks/suest_r_survey_linear_benchmark.csv", row.names = FALSE)
