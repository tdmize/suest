set.seed(10917)
groups <- 80L
periods <- 6L
id <- rep(seq_len(groups), each = periods)
time <- rep(seq_len(periods), groups)
higher <- ceiling(id/5)
x <- rnorm(groups*periods)
z <- rnorm(groups*periods)
common <- rnorm(groups, sd = 0.55)[id]
second <- rnorm(groups, sd = 0.35)[id]
y1 <- rpois(length(id), exp(0.25 + 0.35*x - 0.20*z + common))
y2 <- rpois(length(id), exp(0.10 + 0.20*x + 0.25*z + 0.7*common + second))

dat <- data.frame(id, time, higher, x, z, y1, y2)
dat$y1_partial <- ifelse(time == 6 & id %% 4 == 0, NA, y1)
dat$y2_partial <- ifelse(time == 1 & id %% 5 == 0, NA, y2)
dat$y1_left <- ifelse(id <= 40, y1, NA)
dat$y2_right <- ifelse(id > 40, y2, NA)

write.csv(dat, "suest_r_glmmtmb_poisson_ri_crosslang_benchmark.csv",
  row.names = FALSE, na = "")
cat("GLMMTMB_POISSON_RI_DATA_COMPLETE=1\n")
