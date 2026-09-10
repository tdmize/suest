set.seed(10916)
groups <- 80L
periods <- 6L
id <- rep(seq_len(groups), each = periods)
time <- rep(seq_len(periods), groups)
higher <- ceiling(id/5)
x <- rnorm(length(id))
z <- rnorm(length(id))
common <- rnorm(groups, sd = 1.05)[id]
second <- rnorm(groups, sd = 0.60)[id]
y1 <- rbinom(length(id), 1,
  plogis(-0.30 + 0.65*x - 0.25*z + common))
y2 <- rbinom(length(id), 1,
  plogis(0.10 + 0.45*x + 0.20*z + 0.70*common + second))

dat <- data.frame(id, time, higher, x, z, y1, y2)
dat$y1_partial <- dat$y1
dat$y2_partial <- dat$y2
dat$y1_partial[dat$time == 6 & dat$id %% 4 == 0] <- NA
dat$y2_partial[dat$time == 1 & dat$id %% 5 == 0] <- NA
dat$y1_left <- ifelse(dat$id <= groups/2, dat$y1, NA)
dat$y2_right <- ifelse(dat$id > groups/2, dat$y2, NA)

write.csv(dat, "suest_r_glmmtmb_logit_ri_crosslang_benchmark.csv",
  row.names = FALSE, na = "")
