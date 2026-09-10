set.seed(29107)
groups <- 80L
periods <- 6L
id <- rep(seq_len(groups), each = periods)
time <- rep(seq_len(periods), groups)
higher <- ceiling(id / 5)
x <- rnorm(length(id))
z <- rnorm(length(id))
common <- rnorm(groups, sd = 0.85)[id]
second <- rnorm(groups, sd = 0.45)[id]
y1 <- rbinom(length(id), 1, pnorm(-0.25 + 0.55*x - 0.30*z + common))
y2 <- rbinom(length(id), 1, pnorm(0.15 + 0.40*x + 0.25*z +
  0.70*common + second))

dat <- data.frame(id, time, higher, x, z, y1, y2)
dat$y1_partial <- dat$y1
dat$y2_partial <- dat$y2
dat$y1_partial[dat$time == 6 & dat$id %% 4 == 0] <- NA
dat$y2_partial[dat$time == 1 & dat$id %% 5 == 0] <- NA
dat$y1_left <- ifelse(dat$id <= groups/2, dat$y1, NA)
dat$y2_right <- ifelse(dat$id > groups/2, dat$y2, NA)

write.csv(
  dat,
  "suest_r_panel_probit_re_crosslang_benchmark.csv",
  row.names = FALSE,
  na = ""
)
