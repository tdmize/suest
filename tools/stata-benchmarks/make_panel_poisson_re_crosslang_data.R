set.seed(41903)
groups <- 80L
periods <- 6L
id <- rep(seq_len(groups), each = periods)
time <- rep(seq_len(periods), groups)
higher <- ceiling(id / 5)
x <- rnorm(length(id))
z <- rnorm(length(id))
common <- rgamma(groups, shape = 2.2, rate = 2.2)[id]
second <- rgamma(groups, shape = 3.0, rate = 3.0)[id]
mu1 <- exp(0.20 + 0.35*x - 0.20*z)*common
mu2 <- exp(-0.10 + 0.25*x + 0.20*z)*sqrt(common)*second
y1 <- rpois(length(id), mu1)
y2 <- rpois(length(id), mu2)

dat <- data.frame(id, time, higher, x, z, y1, y2)
dat$y1_partial <- dat$y1
dat$y2_partial <- dat$y2
dat$y1_partial[dat$time == 6 & dat$id %% 4 == 0] <- NA
dat$y2_partial[dat$time == 1 & dat$id %% 5 == 0] <- NA
dat$y1_left <- ifelse(dat$id <= groups/2, dat$y1, NA)
dat$y2_right <- ifelse(dat$id > groups/2, dat$y2, NA)

write.csv(
  dat,
  "suest_r_panel_poisson_re_crosslang_benchmark.csv",
  row.names = FALSE,
  na = ""
)
