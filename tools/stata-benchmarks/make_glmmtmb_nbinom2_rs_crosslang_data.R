glmmtmb_nbinom2_rs_data <- function(groups = 100L, periods = 10L) {
  set.seed(2432)
  id <- rep(seq_len(groups), each = periods)
  x <- runif(length(id), -1.5, 1.5); z <- rnorm(length(id))
  a <- rnorm(groups, sd = .6); b <- .35*a + rnorm(groups, sd = .45)
  a2 <- .6*a + rnorm(groups, sd = .45); b2 <- .6*b + rnorm(groups, sd = .4)
  woman <- rep(rep(0:1, length.out = groups), each = periods)
  data.frame(id = factor(id), time = rep(seq_len(periods), groups), higher = ceiling(id/5),
    x = x, z = z, woman = factor(woman),
    y1 = rnbinom(length(id), size = 2.5, mu = exp(.8 + .3*x - .2*z + .2*woman + a[id] + b[id]*x)),
    y2 = rnbinom(length(id), size = 3.5, mu = exp(.7 + .2*x + .15*z - .1*woman + a2[id] + b2[id]*x)))
}

dat <- glmmtmb_nbinom2_rs_data(); dat$id <- as.integer(dat$id)
dat$y1_partial <- ifelse(dat$time == 10 & dat$id %% 4 == 0, NA, dat$y1)
dat$y2_partial <- ifelse(dat$time == 1 & dat$id %% 5 == 0, NA, dat$y2)
dat$y1_left <- ifelse(dat$id %% 2 == 1, dat$y1, NA)
dat$y2_right <- ifelse(dat$id %% 2 == 0, dat$y2, NA)
write.csv(dat, "suest_r_glmmtmb_nbinom2_rs_crosslang_benchmark.csv", row.names = FALSE, na = "")
cat("NB2_RS_FIXED_DATA_GENERATED=1\n")
