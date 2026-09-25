glmmtmb_logit_rs_data <- function(groups = 100L, periods = 12L) {
  set.seed(2437)
  id <- rep(seq_len(groups), each = periods); x <- runif(length(id), -2, 2); z <- rnorm(length(id))
  a <- rnorm(groups, sd = .9); b <- .35*a + rnorm(groups, sd = .7)
  a2 <- .55*a + rnorm(groups, sd = .75); b2 <- .5*b + rnorm(groups, sd = .75)
  woman <- rep(rep(0:1, length.out = groups), each = periods)
  data.frame(id = factor(id), time = rep(seq_len(periods), groups), higher = ceiling(id/5),
    x = x, z = z, woman = factor(woman),
    y1 = rbinom(length(id), 1, plogis(-.3 + .55*x - .25*z + .2*woman + a[id] + b[id]*x)),
    y2 = rbinom(length(id), 1, plogis(-.15 + .4*x + .2*z - .15*woman + a2[id] + b2[id]*x)))
}

dat <- glmmtmb_logit_rs_data(); dat$id <- as.integer(dat$id)
dat$y1_partial <- ifelse(dat$time == 12 & dat$id %% 4 == 0, NA, dat$y1)
dat$y2_partial <- ifelse(dat$time == 1 & dat$id %% 5 == 0, NA, dat$y2)
dat$y1_left <- ifelse(dat$id %% 2 == 1, dat$y1, NA)
dat$y2_right <- ifelse(dat$id %% 2 == 0, dat$y2, NA)
write.csv(dat, "suest_r_glmmtmb_logit_rs_crosslang_benchmark.csv", row.names = FALSE, na = "")
cat("LOGIT_RS_FIXED_DATA_GENERATED=1\n")
