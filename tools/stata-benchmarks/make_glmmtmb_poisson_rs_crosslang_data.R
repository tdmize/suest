set.seed(2431)
groups <- 100L; periods <- 8L
id <- rep(seq_len(groups), each = periods)
x <- runif(length(id), -1.5, 1.5); z <- rnorm(length(id))
a <- rnorm(groups, sd = .5); b <- .3*a + rnorm(groups, sd = .4)
a2 <- .6*a + rnorm(groups, sd = .4)
b2 <- .5*b + rnorm(groups, sd = .35)
dat <- data.frame(id = id, time = rep(seq_len(periods), groups), higher = ceiling(id/5),
  x = x, z = z,
  y1 = rpois(length(id), exp(.6 + .3*x - .2*z + a[id] + b[id]*x)),
  y2 = rpois(length(id), exp(.4 + .2*x + .15*z + a2[id] + b2[id]*x)))
dat$y1_partial <- dat$y1; dat$y2_partial <- dat$y2
dat$y1_partial[dat$time == 8 & dat$id %% 4 == 0] <- NA
dat$y2_partial[dat$time == 1 & dat$id %% 5 == 0] <- NA
dat$y1_left <- ifelse(dat$id %% 2 == 1, dat$y1, NA)
dat$y2_right <- ifelse(dat$id %% 2 == 0, dat$y2, NA)
write.csv(dat, "suest_r_glmmtmb_poisson_rs_crosslang_benchmark.csv", row.names = FALSE, na = "")
cat("GLMMTMB_POISSON_RS_DATA_COMPLETE=1\n")
