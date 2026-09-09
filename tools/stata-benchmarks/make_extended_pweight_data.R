args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L)
  stop("Supply the input .csv.gz and output .csv paths.", call. = FALSE)

dat <- utils::read.csv(gzfile(args[1L]))
dat$y_nb <- dat$y_count + 3L*(dat$id %% 11L == 0L) +
  5L*(dat$id %% 37L == 0L)
ordered_rank <- rank(dat$y_linear, ties.method = "first")
dat$y_ord <- as.integer(cut(
  ordered_rank,
  breaks = c(0, 600, 1200, 1800, 2400),
  labels = FALSE
))
utils::write.csv(dat, args[2L], row.names = FALSE)
