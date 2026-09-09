args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L)
  stop("Supply at most one output .csv path.", call. = FALSE)
output <- if (length(args)) args[[1L]] else
  "suest_r_panel_fe_crosslang_benchmark.csv"

set.seed(7621)
groups <- 120L
periods <- 5L
id <- rep(seq_len(groups), each = periods)
time <- rep(seq_len(periods), groups)
higher <- ceiling(id / 4)
person1 <- rnorm(groups)[id]
person2 <- rnorm(groups)[id]
common <- rnorm(max(higher))[higher]
x <- rnorm(length(id)) + 0.2*person1
z <- rnorm(length(id))
y1 <- 1.1 + 0.55*x - 0.30*z + 0.15*time + person1 + 0.4*common + rnorm(length(id))
y2 <- -0.4 + 0.35*x + 0.20*z - 0.08*time + person2 + 0.3*person1 + rnorm(length(id))
y1_partial <- y1
y2_partial <- y2
y1_partial[time == 5 & id %% 4 == 0] <- NA
y2_partial[time == 1 & id %% 5 == 0] <- NA

data <- data.frame(id, time, higher, x, z, y1, y2, y1_partial, y2_partial)
utils::write.csv(data, file = output, row.names = FALSE, quote = FALSE, na = "")
