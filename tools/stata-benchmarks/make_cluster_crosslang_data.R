args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L)
  stop("Supply at most one output .csv path.", call. = FALSE)
output <- if (length(args)) args[[1L]] else
  "suest_r_cluster_crosslang_benchmark.csv"

set.seed(6109)
n <- 1200L
id <- seq_len(n)
cluster4 <- ceiling(id / 4)
cluster2 <- ceiling(id / 2)
side <- id %% 2L
cluster_effect <- rnorm(max(cluster4))[cluster4]
x <- rnorm(n) + 0.25*cluster_effect
z <- rnorm(n)
z1 <- rnorm(n)
z2 <- rnorm(n)
disturbance <- rnorm(n)
endogenous <- 0.75*z1 + 0.20*z2 + 0.35*disturbance + rnorm(n)
y <- 1 + 0.45*x - 0.30*z + 0.60*cluster_effect + rnorm(n)
y_iv <- 0.5 + 0.6*x + 1.1*endogenous + disturbance
binary_probability <- plogis(
  -0.25 + 0.50*x - 0.20*z + 0.45*cluster_effect
)
b <- rbinom(n, 1, binary_probability)
pw <- exp(0.20*z + 0.10*cluster_effect)

data <- data.frame(
  id = id,
  cluster4 = cluster4,
  cluster2 = cluster2,
  side = side,
  x = x,
  z = z,
  z1 = z1,
  z2 = z2,
  endogenous = endogenous,
  y = y,
  y_iv = y_iv,
  b = b,
  pw = pw
)

utils::write.csv(
  data,
  file = output,
  row.names = FALSE,
  quote = FALSE
)
