args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L)
  stop("Supply at most one output .csv path.", call. = FALSE)
output <- if (length(args)) args[[1L]] else
  "suest_r_iv_ancillary_benchmark.csv"

set.seed(3371)
n <- 1200L
x <- rnorm(n)
z1 <- rnorm(n)
z2 <- rnorm(n)
disturbance <- rnorm(n)
endogenous <- 0.7*z1 + 0.25*z2 + 0.4*disturbance + rnorm(n)
y <- 0.45 + 0.65*x + 1.15*endogenous + disturbance
hetero_probability <- pnorm((0.2 + 0.7*x - 0.25*z2)/exp(0.35*z1))
y_het <- rbinom(n, 1, hetero_probability)

data <- data.frame(
  id = seq_len(n),
  y = y,
  x = x,
  endogenous = endogenous,
  z1 = z1,
  z2 = z2,
  y_het = y_het
)

utils::write.csv(
  data,
  file = output,
  row.names = FALSE,
  quote = FALSE
)
