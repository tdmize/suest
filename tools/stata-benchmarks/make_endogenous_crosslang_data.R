args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L)
  stop("Supply at most one output .csv path.", call. = FALSE)
output <- if (length(args)) args[[1L]] else
  "suest_r_endogenous_crosslang_benchmark.csv"

set.seed(11843)
n <- 600L
x <- rnorm(n)
z <- rnorm(n)
z1 <- rnorm(n)
z2 <- rnorm(n)

bp_common <- rnorm(n)
bp1 <- 0.45*bp_common + sqrt(1 - 0.45^2)*rnorm(n)
bp2 <- 0.45*bp_common + sqrt(1 - 0.45^2)*rnorm(n)
binary1 <- as.integer(0.20 + 0.65*x - 0.20*z + bp1 > 0)
binary2 <- as.integer(-0.15 + 0.35*x + 0.25*z + bp2 > 0)

iv_common <- rnorm(n)
outcome_error <- 0.35*iv_common + sqrt(1 - 0.35^2)*rnorm(n)
first_stage_error <- iv_common
endogenous <- 0.25*x + 0.70*z1 + 0.20*z2 + first_stage_error
binary_iv <- as.integer(
  0.15 + 0.50*x + 0.65*endogenous + outcome_error > 0
)

write.csv(
  data.frame(
    observation = seq_len(n), x, z, z1, z2,
    binary1, binary2, endogenous, binary_iv
  ),
  output,
  row.names = FALSE,
  quote = FALSE
)
