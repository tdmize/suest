# Independent joint IV-probit likelihood, not Rchoice's score/Hessian methods.
pkgload::load_all(quiet = TRUE)
fits <- readRDS("../endogenous_comparison_fits.rds")
for (model in fits$iv$models) {
  X <- model.matrix(model$formula, model$mf, rhs = 1)
  Z <- model.matrix(model$formula, model$mf, rhs = 2)
  y <- model$mf[[1]]
  endogenous <- X[, "endogenous"]
  k <- ncol(X)
  p <- ncol(Z)
  parameters <- coef(model)
  likelihood <- function(b) {
    sigma <- exp(b[k + p + 1])
    rho <- tanh(b[k + p + 2])
    u <- (endogenous - drop(Z %*% b[k + seq_len(p)]))/sigma
    index <- (drop(X %*% b[seq_len(k)]) + rho*u)/sqrt(1 - rho^2)
    dnorm(u, log = TRUE) - log(sigma) + pnorm((2*y - 1)*index, log.p = TRUE)
  }
  U <- sapply(seq_along(parameters), function(j) {
    plus <- minus <- parameters
    plus[j] <- plus[j] + 1e-5
    minus[j] <- minus[j] - 1e-5
    (likelihood(plus) - likelihood(minus))/2e-5
  })
  cat("LL difference:", sum(likelihood(parameters)) - as.numeric(logLik(model)), "\n")
  cat("Score difference:", max(abs(U - sandwich::estfun(model))), "\n")
  for (step in c(1e-4, 1e-5)) {
    H <- optimHess(parameters, function(b) -sum(likelihood(b)),
      control = list(ndeps = rep(step, length(parameters))))
    cat("Hessian difference:", step, max(abs(H + model$hessian)), "\n")
    B <- solve(H)
    cat("Bread difference:", max(abs(B - sandwich::bread(model)/nrow(X))), "\n")
  }
}
