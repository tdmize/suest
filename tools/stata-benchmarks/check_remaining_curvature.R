pkgload::load_all(quiet = TRUE)
results <- readRDS("../returned_log_comparisons.rds")
for (label in c("ML FALSE", "ML TRUE")) {
  fit <- results[[label]]$fit
  for (model in fit$models) {
    component <- suest:::.suest_lme_ml_components(model)
    frame <- suest:::.suest_model_frame(model, "nlme::lme")
    X <- model.matrix(delete.response(terms(model)), frame)
    y <- model.response(frame)
    panels <- split(seq_along(y), as.character(model$groups[[1]]))
    p <- ncol(X)
    likelihood <- function(b) {
      residual <- y - drop(X %*% b[seq_len(p)])
      su2 <- b[p + 1]^2
      se2 <- b[p + 2]^2
      sum(vapply(panels, function(rows) {
        r <- residual[rows]
        n <- length(rows)
        denominator <- se2 + n*su2
        -0.5*(n*log(2*pi) + (n - 1)*log(se2) + log(denominator) +
          sum(r^2)/se2 - su2*sum(r)^2/(se2*denominator))
      }, numeric(1)))
    }
    b <- component$parameters
    H <- optimHess(b, function(x) -likelihood(x),
      control = list(ndeps = rep(1e-4, length(b))))
    B <- solve(H)
    S <- rowsum(component$score, model$groups[[1]])
    analytic <- component$bread/nrow(frame)
    cat(label, as.character(formula(model))[2],
      "LL gap", likelihood(b) - as.numeric(logLik(model)),
      "bread gap", max(abs(B - analytic)),
      "raw robust gap", max(abs(B %*% crossprod(S) %*% B -
        analytic %*% crossprod(S) %*% analytic)), "\n")
  }
}
fit <- readRDS("../endogenous_comparison_fits.rds")$bp
influences <- list()
for (model in fit$models) {
  frame <- suest:::.suest_model_frame(model, "mvProbit::mvProbit")
  X <- model.matrix(delete.response(terms(formula(model))), frame)
  y <- model.response(frame)
  signs <- 2*y - 1
  k <- ncol(X)
  component <- suest:::.suest_biprobit_components(model)
  score <- function(b) {
    e1 <- drop(X %*% b[seq_len(k)])
    e2 <- drop(X %*% b[k + seq_len(k)])
    rho <- tanh(tail(b, 1))
    s <- sqrt(1 - rho^2)
    probability <- vapply(seq_along(e1), function(i) {
      r <- prod(signs[i, ])*rho
      as.numeric(mvtnorm::pmvnorm(upper = signs[i, ]*c(e1[i], e2[i]),
        corr = matrix(c(1, r, r, 1), 2), algorithm = mvtnorm::TVPACK()))
    }, numeric(1))
    a <- signs[, 1]*dnorm(e1)*pnorm(signs[, 2]*(e2 - rho*e1)/s)/probability
    b2 <- signs[, 2]*dnorm(e2)*pnorm(signs[, 1]*(e1 - rho*e2)/s)/probability
    density <- exp(-(e1^2 - 2*rho*e1*e2 + e2^2)/(2*s^2))/(2*pi*s)
    cbind(X*a, X*b2, signs[, 1]*signs[, 2]*density*s^2/probability)
  }
  b <- component$parameters
  H <- sapply(seq_along(b), function(j) {
    plus <- minus <- b
    plus[j] <- plus[j] + 1e-5
    minus[j] <- minus[j] - 1e-5
    (colSums(score(plus)) - colSums(score(minus)))/2e-5
  })
  B <- solve(-(H + t(H))/2)
  cat("biprobit", k, "score gap", max(abs(score(b) - component$score)),
    "bread gap", max(abs(B - component$bread/nrow(X))), "\n")
  influences[[length(influences) + 1L]] <- component$score %*% B
}
I <- do.call(cbind, influences)
V <- crossprod(I) * nrow(I)/(nrow(I) - 1)
order <- c(2, 1, 4, 3, 5, 7, 8, 6, 10, 11, 9, 12)
cat("Independent biprobit full V versus Stata:",
  max(abs(V[order, order] - results$biprobit$Stata$V)), "\n")
