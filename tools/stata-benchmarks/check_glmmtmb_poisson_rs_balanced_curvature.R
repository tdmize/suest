# Run from tools/stata-benchmarks. Optional arguments: fixture, output RDS.
args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args)) args[1] else
  "../../tests/testthat/fixtures/glmmtmb-poisson-rs-stata-balanced-v4.rds"
s <- readRDS(input)
f <- readRDS("../../tests/testthat/fixtures/glmmtmb-poisson-rs-stata-components.rds")
r <- readRDS("glmmtmb_poisson_rs_r_reference.rds")$systems$balanced
d <- f$data; groups <- split(seq_len(nrow(d)), d$id)
X <- cbind(1, d$x, d$z); Z <- cbind(1, d$x)
stopifnot(s$nobs == 800L, s$n_clusters == 100L, length(groups) == 100L,
  identical(names(s$coefficients), names(r$coefficients)),
  !is.null(s$modelbased), all(is.finite(s$modelbased)))
# Keep this evaluator independent of the package's TMB score implementation.
# Independent conditional-mode solver and bivariate Laplace group likelihood.
group_ll <- function(b,y) {
  sd <- exp(b[4:5]); rho <- tanh(b[6]); S <- outer(sd,sd)*matrix(c(1,rho,rho,1),2)
  P <- solve(S); eta <- drop(X%*%b[1:3]); logdet <- as.numeric(determinant(S,logarithm=TRUE)$modulus)
  vapply(groups,function(idx) {
    zi <- Z[idx,,drop=FALSE]; yi <- y[idx]; ei <- eta[idx]
    objective <- function(u) sum(exp(ei+drop(zi%*%u))-yi*(ei+drop(zi%*%u)))+sum(u*(P%*%u))/2
    u <- c(0,0)
    for(k in 1:100) {
      mu <- exp(ei+drop(zi%*%u)); g <- drop(crossprod(zi,mu-yi)+P%*%u)
      if(max(abs(g)) < 1e-10) break
      H <- crossprod(zi,zi*mu)+P; step <- drop(solve(H,g)); scale <- 1
      # Permit rounding-level equality near the mode; gradient verifies accuracy.
      while(objective(u-scale*step) > objective(u)+1e-12 && scale>1e-8) scale <- scale/2
      u <- u-scale*step
    }
    mu <- exp(ei+drop(zi%*%u)); g <- drop(crossprod(zi,mu-yi)+P%*%u)
    stopifnot(max(abs(g))<1e-8)
    H <- crossprod(zi,zi*mu)+P
    sum(dpois(yi,mu,log=TRUE))-sum(u*(P%*%u))/2-logdet/2-as.numeric(determinant(H,logarithm=TRUE)$modulus)/2
  },numeric(1))
}
# Evaluate at the returned JOINT parameters, so optimizer differences cannot
# explain the decomposition. Independent curvature is block diagonal because
# all cross-equation latent covariances are fixed at zero.
results <- lapply(1:2, function(i) {
  b <- s$coefficients[(i-1L)*6L+1:6]; y <- d[[c("y1", "y2")[i]]]
  scores <- lapply(c(1e-5, 5e-6), function(h) sapply(seq_along(b), function(j) {
    plus <- minus <- b; plus[j] <- plus[j]+h; minus[j] <- minus[j]-h
    (group_ll(plus, y)-group_ll(minus, y))/(2*h)
  }))
  hessian <- function(h) {
    H <- matrix(0, 6, 6); f0 <- sum(group_ll(b, y))
    for (j in 1:6) for (k in j:6) {
      ej <- ek <- numeric(6); ej[j] <- h; ek[k] <- h
      val <- if (j == k) {
        (sum(group_ll(b+ej, y))-2*f0+sum(group_ll(b-ej, y)))/h^2
      } else {
        (sum(group_ll(b+ej+ek, y))-sum(group_ll(b+ej-ek, y))-
          sum(group_ll(b-ej+ek, y))+sum(group_ll(b-ej-ek, y)))/(4*h^2)
      }
      H[j, k] <- H[k, j] <- val
    }
    H
  }
  # Cancel the leading O(h^2) error, then halve both input steps to check it.
  H <- lapply(c(.002, .001, .0005), hessian)
  bread1 <- solve(-(4*H[[2]]-H[[1]])/3)
  bread2 <- solve(-(4*H[[3]]-H[[2]])/3)
  score_change <- max(abs(scores[[1]]-scores[[2]]))
  bread_change <- max(abs(bread1-bread2))
  stopifnot(score_change < 1e-7, bread_change < 1e-7,
    min(eigen(bread2, symmetric = TRUE)$values) > 0)
  cat("Outcome", i, "score step change=", score_change,
    "bread step change=", bread_change, "\n")
  # Away from an exact optimum, observed information depends on coordinates.
  # Sensitivity: transform natural (beta, A, B, C) information to these
  # log-SD/Fisher coordinates, retaining the total-score chain-rule term.
  g <- colSums(scores[[2]]); rho <- tanh(b[6]); q <- rho/(1-rho^2)
  T <- matrix(0, 6, 6)
  T[4, 4] <- 2*g[4]-q*g[6]; T[5, 5] <- 2*g[5]-q*g[6]
  T[4, 5] <- T[5, 4] <- q*g[6]
  T[4, 6] <- T[6, 4] <- T[5, 6] <- T[6, 5] <- g[6]
  T[6, 6] <- -2*rho*g[6]
  bread_natural_delta <- solve(solve(bread2)+T)
  list(scores = scores[[2]], bread = bread2, gradient = g,
    bread_natural_delta = bread_natural_delta,
    logLik = sum(group_ll(b, y)), score_step_change = score_change,
    bread_step_change = bread_change)
})
S <- cbind(results[[1]]$scores, results[[2]]$scores)
B <- matrix(0, 12, 12)
B[1:6, 1:6] <- results[[1]]$bread; B[7:12, 7:12] <- results[[2]]$bread
Bn <- matrix(0, 12, 12)
Bn[1:6, 1:6] <- results[[1]]$bread_natural_delta
Bn[7:12, 7:12] <- results[[2]]$bread_natural_delta
correction <- s$n_clusters/(s$n_clusters-1)
M <- crossprod(S); Bs <- s$modelbased; Vs <- s$vcov
stopifnot(min(eigen(Bs, symmetric = TRUE)$values) > 0)
# V = correction * B M B', with symmetric bread B. Undo both bread factors.
Ms <- solve(Bs, Vs) %*% solve(Bs)/correction
Vind <- B %*% M %*% B*correction
Vhybrid <- Bs %*% M %*% Bs*correction
Vreplace <- B %*% Ms %*% B*correction
metrics <- function(A, reference) c(max_abs_covariance = max(abs(A-reference)),
  max_relative_diagonal_sd = max(abs(sqrt(diag(A)/diag(reference))-1)),
  max_scaled_covariance = max(abs(A-reference)/sqrt(outer(diag(reference), diag(reference)))))
comparison <- rbind(
  bread_stata_vs_independent = metrics(Bs, B),
  meat_recovered_stata_vs_independent = metrics(Ms, M),
  sandwich_independent_vs_R = metrics(Vind, r$vcov),
  stata_bread_independent_meat_vs_stata = metrics(Vhybrid, Vs),
  independent_bread_stata_meat_vs_independent = metrics(Vreplace, Vind),
  raw_stata_vs_independent = metrics(Vs, Vind),
  coordinate_bread_sensitivity = metrics(Bn, B),
  coordinate_sandwich_sensitivity = metrics(Bn %*% M %*% Bn*correction, Vind))
print(comparison, digits = 10)
ll_difference <- sum(vapply(results, `[[`, numeric(1), "logLik"))-s$logLik
cat("INDEPENDENT_MINUS_STATA_LOGLIK=", format(ll_difference, digits = 16), "\n")
cat("STATA_MODELBASED_MAX_CROSS_BLOCK=", max(abs(Bs[1:6, 7:12])), "\n")
# The independent audit retains the earlier 1e-6 covariance / 1e-4 SE limits.
# Raw joint GSEM covariance is a diagnostic, not an accepted equality fixture.
stopifnot(abs(ll_difference) < 1e-5,
  max(abs(Bs %*% Ms %*% Bs*correction-Vs)) < 1e-12,
  comparison["sandwich_independent_vs_R", "max_abs_covariance"] < 1e-6,
  comparison["sandwich_independent_vs_R", "max_relative_diagonal_sd"] < 1e-4,
  comparison["stata_bread_independent_meat_vs_stata", "max_abs_covariance"] < 1e-6,
  comparison["stata_bread_independent_meat_vs_stata", "max_relative_diagonal_sd"] < 1e-4)
if (length(args) >= 2L) saveRDS(list(components = results, comparison = comparison,
  independent_bread = B, natural_delta_bread = Bn, independent_meat = M, stata_bread = Bs,
  recovered_stata_meat = Ms, independent_vcov = Vind,
  stata_bread_independent_meat_vcov = Vhybrid,
  independent_bread_stata_meat_vcov = Vreplace,
  metadata = list(source_fixture = basename(input),
    source_md5 = unname(tools::md5sum(input)), correction = correction,
    logLik_difference = ll_difference)), args[2])
cat("POISSON_RS_BALANCED_CURVATURE_AUDIT_COMPLETE=1\n")
