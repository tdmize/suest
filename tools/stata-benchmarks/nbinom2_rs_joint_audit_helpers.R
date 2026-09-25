# Independent joint covariance audit; no package/TMB score code is used here.
audit_nbinom2_rs_joint <- function(s, r, d, outcomes, cluster, components, evaluate_at = c("Stata", "R")) {
  evaluate_at <- match.arg(evaluate_at)
  evaluation_parameters <- if (evaluate_at == "Stata") s$coefficients else r$coefficients
  groups <- split(seq_len(nrow(d)), d$id)
  X <- cbind(1, d$x, d$z); Z <- cbind(1, d$x)
  stopifnot(s$nobs == nrow(d), identical(names(s$coefficients), names(r$coefficients)),
    !is.null(s$modelbased), all(is.finite(s$modelbased)),
    all(vapply(groups, function(idx) length(unique(d[[cluster]][idx])) == 1L, logical(1))))
  # Keep this evaluator independent of the package's TMB score implementation.
  # Independent conditional-mode solver and bivariate Laplace group likelihood.
  group_ll <- function(b,y) {
    phi <- exp(b[4]); sd <- exp(b[5:6]); rho <- tanh(b[7]); S <- outer(sd,sd)*matrix(c(1,rho,rho,1),2)
    P <- solve(S); eta <- drop(X%*%b[1:3]); logdet <- as.numeric(determinant(S,logarithm=TRUE)$modulus)
    vapply(groups,function(idx) {
      idx <- idx[!is.na(y[idx])]
      if (!length(idx)) return(0)
      zi <- Z[idx,,drop=FALSE]; yi <- y[idx]; ei <- eta[idx]
      objective <- function(u) -sum(dnbinom(yi, size=phi, mu=exp(ei+drop(zi%*%u)), log=TRUE))+sum(u*(P%*%u))/2
      u <- c(0,0)
      for(k in 1:100) {
        mu <- exp(ei+drop(zi%*%u)); g <- drop(crossprod(zi,phi*(mu-yi)/(phi+mu))+P%*%u)
        if(max(abs(g)) < 1e-10) break
        H <- crossprod(zi,zi*as.numeric(phi*mu*(phi+yi)/(phi+mu)^2))+P; step <- drop(solve(H,g)); scale <- 1
        # Permit rounding-level equality near the mode; gradient verifies accuracy.
        while(objective(u-scale*step) > objective(u)+1e-12 && scale>1e-8) scale <- scale/2
        u <- u-scale*step
      }
      mu <- exp(ei+drop(zi%*%u)); g <- drop(crossprod(zi,phi*(mu-yi)/(phi+mu))+P%*%u)
      stopifnot(max(abs(g))<1e-8)
      H <- crossprod(zi,zi*as.numeric(phi*mu*(phi+yi)/(phi+mu)^2))+P
      sum(dnbinom(yi,size=phi,mu=mu,log=TRUE))-sum(u*(P%*%u))/2-logdet/2-as.numeric(determinant(H,logarithm=TRUE)$modulus)/2
    },numeric(1))
  }
  # Audit at Stata's joint parameters for the decomposition, then at R's
  # parameters for a like-for-like adapter check. Curvature is block diagonal because
  # all cross-equation latent covariances are fixed at zero.
  results <- lapply(1:2, function(i) {
    b <- evaluation_parameters[(i-1L)*7L+1:7]; y <- d[[outcomes[i]]]
    scores <- lapply(c(1e-5, 5e-6), function(h) sapply(seq_along(b), function(j) {
      plus <- minus <- b; plus[j] <- plus[j]+h; minus[j] <- minus[j]-h
      (group_ll(plus, y)-group_ll(minus, y))/(2*h)
    }))
    hessian <- function(h) {
      H <- matrix(0, 7, 7); f0 <- sum(group_ll(b, y))
      for (j in 1:7) for (k in j:7) {
        ej <- ek <- numeric(7); ej[j] <- h; ek[k] <- h
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
    # Sensitivity: transform natural (beta, lnalpha, A, B, C) information to these
    # log-SD/Fisher coordinates, retaining the total-score chain-rule term.
    g <- colSums(scores[[2]]); rho <- tanh(b[7]); q <- rho/(1-rho^2)
    T <- matrix(0, 7, 7)
    T[5, 5] <- 2*g[5]-q*g[7]; T[6, 6] <- 2*g[6]-q*g[7]
    T[5, 6] <- T[6, 5] <- q*g[7]
    T[5, 7] <- T[7, 5] <- T[6, 7] <- T[7, 6] <- g[7]
    T[7, 7] <- -2*rho*g[7]
    bread_natural_delta <- solve(solve(bread2)+T)
    # Native glmmTMB correlation coordinate t = sinh(z), z = atanh(rho).
    # The delta-transformed inverse information is (I_z + g_z*rho*E77)^-1:
    # d2z/dt2 = -rho*(dz/dt)^2. Retain this score term away from an exact optimum.
    Tnative <- matrix(0, 7, 7); Tnative[7, 7] <- g[7]*rho
    bread_native_delta <- solve(solve(bread2)+Tnative)
    list(scores = scores[[2]], bread = bread2, gradient = g,
      bread_natural_delta = bread_natural_delta,
      bread_native_delta = bread_native_delta,
      logLik = sum(group_ll(b, y)), score_step_change = score_change,
      bread_step_change = bread_change)
  })
  S <- cbind(results[[1]]$scores, results[[2]]$scores)
  group_clusters <- vapply(groups, function(idx) as.character(d[[cluster]][idx[1L]]), character(1))
  S <- rowsum(S, group_clusters, reorder = FALSE)
  stopifnot(nrow(S) == s$n_clusters)
  component_clusters <- vapply(outcomes, function(y) length(unique(d[[cluster]][!is.na(d[[y]])])), integer(1))
  stopifnot(length(unique(component_clusters)) == 1L)
  r_correction <- component_clusters[1L]/(component_clusters[1L]-1)
  B <- matrix(0, 14, 14)
  B[1:7, 1:7] <- results[[1]]$bread; B[8:14, 8:14] <- results[[2]]$bread
  Bn <- matrix(0, 14, 14)
  Bn[1:7, 1:7] <- results[[1]]$bread_natural_delta
  Bn[8:14, 8:14] <- results[[2]]$bread_natural_delta
  Bt <- matrix(0, 14, 14)
  Bt[1:7, 1:7] <- results[[1]]$bread_native_delta
  Bt[8:14, 8:14] <- results[[2]]$bread_native_delta
  correction <- s$n_clusters/(s$n_clusters-1)
  M <- crossprod(S); Bs <- s$modelbased; Vs <- s$vcov
  stopifnot(min(eigen(Bs, symmetric = TRUE)$values) > 0)
  # V = correction * B M B', with symmetric bread B. Undo both bread factors.
  Ms <- solve(Bs, Vs) %*% solve(Bs)/correction
  Vind <- B %*% M %*% B*correction
  Vhybrid <- Bs %*% M %*% Bs*correction
  Vreplace <- B %*% Ms %*% B*correction
  normalization <- r_correction/correction
  Bc <- matrix(0, 14, 14)
  Bc[1:7, 1:7] <- components[[outcomes[1]]]$vcov
  Bc[8:14, 8:14] <- components[[outcomes[2]]]$vcov
  Vcomponent <- Bc %*% Ms %*% Bc*r_correction
  metrics <- function(A, reference) c(max_abs_covariance = max(abs(A-reference)),
    max_relative_diagonal_sd = max(abs(sqrt(diag(A)/diag(reference))-1)),
    max_scaled_covariance = max(abs(A-reference)/sqrt(outer(diag(reference), diag(reference)))))
  comparison <- rbind(
    bread_stata_vs_independent = metrics(Bs, B),
    meat_recovered_stata_vs_independent = metrics(Ms, M),
    sandwich_independent_vs_R = metrics(Vind*normalization, r$vcov),
    native_coordinate_sandwich_vs_R = metrics(Bt %*% M %*% Bt*r_correction, r$vcov),
    raw_stata_vs_R = metrics(Vs, r$vcov),
    normalized_stata_vs_R = metrics(Vs*normalization, r$vcov),
    component_bread_stata_meat_vs_R = metrics(Vcomponent, r$vcov),
    stata_bread_independent_meat_vs_stata = metrics(Vhybrid, Vs),
    independent_bread_stata_meat_vs_independent = metrics(Vreplace, Vind),
    raw_stata_vs_independent = metrics(Vs, Vind),
    coordinate_bread_sensitivity = metrics(Bn, B),
    coordinate_sandwich_sensitivity = metrics(Bn %*% M %*% Bn*correction, Vind))
  cat("CLUSTER_NORMALIZATION=", normalization, "\n")
  print(comparison, digits = 10)
  ll_difference <- sum(vapply(results, `[[`, numeric(1), "logLik"))-s$logLik
  cat("INDEPENDENT_MINUS_STATA_LOGLIK=", format(ll_difference, digits = 16), "\n")
  cat("STATA_MODELBASED_MAX_CROSS_BLOCK=", max(abs(Bs[1:7, 8:14])), "\n")
  # Equality tolerances are applied by the caller at the matching parameter
  # point. Raw GSEM/component-OIM and different-point comparisons are diagnostic.
  stopifnot(abs(ll_difference) < 1e-5,
    max(abs(Bs %*% Ms %*% Bs*correction-Vs)) < 1e-12)
  list(components = results, comparison = comparison, independent_bread = B,
    natural_delta_bread = Bn, independent_meat = M, stata_bread = Bs,
    native_delta_bread = Bt, native_coordinate_vcov = Bt %*% M %*% Bt*r_correction,
    recovered_stata_meat = Ms, independent_vcov = Vind*normalization,
    stata_bread_independent_meat_vcov = Vhybrid,
    independent_bread_stata_meat_vcov = Vreplace*normalization,
    component_bread_stata_meat_vcov = Vcomponent,
    metadata = list(evaluate_at = evaluate_at, outcomes = outcomes, cluster = cluster, n_clusters = s$n_clusters,
      component_clusters = component_clusters, stata_correction = correction,
      r_correction = r_correction, normalization = normalization,
      logLik_difference = ll_difference))
}
