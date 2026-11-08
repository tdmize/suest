# Gompertz and generalized gamma survival models from flexsurv::flexsurvreg()
# (Stata's streg with distribution(gompertz) or distribution(ggamma)).
# Parameters are flexsurv's coefficients: the baseline value of each
# distribution parameter on flexsurv's scale (log rate and shape for Gompertz;
# mu, log sigma, and Q for the generalized gamma), then covariate effects.

.suest_flexsurv_adapter <- function(model) {
  dist <- model$dlist$name
  if (!dist %in% c("gompertz", "gengamma"))
    stop(paste0("flexsurv::flexsurvreg() models are supported for the Gompertz ",
      "(dist = \"gompertz\") and generalized gamma (dist = \"gengamma\") ",
      "distributions; fit the other parametric survival models with ",
      "survival::survreg()."), call. = FALSE)
  if (!identical(model$opt$convergence, 0L) && !identical(model$opt$convergence, 0))
    stop("The flexsurv::flexsurvreg() model must converge.", call. = FALSE)
  if (length(model$fixedpars))
    stop("flexsurv::flexsurvreg() models with fixed parameters are not supported.",
      call. = FALSE)
  for (argument in c("weights", "bhazard", "rtrunc", "subset"))
    if (!is.null(model$call[[argument]]))
      stop(sprintf("flexsurv::flexsurvreg() models with %s are not supported.",
        argument), call. = FALSE)
  Y <- model$data$Y
  if (!identical(attr(Y, "type"), "right") && !identical(attr(Y, "type"), "counting"))
    stop("flexsurv::flexsurvreg() models must have right-censored outcomes.", call. = FALSE)
  if (any(is.finite(Y[, "time2"]) & Y[, "status"] != 1) ||
      any(!Y[, "status"] %in% c(0, 1)))
    stop("flexsurv::flexsurvreg() models must have right-censored outcomes.", call. = FALSE)
  if (!is.null(stats::model.offset(model$data$m)))
    stop("Offsets are not supported for flexsurv::flexsurvreg() models.", call. = FALSE)
  if (!requireNamespace("flexsurv", quietly = TRUE))
    stop("Package 'flexsurv' is required for flexsurv::flexsurvreg() models.", call. = FALSE)
  if (dist == "gengamma" && !requireNamespace("numDeriv", quietly = TRUE))
    stop("Package 'numDeriv' is required for generalized gamma survival models.",
      call. = FALSE)
  list(engine = "flexsurv::flexsurvreg", type = "flexsurv")
}

# The environment of the model formula (flexsurvreg has no formula method).
.suest_flexsurv_environment <- function(model) {
  env <- environment(attr(model$data$m, "terms"))
  if (is.environment(env)) env else parent.frame()
}

.suest_flexsurv_parameters <- function(model) {
  if (!is.null(model$suest_parameters)) return(model$suest_parameters)
  stats::setNames(as.numeric(model$coefficients), names(model$coefficients))
}

# Positions of each distribution parameter's coefficients: its baseline value,
# then its covariate effects, matching the columns of its design matrix.
.suest_flexsurv_index <- function(model) {
  pars <- model$dlist$pars
  stats::setNames(lapply(seq_along(pars), function(j)
    c(model$basepars[j], model$covpars[model$mx[[pars[j]]]])), pars)
}

.suest_flexsurv_design <- function(model, newdata = NULL) {
  pars <- model$dlist$pars
  if (is.null(newdata)) {
    X <- lapply(pars, function(p) model$data$mml[[p]])
    names(X) <- pars
  } else {
    X <- lapply(pars, function(p) {
      form <- model$all.formulae[[p]]
      if (is.null(form)) return(matrix(1, nrow(newdata), 1, dimnames = list(NULL, "(Intercept)")))
      terms <- stats::delete.response(stats::terms(form))
      frame <- stats::model.frame(terms, newdata, na.action = stats::na.pass,
        xlev = model$covdata$xlev)
      stats::model.matrix(terms, frame, contrasts.arg = model$all.contrasts[[p]])
    })
    names(X) <- pars
  }
  index <- .suest_flexsurv_index(model)
  n <- if (is.null(newdata)) nrow(model$data$m) else nrow(newdata)
  for (p in pars) {
    if (is.null(X[[p]])) X[[p]] <- matrix(1, n, 1, dimnames = list(NULL, "(Intercept)"))
    if (ncol(X[[p]]) != length(index[[p]]))
      stop("Unable to construct the flexsurv design matrices.", call. = FALSE)
  }
  X
}

# Distribution parameters for every observation.
.suest_flexsurv_values <- function(model, theta, X) {
  index <- .suest_flexsurv_index(model)
  pars <- model$dlist$pars
  stats::setNames(lapply(seq_along(pars), function(j) {
    lp <- as.numeric(X[[pars[j]]] %*% theta[index[[pars[j]]]])
    model$dlist$inv.transforms[[j]](lp)
  }), pars)
}

# Log density of the generalized gamma (Prentice's parameterization, as in
# flexsurv::dgengamma and Stata's ggamma). flexsurv's density loses accuracy
# when Q is near 0 (an error of 1e-6 at Q = 1e-4 and 0.3 at Q = 1e-7), which
# makes numerical scores noisy when Q depends on covariates and crosses 0. For
# |Q| < 0.1 this uses Stirling's series for lgamma(Q^-2) and computes
# (exp(a) - 1 - a)/a^2 by its series, which are exact to rounding there and
# give the lognormal density at Q = 0.
.suest_gengamma_log_density <- function(t, mu, sigma, Q) {
  n <- max(length(t), length(mu), length(sigma), length(Q))
  t <- rep_len(t, n); mu <- rep_len(mu, n); sigma <- rep_len(sigma, n); Q <- rep_len(Q, n)
  w <- (log(t) - mu)/sigma
  out <- numeric(n)
  small <- abs(Q) < 0.1
  if (any(!small)) {
    q <- Q[!small]; g <- q^-2; a <- q*w[!small]
    out[!small] <- log(abs(q)) - log(sigma[!small]) - log(t[!small]) + g*log(g) -
      lgamma(g) + g*(a - exp(a))
  }
  if (any(small)) {
    q <- Q[small]; a <- q*w[small]
    h <- numeric(length(a))
    near <- abs(a) < 0.5
    term <- rep(1/2, sum(near)); s <- numeric(sum(near))
    for (k in 0:24) { s <- s + term; term <- term*a[near]/(k + 3) }
    h[near] <- s
    h[!near] <- (expm1(a[!near]) - a[!near])/a[!near]^2
    out[small] <- -0.5*log(2*pi) - log(sigma[small]) - log(t[small]) - q^2/12 + q^6/360 -
      q^10/1260 + q^14/1680 - q^18/1188 - w[small]^2*h
  }
  out
}

.suest_flexsurv_loglik <- function(model, theta, X, Y) {
  v <- .suest_flexsurv_values(model, theta, X)
  time <- Y[, "stop"]
  start <- Y[, "start"]
  event <- Y[, "status"] == 1
  if (model$dlist$name == "gompertz") {
    log_f <- flexsurv::dgompertz(time, v$shape, v$rate, log = TRUE)
    log_s <- flexsurv::pgompertz(time, v$shape, v$rate, lower.tail = FALSE, log.p = TRUE)
    log_s0 <- flexsurv::pgompertz(start, v$shape, v$rate, lower.tail = FALSE, log.p = TRUE)
  } else {
    log_f <- .suest_gengamma_log_density(time, v$mu, v$sigma, v$Q)
    log_s <- flexsurv::pgengamma(time, v$mu, v$sigma, v$Q, lower.tail = FALSE, log.p = TRUE)
    log_s0 <- ifelse(start > 0, flexsurv::pgengamma(pmax(start, .Machine$double.xmin),
      v$mu, v$sigma, v$Q, lower.tail = FALSE, log.p = TRUE), 0)
  }
  ifelse(event, log_f, log_s) - ifelse(start > 0, log_s0, 0)
}

# Gompertz scores in (linear predictor of log rate, shape): exact.
.suest_gompertz_scores <- function(model, theta, X, Y) {
  v <- .suest_flexsurv_values(model, theta, X)
  gamma <- v$shape
  rate <- v$rate
  g <- function(t) {
    gt <- gamma*t
    ifelse(abs(gt) < 1e-3, t + gamma*t^2/2 + gamma^2*t^3/6 + gamma^3*t^4/24,
      expm1(gt)/ifelse(gamma == 0, 1, gamma))
  }
  dg <- function(t) {
    gt <- gamma*t
    ifelse(abs(gt) < 1e-3, t^2/2 + gamma*t^3/3 + gamma^2*t^4/8 + gamma^3*t^5/30,
      (gt*exp(gt) - expm1(gt))/ifelse(gamma == 0, 1, gamma^2))
  }
  time <- Y[, "stop"]
  start <- Y[, "start"]
  event <- as.numeric(Y[, "status"] == 1)
  d_rate <- event - rate*g(time) + rate*g(start)
  d_shape <- event*time - rate*dg(time) + rate*dg(start)
  index <- .suest_flexsurv_index(model)
  U <- matrix(0, nrow(Y), length(theta), dimnames = list(NULL, names(theta)))
  U[, index$rate] <- X$rate*d_rate
  U[, index$shape] <- X$shape*d_shape
  U
}

.suest_flexsurv_components <- function(model) {
  theta <- .suest_flexsurv_parameters(model)
  X <- .suest_flexsurv_design(model)
  Y <- model$data$Y
  loglik <- function(b) .suest_flexsurv_loglik(model, b, X, Y)
  check <- sum(loglik(theta))
  if (!is.finite(check) ||
      abs(check - as.numeric(model$loglik)) > 1e-6*max(1, abs(as.numeric(model$loglik))))
    stop("Could not reproduce the flexsurv::flexsurvreg() log likelihood.", call. = FALSE)
  if (model$dlist$name == "gompertz") {
    score <- function(b) .suest_gompertz_scores(model, b, X, Y)
    U <- score(theta)
    information <- .suest_score_information(theta, function(b) colSums(score(b)))
  } else {
    U <- numDeriv::jacobian(loglik, theta, method = "Richardson")
    information <- -numDeriv::hessian(function(b) sum(loglik(b)), theta)
    information <- (information + t(information))/2
  }
  colnames(U) <- names(theta)
  rownames(U) <- rownames(model$data$m)
  B <- nrow(U)*.suest_information_inverse(information)
  dimnames(B) <- list(names(theta), names(theta))
  if (any(!is.finite(U)) || any(!is.finite(B)))
    stop("Unable to construct finite flexsurv scores and bread.", call. = FALSE)
  list(score = U, bread = B, parameters = theta)
}

# Response: median survival time. Link: linear predictor of the location
# parameter (log rate for Gompertz, mu for the generalized gamma).
.suest_flexsurv_predict <- function(model, newdata, type) {
  theta <- .suest_flexsurv_parameters(model)
  X <- .suest_flexsurv_design(model, newdata)
  location <- model$dlist$location
  if (identical(type, "link")) {
    index <- .suest_flexsurv_index(model)
    return(as.numeric(X[[location]] %*% theta[index[[location]]]))
  }
  v <- .suest_flexsurv_values(model, theta, X)
  if (model$dlist$name == "gompertz")
    flexsurv::qgompertz(0.5, v$shape, v$rate) else
    flexsurv::qgengamma(0.5, v$mu, v$sigma, v$Q)
}
