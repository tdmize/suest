# Release validation gate 1. Run from the suest repository root:
# source("suest_release_gate1.R")  # after extracting the review ZIP here
# Or: source("tools/release-validation/suest_release_gate1.R")
# This script loads local package source. It does not install or edit anything.
# Prepared and executed in R 4.5.2 on 2026-09-28.

run_suest_release_gate1 <- function(root = getwd(), bootstrap_reps = 999L) {
  root <- normalizePath(root, mustWork = TRUE)
  log_path <- file.path(root, "suest_release_gate1.log")
  con <- file(log_path, "wt")
  old_message <- sink.number(type = "message")
  sink(con, split = TRUE)
  sink(con, type = "message")
  on.exit({
    sink(if (old_message == 2L) NULL else getConnection(old_message), type = "message")
    sink()
    close(con)
  }, add = TRUE)
  old_options <- options(warn = 1, digits = 15, width = 180)
  on.exit(options(old_options), add = TRUE)
  cat("SUEST RELEASE GATE 1 - diagnostic run, not full certification\n")
  cat("Started:", format(Sys.time(), tz = "UTC"), "UTC\nRoot:", root, "\n")
  setup_ok <- tryCatch({
    if (!file.exists(file.path(root, "R", "suest.R")))
      stop("Run this script from the suest source repository root.")
    required <- c("pkgload", "sandwich", "marginaleffects", "nnet")
    missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) stop("Install missing packages first: ", paste(missing, collapse = ", "))
    if (length(bootstrap_reps) != 1L || !is.finite(bootstrap_reps) ||
        bootstrap_reps < 199L || bootstrap_reps != as.integer(bootstrap_reps))
      stop("bootstrap_reps must be an integer of at least 199.")
    pkgload::load_all(root, quiet = TRUE)
    TRUE
  }, error = function(e) {
    cat("SETUP_ERROR:", conditionMessage(e), "\nGATE1_COMPLETED=0\n")
    FALSE
  })
  if (!setup_ok) return(invisible(list(status = "SETUP_ERROR", log = log_path)))
  cat("Version:", read.dcf(file.path(root, "DESCRIPTION"))[1, "Version"], "\n")
  if (nzchar(Sys.which("git"))) {
    cat("Git HEAD:", system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE), "\n")
    cat("Git working tree:\n")
    cat(system2("git", c("-C", shQuote(root), "status", "--short"), stdout = TRUE), sep = "\n")
  }
  cat("\nSource fingerprints:\n")
  print(tools::md5sum(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)))
  cat("Loaded source:", getNamespaceInfo(asNamespace("suest"), "path"), "\n")
  results <- list()
  run_case <- function(name, fun) {
    cat("\nCASE:", name, "\n")
    warnings <- character()
    ans <- tryCatch(withCallingHandlers(fun(), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
    }), error = function(e) list(status = "ERROR", error = conditionMessage(e)))
    if (is.null(ans$status)) ans$status <- "REVIEW"
    if (length(warnings) && ans$status %in% c("PASS", "SCREEN_CLEAR"))
      ans$status <- "REVIEW_WARNING"
    ans$warnings <- unique(warnings)
    results[[name]] <<- ans
    print(ans)
    invisible(ans)
  }
  rejection_probe <- function(model, expected_message) {
    value <- tryCatch(suest::suest(model, model, model_names = c("A", "B")), error = identity)
    if (inherits(value, "error")) return(list(
      status = if (grepl(expected_message, conditionMessage(value), ignore.case = TRUE))
        "PASS_REJECTION" else "REVIEW_REJECTION",
      message = conditionMessage(value)))
    list(status = "ACCEPTED_REQUIRES_REVIEW", finite_covariance = all(is.finite(vcov(value))),
      max_standard_error = max(sqrt(diag(vcov(value)))))
  }

  set.seed(20260928)
  n <- 1200L
  dat <- data.frame(id = seq_len(n), cluster = rep(seq_len(n/4L), each = 4L),
    x = runif(n, -1, 1), z = runif(n, -1, 1))
  shared <- rep(rnorm(n/4L), each = 4L)
  e1 <- .35*shared + rnorm(n)
  e2 <- .35*shared + .4*e1 + rnorm(n)
  eta1 <- -.2 + .65*dat$x - .3*dat$z
  eta2 <- .1 + .35*dat$x + .25*dat$z
  dat$binary1 <- as.integer(eta1 + e1 > 0)
  dat$binary2 <- as.integer(eta2 + e2 > 0)
  dat$normal1 <- eta1 + e1
  dat$normal2 <- eta2 + e2
  dat$positive1 <- exp(.2 + .25*eta1 + .45*e1)
  dat$positive2 <- exp(.1 + .25*eta2 + .45*e2)
  dat$fraction1 <- plogis(eta1 + .65*e1)
  dat$fraction2 <- plogis(eta2 + .65*e2)
  dat$count1 <- qpois(pnorm(e1/sqrt(1+.35^2)), exp(.4+.3*dat$x-.2*dat$z))
  dat$count2 <- qpois(pnorm(e2/sqrt(1+.4^2+(.35*1.4)^2)), exp(.3+.2*dat$x+.1*dat$z))

  run_case("Converged ordinary GLM control", function() {
    good <- glm(binary1 ~ x + z, dat, family = binomial(), control = glm.control(epsilon = 1e-10))
    stopifnot(isTRUE(good$converged))
    fit <- suest::suest(good, good)
    list(status = if (all(is.finite(vcov(fit)))) "PASS" else "FAIL")
  })
  run_case("Actual nonconverged ordinary GLM", function() {
    bad <- glm(binary1 ~ x + z, dat, family = binomial(), control = glm.control(maxit = 1))
    if (isTRUE(bad$converged)) stop("Probe did not create a nonconverged fit.")
    c(list(native_converged = bad$converged), rejection_probe(bad, "converg"))
  })
  run_case("Multinomial decay=0 control", function() {
    dat$category <- factor(ifelse(eta1 + e1 < -.4, "A", ifelse(eta2 + e2 > .6, "C", "B")))
    fit <- nnet::multinom(category ~ x + z, dat, decay = 0, trace = FALSE, Hess = TRUE, maxit = 500)
    stopifnot(fit$convergence == 0)
    combined <- suest::suest(fit, fit)
    list(status = if (all(is.finite(vcov(combined)))) "PASS" else "FAIL")
  })
  run_case("Penalized multinomial decay=2", function() {
    dat$category <- factor(ifelse(eta1 + e1 < -.4, "A", ifelse(eta2 + e2 > .6, "C", "B")))
    fit <- nnet::multinom(category ~ x + z, dat, decay = 2, trace = FALSE, Hess = TRUE, maxit = 500)
    stopifnot(fit$convergence == 0)
    X <- model.matrix(~x + z, dat)
    P <- fitted(fit)
    Y <- sapply(fit$lev, function(level) as.numeric(dat$category == level))
    ordinary_score <- crossprod(X, Y[, -1, drop = FALSE] - P[, -1, drop = FALSE])
    c(list(decay = fit$decay, max_unpenalized_score_sum = max(abs(ordinary_score))),
      rejection_probe(fit, "penal|decay|regulariz"))
  })

  # Same-family stacked fits estimate separate equation coefficients. The
  # cluster sandwich supplies all cross-equation blocks independently of suest.
  # Every case explicitly clusters, aligning G/(G-1) without extra HC1 factors.
  loglog <- structure(list(linkfun = function(mu) -log(-log(mu)),
    linkinv = function(eta) exp(-exp(-eta)),
    mu.eta = function(eta) exp(-eta-exp(-eta)),
    valideta = function(eta) TRUE, name = "loglog"), class = "link-glm")
  specs <- list(
    gaussian = list(outcome = "normal", family = gaussian()),
    logit = list(outcome = "binary", family = binomial()),
    probit = list(outcome = "binary", family = binomial("probit")),
    cloglog = list(outcome = "binary", family = binomial("cloglog")),
    poisson = list(outcome = "count", family = poisson()),
    gamma = list(outcome = "positive", family = Gamma("log")),
    fractional_logit = list(outcome = "fraction", family = quasibinomial()),
    fractional_probit = list(outcome = "fraction", family = quasibinomial("probit")),
    fractional_cloglog = list(outcome = "fraction", family = quasibinomial("cloglog")),
    fractional_loglog = list(outcome = "fraction", family = quasibinomial(loglog)))
  patterns <- list(identical = list(rep(TRUE, n), rep(TRUE, n)),
    partial = list(dat$id %% 4 != 0, dat$id %% 4 != 1),
    disjoint_shared_clusters = list(dat$id %% 2 == 0, dat$id %% 2 == 1),
    disjoint_clusters = list(dat$cluster <= 150, dat$cluster > 150))
  nd <- dat[, c("x", "z")]
  nd0 <- nd1 <- nd
  nd0$x <- -.5; nd1$x <- .5
  X <- model.matrix(~x + z, nd)
  X0 <- model.matrix(~x + z, nd0)
  X1 <- model.matrix(~x + z, nd1)
  estimand <- function(beta, family, kind) {
    if (kind == "prediction") return(mean(family$linkinv(drop(X %*% beta))))
    mean(family$linkinv(drop(X1 %*% beta)) - family$linkinv(drop(X0 %*% beta)))
  }
  gradient <- function(beta, family, kind) {
    if (kind == "prediction") return(colMeans(X * family$mu.eta(drop(X %*% beta))))
    colMeans(X1 * family$mu.eta(drop(X1 %*% beta)) -
      X0 * family$mu.eta(drop(X0 %*% beta)))
  }
  fit_glm <- function(d, family) {
    m <- glm(y ~ x + z, data = d, family = family,
      control = glm.control(epsilon = 1e-11, maxit = 100))
    if (!isTRUE(m$converged) || any(!is.finite(coef(m)))) stop("Reference fit failed to converge.")
    m
  }
  bootstrap_case <- NULL
  for (name in names(specs)) for (pattern in names(patterns)) {
    run_case(paste("Stacked", name, pattern), function() {
      spec <- specs[[name]]
      a <- dat[patterns[[pattern]][[1]], ]; b <- dat[patterns[[pattern]][[2]], ]
      a$y <- a[[paste0(spec$outcome, "1")]]; b$y <- b[[paste0(spec$outcome, "2")]]
      ma <- fit_glm(a, spec$family); mb <- fit_glm(b, spec$family)
      joint <- suest::suest(ma, mb, model_names = c("A", "B"),
        observation_id = list(a$id, b$id), cluster = list(a$cluster, b$cluster))
      stack <- rbind(data.frame(a[, c("y", "cluster", "x", "z")], equation = "A"),
        data.frame(b[, c("y", "cluster", "x", "z")], equation = "B"))
      stack$a0 <- as.numeric(stack$equation == "A"); stack$b0 <- 1-stack$a0
      stack$ax <- stack$a0*stack$x; stack$az <- stack$a0*stack$z
      stack$bx <- stack$b0*stack$x; stack$bz <- stack$b0*stack$z
      stacked <- glm(y ~ 0 + a0 + ax + az + b0 + bx + bz, stack,
        family = spec$family, control = glm.control(epsilon = 1e-11, maxit = 100))
      stopifnot(isTRUE(stacked$converged), length(coef(joint)) == 6L,
        identical(names(coef(stacked)), c("a0", "ax", "az", "b0", "bx", "bz")))
      V <- unname(sandwich::vcovCL(stacked, cluster = stack$cluster, type = "HC0", cadjust = TRUE))
      information_convention <- "native GLM expected information"
      convention_gap <- 0
      if (name == "probit") {
        # Independent Bernoulli log likelihood: log Phi(s*eta), s = 2*y-1.
        # Its score is s*lambda and negative second derivative lambda*(t+lambda).
        # Ordinary binary probit in suest uses this observed-information bread.
        V_expected <- V
        Z <- model.matrix(stacked)
        sign_y <- 2*stack$y-1
        signed_eta <- sign_y*drop(Z %*% coef(stacked))
        mills <- exp(dnorm(signed_eta, log = TRUE)-pnorm(signed_eta, log.p = TRUE))
        I <- crossprod(Z, Z * (mills*(signed_eta+mills)))
        inverse <- solve(I)
        U <- rowsum(Z * (sign_y*mills), stack$cluster, reorder = FALSE)
        G <- nrow(U)
        V <- inverse %*% crossprod(U) %*% t(inverse) * G/(G-1)
        convention_gap <- max(abs(V-V_expected))
        information_convention <- "independent observed Bernoulli-probit information"
      }
      W <- unname(vcov(joint))
      stopifnot(all(is.finite(V)), all(is.finite(W)), all(diag(V) > 0))
      coefficient_gap <- max(abs(unname(coef(joint)) - unname(coef(stacked))))
      covariance_gap <- max(abs(W - V))
      standardized_covariance_gap <- max(abs(W - V)/sqrt(outer(diag(V), diag(V))))
      out <- list(coefficient_gap = coefficient_gap, covariance_gap = covariance_gap,
        standardized_covariance_gap = standardized_covariance_gap,
        reference_crossblock_max = max(abs(V[1:3, 4:6])),
        information_convention = information_convention,
        expected_vs_observed_covariance_gap = convention_gap)
      passed <- coefficient_gap < 2e-6 && standardized_covariance_gap < 1e-5
      for (kind in c("prediction", "comparison")) {
        tab <- if (kind == "prediction") marginaleffects::avg_predictions(joint, newdata = nd,
          numderiv = list("fdcenter", eps = 1e-5)) else marginaleffects::avg_comparisons(joint,
          variables = list(x = c(-.5, .5)), newdata = nd,
          numderiv = list("fdcenter", eps = 1e-5))
        stopifnot(nrow(tab) == 2L, setequal(as.character(tab$group), c("A", "B")))
        signs <- ifelse(as.character(tab$group) == "A", 1, -1)
        contrast <- marginaleffects::hypotheses(tab, hypothesis = matrix(signs, ncol = 1L))
        ba <- coef(stacked)[1:3]; bb <- coef(stacked)[4:6]
        delta <- estimand(ba, spec$family, kind) - estimand(bb, spec$family, kind)
        g <- c(gradient(ba, spec$family, kind), -gradient(bb, spec$family, kind))
        se <- sqrt(drop(g %*% V %*% g))
        estimate_gap <- abs(as.numeric(contrast$estimate) - delta)
        relative_se_gap <- abs(as.numeric(contrast$std.error)/se - 1)
        out[[kind]] <- c(reference_estimate = delta, reference_se = se,
          estimate_gap = estimate_gap, relative_se_gap = relative_se_gap)
        passed <- passed && estimate_gap < 2e-7 && relative_se_gap < 2e-4
        if (name == "logit" && pattern == "identical" && kind == "comparison")
          bootstrap_case <<- list(a = a, b = b, family = spec$family,
            analytic_se = as.numeric(contrast$std.error), estimate = delta)
      }
      out$status <- if (passed) "PASS" else "FAIL"
      out
    })
  }

  # A separate resampling check: both outcomes use the SAME resampled clusters.
  # Evaluation covariates stay fixed, matching the package's conditional target.
  run_case("Paired cluster bootstrap: logit finite-change difference", function() {
    if (is.null(bootstrap_case)) stop("The deterministic logit reference did not complete.")
    bc <- bootstrap_case
    groups <- split(seq_len(nrow(bc$a)), bc$a$cluster)
    stopifnot(identical(bc$a$id, bc$b$id))
    set.seed(20260929)
    draws <- rep(NA_real_, bootstrap_reps)
    failed <- character()
    for (r in seq_len(bootstrap_reps)) {
      rows <- unlist(groups[sample.int(length(groups), length(groups), replace = TRUE)], use.names = FALSE)
      draws[r] <- tryCatch({
        fa <- fit_glm(bc$a[rows, ], bc$family); fb <- fit_glm(bc$b[rows, ], bc$family)
        estimand(coef(fa), bc$family, "comparison") - estimand(coef(fb), bc$family, "comparison")
      }, error = function(e) { failed <<- c(failed, conditionMessage(e)); NA_real_ })
      if (r %% 100L == 0L) cat("Bootstrap replicate", r, "of", bootstrap_reps, "\n")
    }
    successful <- sum(is.finite(draws))
    if (successful < 2L) stop("Fewer than two usable bootstrap draws.")
    bootstrap_se <- sd(draws[is.finite(draws)])
    ratio <- bootstrap_se/bc$analytic_se
    # Approximate MCSE assumes roughly normal bootstrap contrasts. A 10% band
    # is a diagnostic screen, not an equality tolerance or proof of coverage.
    approximate_relative_mcse <- 1/sqrt(2*(successful-1))
    limit <- max(.10, 4*approximate_relative_mcse)
    list(status = if (length(failed) || abs(log(ratio)) > limit) "REVIEW" else "SCREEN_CLEAR",
      requested = bootstrap_reps, successful = successful,
      failed_messages = unique(failed), analytic_se = bc$analytic_se,
      bootstrap_se = bootstrap_se, se_ratio = ratio,
      approximate_relative_mcse = approximate_relative_mcse,
      log_ratio_screen_limit = limit)
  })
  cat("\nSUMMARY\n")
  print(vapply(results, `[[`, character(1), "status"))
  cat("\nSESSION INFO\n"); print(sessionInfo())
  saveRDS(results, file.path(root, "suest_release_gate1_results.rds"))
  cat("\nGATE1_COMPLETED=1\n")
  cat("Completion means the diagnostic run finished; inspect each case status.\n")
  cat("This gate does not certify every family or resolve raw GSEM discrepancies.\n")
  invisible(results)
}

suest_release_gate1_results <- run_suest_release_gate1()
