# Generalized ordered logit (Stata's gologit2) from VGAM::vglm() with the
# cumulative() family. With reverse = TRUE, P(Y > j) = F(eta_j), the gologit2
# parameterization. Constraint matrices give the proportional (parallel) and
# partial proportional odds forms. suest computes the scores analytically and
# the observed information by differentiating them.

.suest_wrap_vglm <- function(model) {
  if (!methods::is(model, "vglm")) return(model)
  if (!requireNamespace("VGAM", quietly = TRUE))
    stop("Package 'VGAM' is required for VGAM::vglm() models.", call. = FALSE)
  formula <- stats::formula(model)
  call <- methods::slot(model, "call")
  data <- try(eval(call$data, envir = environment(formula)), silent = TRUE)
  if (inherits(data, "try-error") || !is.data.frame(data))
    stop("Could not recover the original data for the VGAM::vglm() model; keep its data object available.",
      call. = FALSE)
  terms <- stats::terms(model)
  frame <- stats::model.frame(terms, data = data, na.action = stats::na.omit)
  vlm_X <- stats::model.matrix(model, type = "lm")
  if (nrow(frame) != nrow(vlm_X))
    stop("Could not align the VGAM::vglm() estimation sample to its data.", call. = FALSE)
  factors <- vapply(frame, is.factor, logical(1)) &
    names(frame) %in% all.vars(stats::delete.response(terms))
  structure(list(
    fit = model,
    call = call,
    formula = formula,
    terms = terms,
    frame = frame,
    xlevels = lapply(frame[factors], levels),
    contrasts = attr(vlm_X, "contrasts")
  ), class = "suest_vglm")
}

.suest_vglm_info <- function(wrapper) {
  fit <- wrapper$fit
  family <- methods::slot(fit, "family")
  if (!"cumulative" %in% family@vfamily)
    stop("Only VGAM::vglm() models with the cumulative() family (generalized ordered models) are supported.",
      call. = FALSE)
  links <- unique(unname(fit@misc$link))
  link <- c(logitlink = "logit", probitlink = "probit", clogloglink = "cloglog")[links]
  if (length(links) != 1L || is.na(link))
    stop("Generalized ordered models must use a logit, probit, or cloglog link.", call. = FALSE)
  y <- stats::model.response(wrapper$frame)
  levels <- if (is.factor(y)) levels(y) else as.character(sort(unique(y)))
  list(link = unname(link), reverse = isTRUE(fit@misc$reverse),
    levels = levels, M = length(levels) - 1L)
}

.suest_vglm_adapter <- function(wrapper) {
  info <- .suest_vglm_info(wrapper)
  fit <- wrapper$fit
  w <- methods::slot(fit, "prior.weights")
  if (!is.null(wrapper$call$weights) || (length(w) && any(w != 1)))
    stop("Weighted VGAM::vglm() models are not supported.", call. = FALSE)
  if (!is.null(wrapper$call$offset) || length(attr(wrapper$terms, "offset")))
    stop("Offsets are not supported for VGAM::vglm() models.", call. = FALSE)
  y <- stats::model.response(wrapper$frame)
  if (!is.factor(y) || is.matrix(y) || info$M < 2L)
    stop("Generalized ordered models need a factor response with at least three categories, one row per observation.",
      call. = FALSE)
  list(engine = "VGAM::vglm", type = "gologit")
}

.suest_vglm_parameters <- function(wrapper) {
  if (!is.null(wrapper$suest_parameters)) return(wrapper$suest_parameters)
  methods::slot(wrapper$fit, "coefficients")
}

# Linear predictors (n x M) from an ordinary model matrix and the constraints.
.suest_vglm_eta <- function(wrapper, parameters, X) {
  constraints <- methods::slot(wrapper$fit, "constraints")
  assign <- attr(X, "assign")
  labels <- c("(Intercept)", attr(wrapper$terms, "term.labels"))
  M <- .suest_vglm_info(wrapper)$M
  eta <- matrix(0, nrow(X), M)
  position <- 0L
  for (term in names(constraints)) {
    H <- constraints[[term]]
    columns <- which(labels[assign + 1L] == term)
    for (column in columns) {
      index <- position + seq_len(ncol(H))
      eta <- eta + outer(X[, column], drop(H %*% parameters[index]))
      position <- position + ncol(H)
    }
  }
  if (position != length(parameters))
    stop("Could not map the VGAM::vglm() coefficients to its constraints.", call. = FALSE)
  eta
}

.suest_vglm_cdf <- function(link) switch(link,
  logit = list(F = stats::plogis, f = stats::dlogis),
  probit = list(F = stats::pnorm, f = stats::dnorm),
  cloglog = list(F = function(x) -expm1(-exp(x)), f = function(x) exp(x - exp(x))))

# Category probabilities (n x K) from the linear predictors.
.suest_vglm_probabilities <- function(eta, info) {
  cdf <- .suest_vglm_cdf(info$link)
  upper <- if (info$reverse) cdf$F(eta) else 1 - cdf$F(eta)
  above <- cbind(1, upper, 0)
  probabilities <- above[, -ncol(above), drop = FALSE] - above[, -1L, drop = FALSE]
  colnames(probabilities) <- info$levels
  probabilities
}

.suest_vglm_components <- function(wrapper) {
  info <- .suest_vglm_info(wrapper)
  parameters <- .suest_vglm_parameters(wrapper)
  X <- stats::model.matrix(wrapper$terms, wrapper$frame, contrasts.arg = wrapper$contrasts)
  y <- as.integer(stats::model.response(wrapper$frame))
  n <- nrow(X)
  cdf <- .suest_vglm_cdf(info$link)
  sign <- if (info$reverse) 1 else -1

  # d log P(Y = y_i)/d eta_ij for the two cutpoints bordering y_i.
  score_eta <- function(theta) {
    eta <- .suest_vglm_eta(wrapper, theta, X)
    probabilities <- .suest_vglm_probabilities(eta, info)
    p <- probabilities[cbind(seq_len(n), y)]
    density <- sign*cdf$f(eta)
    out <- matrix(0, n, info$M)
    lower <- y - 1L
    has_lower <- lower >= 1L
    out[cbind(which(has_lower), lower[has_lower])] <- density[cbind(which(has_lower), lower[has_lower])]/p[has_lower]
    has_upper <- y <= info$M
    out[cbind(which(has_upper), y[has_upper])] <- -density[cbind(which(has_upper), y[has_upper])]/p[has_upper]
    list(d = out, loglik = sum(log(p)))
  }
  # Chain rule through the constraints: each parameter's derivative of eta.
  design <- lapply(seq_along(parameters), function(k) {
    unit <- replace(numeric(length(parameters)), k, 1)
    .suest_vglm_eta(wrapper, unit, X)
  })
  scores <- function(theta) {
    d <- score_eta(theta)$d
    vapply(design, function(D) rowSums(d*D), numeric(n))
  }

  total <- score_eta(parameters)$loglik
  if (!isTRUE(all.equal(total, methods::slot(wrapper$fit, "criterion")$loglikelihood,
      tolerance = 1e-8)))
    stop("Could not reproduce the VGAM::vglm() log likelihood for its scores.", call. = FALSE)
  score <- matrix(scores(parameters), n, dimnames = list(rownames(wrapper$frame), names(parameters)))
  step <- 1e-5*pmax(1, abs(parameters))
  information <- -vapply(seq_along(parameters), function(k) {
    up <- down <- parameters
    up[k] <- up[k] + step[k]
    down[k] <- down[k] - step[k]
    (colSums(scores(up)) - colSums(scores(down)))/(2*step[k])
  }, numeric(length(parameters)))
  information <- (information + t(information))/2
  B <- n*.suest_information_inverse(information)
  dimnames(B) <- list(names(parameters), names(parameters))
  list(score = score, bread = B, parameters = parameters)
}

.suest_vglm_predict <- function(wrapper, newdata) {
  info <- .suest_vglm_info(wrapper)
  terms <- stats::delete.response(wrapper$terms)
  frame <- stats::model.frame(terms, newdata, na.action = stats::na.pass,
    xlev = wrapper$xlevels)
  X <- stats::model.matrix(terms, frame, contrasts.arg = wrapper$contrasts)
  .suest_vglm_probabilities(
    .suest_vglm_eta(wrapper, .suest_vglm_parameters(wrapper), X), info)
}

.suest_vglm_frame <- function(wrapper) {
  frame <- as.data.frame(wrapper$frame)
  original <- try(eval(wrapper$call$data, envir = environment(wrapper$formula)),
    silent = TRUE)
  if (!inherits(original, "try-error") && is.data.frame(original)) {
    index <- match(rownames(frame), rownames(original))
    if (!anyNA(index))
      for (variable in setdiff(names(original), names(frame)))
        frame[[variable]] <- original[[variable]][index]
  }
  frame
}
