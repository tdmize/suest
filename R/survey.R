# One-stage survey linearization. The full design supplies the sampling universe;
# each model supplies coefficient influences only for its estimation sample.
.suest_check_survey_design <- function(design) {
  if (!inherits(design, "survey.design2") ||
      !identical(sort(class(design)), sort(c("survey.design2", "survey.design"))))
    stop("Survey support requires an ordinary survey::svydesign() design, not replicate-weight or two-phase designs.", call. = FALSE)
  if (ncol(design$cluster) != 1L || ncol(design$strata) != 1L)
    stop("Survey support currently requires a one-stage design.", call. = FALSE)
  if (!is.null(design$postStrata) || isTRUE(design$pps) || isTRUE(design$fpc$pps))
    stop("Calibrated, raked, post-stratified, and PPS designs are not supported.", call. = FALSE)
  if (!is.data.frame(design$variables) || anyNA(design$cluster) || anyNA(design$strata) ||
      length(design$prob) != nrow(design$variables) ||
      any(!is.finite(design$prob) | design$prob <= 0))
    stop("The survey design must retain data and finite, strictly positive sampling weights.", call. = FALSE)
  invisible(TRUE)
}

.suest_survey_ids <- function(design, observation_id, label) {
  if (!all(observation_id %in% names(design$variables)))
    stop(sprintf("Observation-ID columns are missing from %s.", label), call. = FALSE)
  .suest_id_keys(design$variables[observation_id], nrow(design$variables), label)
}

.suest_survey_vcov <- function(influence, design) {
  survey::svyrecvar(influence, design$cluster, design$strata, design$fpc)
}

.suest_survey_type <- function(model) {
  family <- model$family$family
  link <- model$family$link
  if (identical(family, "gaussian") && identical(link, "identity")) return("survey_lm")
  if (identical(family, "quasibinomial") && link %in% c("logit", "probit", "cloglog"))
    return(link)
  if (identical(family, "quasipoisson") && identical(link, "log")) return("poisson")
  "unsupported"
}

# Direct reconstruction from the retained GLM fitting state (any GLM family), independent of
# survey's optional stored influence. Algebraically, the working residual times
# the final working weight is w_i (y_i-mu_i) mu'_i / V(mu_i), and naive.cov is
# the inverse expected information used by svyglm(). Using the retained values
# is necessary for exact native-covariance reproduction because recomputing
# them at the reported coefficients can differ at the GLM stopping tolerance.
.suest_survey_binary_influence <- function(model, X, y, w, beta, offset) {
  working_residual <- as.numeric(stats::residuals(model, type = "working"))
  working_weight <- as.numeric(model$weights)
  bread <- model$naive.cov
  if (length(working_residual) != nrow(X) || length(working_weight) != nrow(X) ||
      !is.matrix(bread) || any(dim(bread) != ncol(X)) ||
      any(!is.finite(c(working_residual, working_weight, bread))) ||
      any(working_weight <= 0))
    stop("Survey binary-model linearization is numerically undefined at the fitted values.", call. = FALSE)
  (X * (working_residual * working_weight)) %*% bread
}

# Survey fits suest combines under one design.
.suest_survey_classes <- c("svyglm", "svyolr", "svy_vglm", "svysurvreg", "svyivreg")

.suest_is_survey_model <- function(model) inherits(model, .suest_survey_classes)

# Each engine returns its fitted design, model frame (row names matching the
# design's variables), coefficient influences, and the influences that
# reproduce its own vcov(). svyolr() uses a numerical Hessian and svy_vglm()
# VGAM's expected information; suest uses exact observed information, as Stata
# does, after checking the scores against the fit's own covariance.
.suest_survey_model <- function(model, env) {
  if (inherits(model, "svyglm")) return(.suest_survey_glm(model))
  if (inherits(model, "svyolr")) return(.suest_survey_olr(model))
  if (inherits(model, "svy_vglm")) return(.suest_survey_vglm(model))
  if (inherits(model, "svysurvreg")) return(.suest_survey_survreg(model))
  .suest_survey_ivreg(model, env)
}

.suest_survey_glm <- function(model) {
  type <- .suest_survey_type(model)
  if (type == "unsupported")
    stop("Survey svyglm() support is restricted to Gaussian identity, quasibinomial logit/probit/cloglog, and quasipoisson log models.", call. = FALSE)
  if (!isTRUE(model$converged) || any(!is.finite(model$coefficients)))
    stop("Survey models must converge and have finite, nonaliased coefficients.", call. = FALSE)
  if (!is.null(model$call$weights))
    stop("Additional svyglm() weights are unsupported; specify sampling weights only in svydesign().", call. = FALSE)
  if (!is.null(model$call$offset))
    stop("For survey models, specify offsets in the formula with offset(), not the offset argument.", call. = FALSE)
  mf <- model$model
  if (is.null(mf))
    stop("Refit svyglm() with model = TRUE to retain its estimation sample.", call. = FALSE)
  X <- stats::model.matrix(model)
  y <- if (type == "survey_lm") stats::model.response(mf) else model$y
  w <- as.numeric(model$prior.weights)
  beta <- stats::coef(model)
  basic_bad <- nrow(X) != nrow(mf) || ncol(X) != length(beta) || nrow(X) <= ncol(X) ||
    length(w) != nrow(X) || any(!is.finite(w) | w <= 0)
  if (type == "survey_lm") {
    if (!is.numeric(y) || is.matrix(y) || any(!is.finite(y)) || basic_bad)
      stop("Survey linear models require a full-rank numeric response, positive weights, and residual observations.", call. = FALSE)
  } else if (type == "poisson") {
    if (is.null(y) || !is.numeric(y) || is.matrix(y) || any(!is.finite(y)) ||
        any(y < 0) || basic_bad)
      stop("Survey Poisson models require a full-rank nonnegative numeric response, positive weights, and residual observations.", call. = FALSE)
  } else if (is.null(y) || !is.numeric(y) || is.matrix(y) || any(!is.finite(y)) ||
             any(!y %in% c(0, 1)) || basic_bad) {
    stop("Survey binary models require a full-rank numeric 0/1 response, positive weights, and residual observations.", call. = FALSE)
  }
  offset <- stats::model.offset(mf)
  if (is.null(offset)) offset <- rep(0, nrow(X))
  influence <- if (type == "survey_lm") {
    residual <- y - as.numeric(X %*% beta) - offset
    (X * (w * residual)) %*% solve(crossprod(X, X * w))
  } else {
    .suest_survey_binary_influence(model, X, y, w, beta, offset)
  }
  list(type = type, engine = "survey::svyglm", design = model$survey.design, frame = mf,
    weights = w, parameters = beta, influence = influence, model = model)
}

.suest_survey_olr <- function(model) {
  if (!model$method %in% c("logistic", "probit"))
    stop("Survey ordered models must use the logistic or probit method.", call. = FALSE)
  if (!identical(model$convergence, 0L) || any(!is.finite(c(model$coefficients, model$zeta))))
    stop("The svyolr() model did not converge.", call. = FALSE)
  design <- model$survey.design
  mf <- stats::model.frame(model$terms, design$variables, na.action = stats::na.omit)
  X <- stats::model.matrix(model$terms, mf)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  if (nrow(mf) != nrow(design$variables) || !identical(colnames(X), names(model$coefficients)))
    stop("Could not align the svyolr() model with its fitted survey design.", call. = FALSE)
  w <- 1/design$prob
  parts <- .suest_polr_components(model, X = X, y = as.integer(stats::model.response(mf)),
    weights = w, offset = stats::model.offset(mf))
  model$xlevels <- stats::.getXlevels(model$terms, mf)
  model$contrasts <- attr(X, "contrasts")
  list(type = if (model$method == "logistic") "ologit" else "oprobit",
    engine = "survey::svyolr", design = design, frame = mf, weights = w,
    parameters = parts$parameters, influence = parts$score %*% parts$bread/nrow(X),
    native_influence = parts$score %*% solve(model$Hessian), model = model,
    levels = model$lev)
}

.suest_survey_vglm <- function(model) {
  if (!requireNamespace("VGAM", quietly = TRUE))
    stop("Package 'VGAM' is required for svyVGAM::svy_vglm() models.", call. = FALSE)
  fit <- model$fit
  design <- model$design
  terms <- stats::terms(fit)
  if (length(attr(terms, "offset")) || length(fit@offset) && any(fit@offset != 0))
    stop("Offsets are not supported for survey VGAM models.", call. = FALSE)
  frame <- stats::model.frame(terms, design$variables, na.action = stats::na.omit)
  vlm_X <- stats::model.matrix(fit, type = "lm")
  if (nrow(frame) != nrow(vlm_X) || nrow(frame) != nrow(design$variables))
    stop("Could not align the svy_vglm() model with its fitted survey design.", call. = FALSE)
  factors <- vapply(frame, is.factor, logical(1)) &
    names(frame) %in% all.vars(stats::delete.response(terms))
  wrapper <- structure(list(fit = fit, call = methods::slot(fit, "call"),
    formula = stats::formula(fit), terms = terms, frame = frame,
    xlevels = lapply(frame[factors], levels), contrasts = attr(vlm_X, "contrasts"),
    weights = as.numeric(methods::slot(fit, "prior.weights")), survey = TRUE),
    class = "suest_vglm")
  info <- .suest_vglm_info(wrapper)
  y <- stats::model.response(frame)
  if (info$family != "negbinomial" && (!is.factor(y) || info$M < 2L))
    stop("Survey ordered and multinomial VGAM models need a factor response with at least three categories.",
      call. = FALSE)
  parts <- .suest_vglm_components(wrapper)
  list(type = c(cumulative = "gologit", multinomial = "multinom", negbinomial = "negbin")[[info$family]],
    engine = "svyVGAM::svy_vglm", design = design, frame = frame, weights = wrapper$weights,
    parameters = parts$parameters, influence = parts$score %*% parts$bread/nrow(frame),
    native_influence = parts$score %*% model$naive.var, model = wrapper,
    levels = info$levels)
}

.suest_survey_survreg <- function(model) {
  if (!is.null(model$naive.var) || length(model$scale) != 1L ||
      !is.null(attr(stats::terms(model), "specials")$strata))
    stop("Survey survreg models need one common scale, without robust = TRUE.", call. = FALSE)
  mf <- model$model
  if (is.null(mf) || is.null(model$inv.info))
    stop("Could not recover the svysurvreg() estimation sample.", call. = FALSE)
  fit <- model
  fit$var <- model$inv.info
  class(fit) <- setdiff(class(model), "svysurvreg")
  parts <- .suest_survreg_components(fit)
  # svysurvreg()'s own influences use survival's scale residuals, which are
  # wrong for interval-censored rows; suest's scores recompute them.
  list(type = "survreg", engine = "survey::svysurvreg", design = model$survey.design,
    frame = mf, weights = stats::model.weights(mf), parameters = parts$parameters,
    influence = parts$score %*% parts$bread/nrow(mf),
    native_influence = stats::residuals(fit, "dfbeta", weighted = TRUE), model = model)
}

.suest_survey_ivreg <- function(model, env) {
  if (!requireNamespace("AER", quietly = TRUE))
    stop("Package 'AER' is required for survey::svyivreg() models.", call. = FALSE)
  call <- tryCatch(match.call(survey::svyivreg, model$call), error = function(e) model$call)
  design <- tryCatch(eval(call$design, env), error = function(e) NULL)
  if (!inherits(design, "survey.design2"))
    design <- tryCatch(eval(call$design, environment(stats::formula(model))),
      error = function(e) NULL)
  if (!inherits(design, "survey.design2"))
    stop("Could not recover the design of the svyivreg() model; keep the design object it was fitted with available.",
      call. = FALSE)
  mf <- model$model
  beta <- stats::coef(model)
  if (is.null(mf) || any(!is.finite(beta)))
    stop("Survey instrumental-variable models need finite, nonaliased coefficients and their model frame.",
      call. = FALSE)
  X <- stats::model.matrix(model, component = "projected")
  w <- stats::weights(model)
  score <- X * (w * stats::residuals(model))
  list(type = "ivreg", engine = "survey::svyivreg", design = design, frame = mf,
    weights = w, parameters = beta, influence = score %*% model$invinf, model = model)
}

.suest_survey_system <- function(models, model_names, observation_id,
    design, cluster, weight_type, call, env = parent.frame()) {
  if (!requireNamespace("survey", quietly = TRUE))
    stop("Package 'survey' is required for survey models.", call. = FALSE)
  if (!all(vapply(models, .suest_is_survey_model, logical(1))))
    stop("Survey systems require only survey models: svyglm(), svyolr(), svysurvreg(), svyivreg(), or svyVGAM::svy_vglm().", call. = FALSE)
  if (is.null(design))
    stop("Supply the full common survey design with survey_design =, including observations outside the model samples.", call. = FALSE)
  if (!is.null(cluster) || !is.null(weight_type))
    stop("For survey systems, supply clustering and weights through survey_design, without cluster or weight_type.", call. = FALSE)
  if (!is.character(observation_id) || !length(observation_id) ||
      anyNA(observation_id) || any(observation_id == "") || anyDuplicated(observation_id))
    stop("Survey systems require observation_id as one or more column names in the full design and each fitted design.", call. = FALSE)
  .suest_check_survey_design(design)
  lonely <- getOption("survey.lonely.psu", "fail")
  if (!lonely %in% c("fail", "remove", "certainty") ||
      isTRUE(getOption("survey.adjust.domain.lonely")))
    stop("Initial survey support requires survey.lonely.psu = 'fail', 'remove', or 'certainty' and survey.adjust.domain.lonely = FALSE.", call. = FALSE)

  full_ids <- .suest_survey_ids(design, observation_id, "survey_design")
  strata <- as.character(design$strata[[1L]])
  psu <- as.character(design$cluster[[1L]])
  psu_keys <- paste(encodeString(strata, quote = '"'),
    encodeString(psu, quote = '"'), sep = "\035")
  actual_psus <- vapply(split(psu, strata), function(x) length(unique(x)), integer(1))
  if (any(design$fpc$sampsize[, 1L] != actual_psus[strata]))
    stop("survey_design must retain all sampled PSUs; supply the full design before domain subsetting.", call. = FALSE)

  frames <- sample_ids <- sample_keys <- model_weights <- parameters <- influences <- vector("list", length(models))
  types <- engines <- character(length(models))
  levels <- vector("list", length(models))
  native_error <- numeric(length(models))
  same <- function(x, y) isTRUE(all.equal(unname(x), unname(y), tolerance = 1e-10, check.attributes = FALSE))
  pop <- function(d) if (is.null(d$fpc$popsize)) rep(Inf, nrow(d$variables)) else d$fpc$popsize[, 1L]
  for (i in seq_along(models)) {
    part <- .suest_survey_model(models[[i]], env)
    fitted_design <- part$design
    .suest_check_survey_design(fitted_design)
    ids <- .suest_survey_ids(fitted_design, observation_id, model_names[i])
    full_rows <- match(ids$keys, full_ids$keys)
    if (anyNA(full_rows))
      stop("Every fitted observation ID must occur in survey_design.", call. = FALSE)
    if (!identical(as.character(fitted_design$strata[[1L]]), strata[full_rows]) ||
        !identical(as.character(fitted_design$cluster[[1L]]), psu[full_rows]) ||
        !same(fitted_design$prob, design$prob[full_rows]) ||
        !same(fitted_design$fpc$sampsize[, 1L], design$fpc$sampsize[full_rows, 1L]) ||
        !same(pop(fitted_design), pop(design)[full_rows]))
      stop("Fitted models and survey_design must agree on weights, PSU IDs, strata, and FPCs for every fitted observation.", call. = FALSE)

    mf <- part$frame
    rows <- match(rownames(mf), rownames(fitted_design$variables))
    if (anyNA(rows) || anyDuplicated(rows) || nrow(part$influence) != length(rows))
      stop(sprintf("Could not align the estimation sample of model '%s' with its fitted survey design.", model_names[i]), call. = FALSE)
    raw_weights <- 1 / fitted_design$prob[rows]
    w <- part$weights
    if (!is.null(w) && (length(w) != length(rows) || !same(w / mean(w), raw_weights / mean(raw_weights))))
      stop("The model fitting weights do not match its survey sampling weights.", call. = FALSE)
    align <- function(z) {
      out <- matrix(0, nrow(design$variables), ncol(z))
      out[full_rows[rows], ] <- z
      out
    }
    aligned <- align(part$influence)
    native_influence <- if (is.null(part$native_influence)) aligned else align(part$native_influence)
    native <- stats::vcov(models[[i]])
    reproduced <- .suest_survey_vcov(native_influence, design)
    native_error[i] <- max(abs(.suest_survey_vcov(aligned, design) - native))
    if (!all(is.finite(native)) || !same(reproduced, native))
      stop(sprintf("Could not reproduce the native design-based covariance for model '%s'; check the full design, fitting variance method, and survey options used at estimation.", model_names[i]), call. = FALSE)
    # Preserve all original columns for IDs, averaging weights, and predictions.
    extra <- setdiff(names(fitted_design$variables), names(mf))
    mf[extra] <- fitted_design$variables[rows, extra, drop = FALSE]
    frames[[i]] <- mf
    sample_ids[[i]] <- ids$values[rows, , drop = FALSE]
    sample_keys[[i]] <- ids$keys[rows]
    model_weights[[i]] <- raw_weights
    parameters[[i]] <- part$parameters
    influences[[i]] <- aligned
    types[i] <- part$type
    engines[i] <- part$engine
    levels[i] <- list(part$levels)
    models[[i]] <- part$model
  }
  names(frames) <- names(sample_ids) <- names(sample_keys) <- names(model_weights) <- model_names
  local_names <- lapply(parameters, names)
  joint_names <- unlist(Map(function(nm, terms) paste0(nm, "::", terms), model_names, local_names), use.names = FALSE)
  b <- stats::setNames(unlist(parameters, use.names = FALSE), joint_names)
  U <- do.call(cbind, influences)
  V <- .suest_survey_vcov(U, design)
  dimnames(V) <- list(joint_names, joint_names)
  p <- lengths(parameters)
  index <- Map(function(start, size) start - 1L + seq_len(size), cumsum(c(1L, p[-length(p)])), p)
  names(index) <- model_names
  union_keys <- unique(unlist(sample_keys, use.names = FALSE))
  n_model <- stats::setNames(vapply(frames, nrow, integer(1)), model_names)
  types <- stats::setNames(types, model_names)
  categorical <- types %in% c("ologit", "oprobit", "multinom", "gologit")
  levels[!categorical] <- list(NULL)
  pair <- if (all(types == "survey_lm")) {
    list(key = .suest_pair_key(types), mixed = FALSE, scale = "fitted values")
  } else {
    .suest_pair_info(types)
  }
  out <- list(models = models, model_frames = frames, model_weights = model_weights,
    weight_type = "survey", coefficients = b, vcov = V, model_names = model_names,
    local_names = local_names, index = index, model_types = types,
    model_engines = stats::setNames(engines, model_names),
    category_levels = levels, pair_key = pair$key,
    mixed_models = pair$mixed, comparison_scale = pair$scale, theta = NULL,
    sample_sources = rep("<common survey design>", length(models)),
    sample_rows = lapply(frames, rownames), sample_ids = sample_ids,
    cluster_ids = lapply(sample_keys, function(keys) design$cluster[match(keys, full_ids$keys), , drop = FALSE]),
    cluster_keys = psu_keys[match(union_keys, full_ids$keys)], n_clusters = length(unique(psu_keys)),
    sample_keys = sample_keys, nobs_models = n_model, nobs_union = length(union_keys),
    nobs_overlap = length(Reduce(intersect, sample_keys)),
    nobs_overlap_pairwise = .suest_pairwise_overlap(sample_keys, model_names),
    survey = list(design_df = survey::degf(design), n_strata = length(unique(strata)),
      n_design = nrow(design$variables), lonely_psu = lonely,
      native_vcov_max_abs_error = stats::setNames(native_error, model_names),
      inference = "Asymptotic normal by default; design_df is metadata, not an automatic t adjustment."),
    call = call)
  class(out) <- "suest_model"
  out
}
