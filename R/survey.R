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
  if (identical(family, "quasibinomial") && link %in% c("logit", "probit")) return(link)
  "unsupported"
}

# Direct reconstruction from the retained GLM fitting state, independent of
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

.suest_survey_system <- function(models, model_names, observation_id,
    design, cluster, weight_type, call) {
  if (!requireNamespace("survey", quietly = TRUE))
    stop("Package 'survey' is required for survey models.", call. = FALSE)
  if (!all(vapply(models, inherits, logical(1), "svyglm")))
    stop("Survey systems currently require only survey::svyglm() models.", call. = FALSE)
  types <- vapply(models, .suest_survey_type, character(1))
  if (any(types == "unsupported"))
    stop("Survey support is restricted to Gaussian identity and quasibinomial logit/probit svyglm() models.", call. = FALSE)
  if (any(types == "survey_lm") && any(types != "survey_lm"))
    stop("This survey increment does not combine Gaussian and binary survey models in one system.", call. = FALSE)
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
  native_error <- numeric(length(models))
  for (i in seq_along(models)) {
    model <- models[[i]]
    type <- types[i]
    if (!isTRUE(model$converged) || any(!is.finite(model$coefficients)))
      stop("Survey models must converge and have finite, nonaliased coefficients.", call. = FALSE)
    if (!is.null(model$call$weights))
      stop("Additional svyglm() weights are unsupported; specify sampling weights only in svydesign().", call. = FALSE)
    if (!is.null(model$call$offset))
      stop("For survey models, specify offsets in the formula with offset(), not the offset argument.", call. = FALSE)
    fitted_design <- model$survey.design
    .suest_check_survey_design(fitted_design)
    ids <- .suest_survey_ids(fitted_design, observation_id, model_names[i])
    full_rows <- match(ids$keys, full_ids$keys)
    if (anyNA(full_rows))
      stop("Every fitted observation ID must occur in survey_design.", call. = FALSE)
    same <- function(x, y) isTRUE(all.equal(unname(x), unname(y), tolerance = 1e-10, check.attributes = FALSE))
    pop <- function(d) if (is.null(d$fpc$popsize)) rep(Inf, nrow(d$variables)) else d$fpc$popsize[, 1L]
    if (!identical(as.character(fitted_design$strata[[1L]]), strata[full_rows]) ||
        !identical(as.character(fitted_design$cluster[[1L]]), psu[full_rows]) ||
        !same(fitted_design$prob, design$prob[full_rows]) ||
        !same(fitted_design$fpc$sampsize[, 1L], design$fpc$sampsize[full_rows, 1L]) ||
        !same(pop(fitted_design), pop(design)[full_rows]))
      stop("Fitted models and survey_design must agree on weights, PSU IDs, strata, and FPCs for every fitted observation.", call. = FALSE)

    mf <- model$model
    if (is.null(mf))
      stop("Refit svyglm() with model = TRUE to retain its estimation sample.", call. = FALSE)
    rows <- match(rownames(mf), rownames(fitted_design$variables))
    if (anyNA(rows) || anyDuplicated(rows) || length(rows) != nrow(fitted_design$variables))
      stop("Could not align the retained svyglm() model frame with its fitted survey design.", call. = FALSE)
    X <- stats::model.matrix(model)
    y <- if (type == "survey_lm") stats::model.response(mf) else model$y
    w <- as.numeric(model$prior.weights)
    beta <- stats::coef(model)
    basic_bad <- nrow(X) != nrow(mf) || ncol(X) != length(beta) || nrow(X) <= ncol(X) ||
      length(w) != nrow(X) || any(!is.finite(w) | w <= 0)
    if (type == "survey_lm") {
      if (!is.numeric(y) || is.matrix(y) || any(!is.finite(y)) || basic_bad)
        stop("Survey linear models require a full-rank numeric response, positive weights, and residual observations.", call. = FALSE)
    } else if (is.null(y) || !is.numeric(y) || is.matrix(y) || any(!is.finite(y)) ||
               any(!y %in% c(0, 1)) || basic_bad) {
      stop("Survey binary models require a full-rank numeric 0/1 response, positive weights, and residual observations.", call. = FALSE)
    }
    raw_weights <- 1 / fitted_design$prob[rows]
    if (!same(w / mean(w), raw_weights / mean(raw_weights)))
      stop("The model fitting weights do not match its survey sampling weights.", call. = FALSE)
    offset <- stats::model.offset(mf)
    if (is.null(offset)) offset <- rep(0, nrow(X))
    if (type == "survey_lm") {
      residual <- y - as.numeric(X %*% beta) - offset
      influence <- (X * (w * residual)) %*% solve(crossprod(X, X * w))
    } else {
      influence <- .suest_survey_binary_influence(model, X, y, w, beta, offset)
    }
    aligned <- matrix(0, nrow(design$variables), ncol(X))
    aligned[full_rows[rows], ] <- influence
    native <- stats::vcov(model)
    reproduced <- .suest_survey_vcov(aligned, design)
    native_error[i] <- max(abs(reproduced - native))
    if (!all(is.finite(native)) || !same(reproduced, native))
      stop(sprintf("Could not reproduce the native design-based covariance for model '%s'; check the full design, fitting variance method, and survey options used at estimation.", model_names[i]), call. = FALSE)
    # Preserve all original columns for IDs, averaging weights, and predictions.
    extra <- setdiff(names(fitted_design$variables), names(mf))
    mf[extra] <- fitted_design$variables[rows, extra, drop = FALSE]
    frames[[i]] <- mf
    sample_ids[[i]] <- ids$values[rows, , drop = FALSE]
    sample_keys[[i]] <- ids$keys[rows]
    model_weights[[i]] <- raw_weights
    parameters[[i]] <- beta
    influences[[i]] <- aligned
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
  pair <- if (all(types == "survey_lm")) {
    list(key = .suest_pair_key(types), mixed = FALSE, scale = "fitted values")
  } else {
    .suest_pair_info(types)
  }
  out <- list(models = models, model_frames = frames, model_weights = model_weights,
    weight_type = "survey", coefficients = b, vcov = V, model_names = model_names,
    local_names = local_names, index = index, model_types = types,
    model_engines = stats::setNames(rep("survey::svyglm", length(models)), model_names),
    category_levels = rep(list(NULL), length(models)), pair_key = pair$key,
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
