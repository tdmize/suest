.suest_pair_key <- function(types) {
  paste(sort(types), collapse = "+")
}

.suest_pair_info <- function(types) {
  key <- .suest_pair_key(types)

  same_type <- length(unique(types)) == 1L
  scalar <- types %in% c(
    "lm", "logit", "probit", "cloglog", "poisson", "negbin", "glm",
    "fraclogit", "fracprobit", "fraccloglog", "fracloglog", "survreg", "betareg",
    "zip", "zinb", "truncreg", "censreg", "ivreg", "ivprobit", "hetprobit",
    "hetlogit", "biprobit", "panel_fe", "panel_be", "panel_re", "panel_ml",
    "panel_gee"
  )
  categorical <- types %in% c("ologit", "oprobit", "multinom")

  scale <- if (all(categorical)) {
    "category probabilities"
  } else if (all(types %in% c("poisson", "negbin"))) {
    "expected counts"
  } else if (all(types %in% c(
               "logit", "probit", "cloglog", "ivprobit", "biprobit"
             ))) {
    "predicted probabilities"
  } else if (all(types %in% c("hetprobit", "hetlogit"))) {
    "predicted probabilities"
  } else if (all(types %in% c(
               "fraclogit", "fracprobit", "fraccloglog", "fracloglog"
             ))) {
    "fitted proportions"
  } else if (all(types == "betareg")) {
    "fitted proportions"
  } else if (all(types %in% c("zip", "zinb"))) {
    "expected counts"
  } else if (all(types %in% c(
               "lm", "panel_fe", "panel_be", "panel_re", "panel_ml"
             ))) {
    "fitted values"
  } else {
    "model-specific response values"
  }

  list(
    key = key,
    mixed = !same_type,
    scale = scale
  )
}

.suest_data_source <- function(model) {
  data_call <- model$call$data
  env <- if (inherits(model, "mvProbit")) {
    fit_environment <- environment(model$objectiveFn)
    source_formula <- if (is.environment(fit_environment) &&
                          exists("formula", fit_environment, inherits = FALSE)) {
      get("formula", fit_environment, inherits = FALSE)
    } else {
      stats::formula(model)
    }
    environment(source_formula)
  } else {
    environment(stats::formula(model))
  }
  env_key <- format(env)

  if (is.null(data_call)) {
    paste0("<formula environment>@", env_key)
  } else {
    paste0(
      paste(deparse(data_call, width.cutoff = 500L), collapse = ""),
      "@",
      env_key
    )
  }
}

.suest_id_keys <- function(x, n, model_name) {
  if (is.atomic(x) && is.null(dim(x))) {
    x <- data.frame(id = x, check.names = FALSE)
  } else if (is.matrix(x)) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  }

  if (!is.data.frame(x) || ncol(x) == 0L)
    stop(
      paste0(
        "Observation IDs for model '", model_name, "' must be a vector, ",
        "matrix, or data frame with at least one column."
      ),
      call. = FALSE
    )

  if (nrow(x) != n)
    stop(
      sprintf(
        "Observation IDs for model '%s' must have one row per estimation-sample observation.",
        model_name
      ),
      call. = FALSE
    )

  invalid_column <- vapply(
    x,
    function(column) !is.atomic(column) || !is.null(dim(column)),
    logical(1)
  )
  if (any(invalid_column))
    stop("Observation-ID columns must be atomic vectors.", call. = FALSE)

  if (any(vapply(x, anyNA, logical(1))))
    stop(
      sprintf("Observation IDs for model '%s' cannot be missing.", model_name),
      call. = FALSE
    )

  encoded <- lapply(x, function(column) {
    encodeString(as.character(column), quote = "\"")
  })
  keys <- do.call(paste, c(encoded, sep = "\035"))

  if (anyDuplicated(keys))
    stop(
      sprintf(
        "Observation IDs must be unique within model '%s'.",
        model_name
      ),
      call. = FALSE
    )

  list(
    keys = paste0("<observation_id>\r", keys),
    values = x
  )
}

.suest_observation_ids_from_model <- function(
    model,
    model_frame,
    variables,
    model_name) {
  if (all(variables %in% names(model_frame))) {
    values <- model_frame[, variables, drop = FALSE]
    return(.suest_id_keys(values, nrow(model_frame), model_name))
  }

  data_call <- model$call$data
  if (is.null(data_call))
    stop(
      paste0(
        "Could not retrieve observation-ID columns for model '",
        model_name,
        "'. Supply model-specific ID vectors in a list instead."
      ),
      call. = FALSE
    )

  data_environment <- environment(stats::formula(model))
  data <- try(eval(data_call, envir = data_environment), silent = TRUE)
  if (inherits(data, "try-error") || !is.data.frame(data))
    stop(
      paste0(
        "Could not evaluate the original data for model '",
        model_name,
        "'. Supply model-specific ID vectors in a list instead."
      ),
      call. = FALSE
    )

  missing_variables <- setdiff(variables, names(data))
  if (length(missing_variables))
    stop(
      sprintf(
        "Observation-ID column(s) not found for model '%s': %s.",
        model_name,
        paste(missing_variables, collapse = ", ")
      ),
      call. = FALSE
    )

  index <- match(rownames(model_frame), rownames(data))
  if (anyNA(index))
    stop(
      paste0(
        "Observation-ID columns could not be aligned to the estimation ",
        "sample for model '", model_name, "'. Supply model-specific ID ",
        "vectors in a list instead."
      ),
      call. = FALSE
    )

  values <- data[index, variables, drop = FALSE]
  .suest_id_keys(values, nrow(model_frame), model_name)
}

.suest_observation_ids <- function(
    models,
    model_frames,
    observation_id,
    model_names) {
  if (is.character(observation_id) && length(observation_id) >= 1L) {
    if (anyNA(observation_id) || any(observation_id == "") ||
        anyDuplicated(observation_id))
      stop(
        "'observation_id' column names must be unique and nonempty.",
        call. = FALSE
      )

    return(Map(
      function(model, model_frame, model_name) {
        .suest_observation_ids_from_model(
          model,
          model_frame,
          observation_id,
          model_name
        )
      },
      models,
      model_frames,
      model_names
    ))
  }

  if (!is.list(observation_id) || is.data.frame(observation_id) ||
      length(observation_id) != length(models))
    stop(
      paste0(
        "'observation_id' must be NULL, one or more ID column names, or ",
        "a list containing one ID vector, matrix, or data frame per model."
      ),
      call. = FALSE
    )

  id_names <- names(observation_id)
  if (!is.null(id_names) && any(nzchar(id_names))) {
    if (any(!nzchar(id_names)) || !setequal(id_names, model_names))
      stop(
        "Named 'observation_id' lists must use the model names.",
        call. = FALSE
      )
    observation_id <- observation_id[model_names]
  }

  Map(
    .suest_id_keys,
    observation_id,
    vapply(model_frames, nrow, integer(1)),
    model_names
  )
}

.suest_cluster_keys <- function(x, n, model_name) {
  if (is.atomic(x) && is.null(dim(x))) {
    x <- data.frame(cluster = x, check.names = FALSE)
  } else if (is.matrix(x)) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  }

  if (!is.data.frame(x) || ncol(x) == 0L)
    stop(
      paste0(
        "Cluster IDs for model '", model_name, "' must be a vector, ",
        "matrix, or data frame with at least one column."
      ),
      call. = FALSE
    )

  if (nrow(x) != n)
    stop(
      sprintf(
        "Cluster IDs for model '%s' must have one row per estimation-sample observation.",
        model_name
      ),
      call. = FALSE
    )

  invalid_column <- vapply(
    x,
    function(column) !is.atomic(column) || !is.null(dim(column)),
    logical(1)
  )
  if (any(invalid_column))
    stop("Cluster-ID columns must be atomic vectors.", call. = FALSE)

  if (any(vapply(x, anyNA, logical(1))))
    stop(
      sprintf("Cluster IDs for model '%s' cannot be missing.", model_name),
      call. = FALSE
    )

  encoded <- lapply(x, function(column) {
    encodeString(as.character(column), quote = "\"")
  })
  keys <- do.call(paste, c(encoded, sep = "\035"))

  list(
    keys = paste0("<cluster>\r", keys),
    values = x
  )
}

.suest_clusters_from_model <- function(
    model,
    model_frame,
    variables,
    model_name) {
  if (all(variables %in% names(model_frame))) {
    values <- model_frame[, variables, drop = FALSE]
    return(.suest_cluster_keys(values, nrow(model_frame), model_name))
  }

  data_call <- model$call$data
  if (is.null(data_call))
    stop(
      paste0(
        "Could not retrieve cluster-ID columns for model '", model_name,
        "'. Supply model-specific cluster vectors in a list instead."
      ),
      call. = FALSE
    )

  data_environment <- environment(stats::formula(model))
  data <- try(eval(data_call, envir = data_environment), silent = TRUE)
  if (inherits(data, "try-error") || !is.data.frame(data))
    stop(
      paste0(
        "Could not evaluate the original data for model '", model_name,
        "'. Supply model-specific cluster vectors in a list instead."
      ),
      call. = FALSE
    )

  missing_variables <- setdiff(variables, names(data))
  if (length(missing_variables))
    stop(
      sprintf(
        "Cluster-ID column(s) not found for model '%s': %s.",
        model_name,
        paste(missing_variables, collapse = ", ")
      ),
      call. = FALSE
    )

  index <- match(rownames(model_frame), rownames(data))
  if (anyNA(index))
    stop(
      paste0(
        "Cluster-ID columns could not be aligned to the estimation sample ",
        "for model '", model_name, "'. Supply model-specific cluster ",
        "vectors in a list instead."
      ),
      call. = FALSE
    )

  values <- data[index, variables, drop = FALSE]
  .suest_cluster_keys(values, nrow(model_frame), model_name)
}

.suest_clusters <- function(models, model_frames, cluster, model_names) {
  if (is.character(cluster) && length(cluster) >= 1L) {
    if (anyNA(cluster) || any(cluster == "") || anyDuplicated(cluster))
      stop(
        "'cluster' column names must be unique and nonempty.",
        call. = FALSE
      )

    return(Map(
      function(model, model_frame, model_name) {
        .suest_clusters_from_model(
          model,
          model_frame,
          cluster,
          model_name
        )
      },
      models,
      model_frames,
      model_names
    ))
  }

  if (!is.list(cluster) || is.data.frame(cluster) ||
      length(cluster) != length(models))
    stop(
      paste0(
        "'cluster' must be NULL, one or more cluster column names, or a ",
        "list containing one cluster vector, matrix, or data frame per model."
      ),
      call. = FALSE
    )

  cluster_names <- names(cluster)
  if (!is.null(cluster_names) && any(nzchar(cluster_names))) {
    if (any(!nzchar(cluster_names)) || !setequal(cluster_names, model_names))
      stop(
        "Named 'cluster' lists must use the model names.",
        call. = FALSE
      )
    cluster <- cluster[model_names]
  }

  Map(
    .suest_cluster_keys,
    cluster,
    vapply(model_frames, nrow, integer(1)),
    model_names
  )
}

.suest_panel_id <- function(model) {
  if (inherits(model, "geeglm")) {
    model$id
  } else if (inherits(model, "plm")) {
    plm::index(model)[[1L]]
  } else if (inherits(model, "lme")) {
    model$groups[[1L]]
  } else {
    stop("Internal error: unsupported panel-model object.", call. = FALSE)
  }
}

.suest_panel_clusters <- function(models, model_frames, model_names) {
  Map(
    function(model, model_frame, model_name) {
      panel <- .suest_panel_id(model)
      .suest_cluster_keys(panel, nrow(model_frame), model_name)
    },
    models,
    model_frames,
    model_names
  )
}

.suest_validate_panel_cluster_nesting <- function(
    models,
    cluster_info,
    model_names) {
  for (i in seq_along(models)) {
    panel <- as.character(.suest_panel_id(models[[i]]))
    cluster <- cluster_info[[i]]$keys
    mapping <- unique(data.frame(panel = panel, cluster = cluster))
    if (anyDuplicated(mapping$panel))
      stop(
        sprintf(
          "Panel IDs must be nested within clusters for model '%s'.",
          model_names[i]
        ),
        call. = FALSE
      )
  }

  invisible(TRUE)
}

.suest_union_clusters <- function(
    cluster_info,
    sample_keys,
    union_keys,
    model_names) {
  union_clusters <- rep.int(NA_character_, length(union_keys))

  for (i in seq_along(sample_keys)) {
    index <- match(sample_keys[[i]], union_keys)
    candidate <- cluster_info[[i]]$keys
    conflict <- !is.na(union_clusters[index]) &
      union_clusters[index] != candidate

    if (any(conflict))
      stop(
        sprintf(
          paste0(
            "Cluster IDs disagree for observations shared with model '%s'. ",
            "The first conflict is at overlapping row '%s'."
          ),
          model_names[i],
          sub("^[^\r]*\r", "", sample_keys[[i]][which(conflict)[1L]])
        ),
        call. = FALSE
      )

    missing <- is.na(union_clusters[index])
    union_clusters[index[missing]] <- candidate[missing]
  }

  if (anyNA(union_clusters))
    stop("Internal error while aligning cluster IDs.", call. = FALSE)

  union_clusters
}

.suest_extract_parameters <- function(model, type, engine) {
  engine <- .suest_engine_name(engine)
  if (type == "panel_fe") {
    .suest_plm_fe_parameters(model)
  } else if (type == "panel_ml") {
    .suest_lme_parameters(model)
  } else if (type == "biprobit") {
    raw <- stats::coef(model)
    c(raw[-length(raw)], athrho = atanh(unname(raw[length(raw)])))
  } else if (type == "ivprobit") {
    parameters <- stats::coef(model)
    names(parameters)[names(parameters) == "atanhrho"] <- "athrho"
    parameters
  } else if (identical(engine, "ordinal::clm")) {
    stats::coef(model)
  } else if (type %in% c("ologit", "oprobit")) {
    c(model$coefficients, model$zeta)
  } else if (type == "multinom") {
    cf <- stats::coef(model)
    b <- as.vector(t(cf))
    names(b) <- as.vector(t(outer(
      rownames(cf),
      colnames(cf),
      paste,
      sep = ":"
    )))
    b
  } else if (type == "negbin") {
    c(stats::coef(model), ln_theta = log(model$theta))
  } else if (type == "survreg") {
    parameters <- stats::coef(model)
    covariance_names <- colnames(stats::vcov(model))
    if ("Log(scale)" %in% covariance_names)
      parameters <- c(parameters, "Log(scale)" = log(model$scale))
    parameters
  } else if (type == "zinb") {
    c(stats::coef(model), ln_theta = log(model$theta))
  } else {
    stats::coef(model)
  }
}

.suest_lme_parameters <- function(model) {
  sigma_e <- as.numeric(model$sigma)
  variance <- nlme::pdMatrix(model$modelStruct$reStruct[[1L]])
  sigma_u <- sigma_e * sqrt(variance[1L, 1L])
  c(model$coefficients$fixed, sigma_u = sigma_u, sigma_e = sigma_e)
}

.suest_plm_fe_parameters <- function(model) {
  beta <- stats::coef(model)
  X <- stats::model.matrix(model, model = "pooling")
  y <- as.numeric(plm::pmodel.response(model, model = "pooling"))
  missing <- setdiff(names(beta), colnames(X))
  if (length(missing))
    stop(
      "Unable to align the plm fixed-effects coefficients to its model matrix.",
      call. = FALSE
    )

  means <- colMeans(X[, names(beta), drop = FALSE])
  intercept <- mean(y) - sum(means * beta)
  c(`(Intercept)` = intercept, beta)
}


.suest_validate_weight_type <- function(weight_type) {
  if (is.null(weight_type))
    return(NULL)

  if (length(weight_type) != 1L || is.na(weight_type))
    stop(
      "'weight_type' must be NULL or \"pweight\".",
      call. = FALSE
    )

  weight_type <- as.character(weight_type)
  if (!identical(weight_type, "pweight"))
    stop(
      "'weight_type' must be NULL or \"pweight\".",
      call. = FALSE
    )

  weight_type
}

.suest_model_weights <- function(model_frame, model = NULL, engine = NULL) {
  engine <- if (is.null(engine)) NULL else .suest_engine_name(engine)
  weights <- if (engine %in% c(
                   "nnet::multinom", "betareg::betareg", "pscl::zeroinfl",
                   "fixest::feols"
                 ) &&
                 !is.null(model$weights)) {
    model$weights
  } else {
    stats::model.weights(model_frame)
  }
  if (is.null(weights))
    weights <- rep.int(1, nrow(model_frame))

  weights <- as.numeric(weights)
  if (length(weights) != nrow(model_frame))
    stop(
      "The model weights do not align with the estimation sample.",
      call. = FALSE
    )

  weights
}

.suest_validate_pweights <- function(weights, model_name) {
  if (anyNA(weights) || any(!is.finite(weights)))
    stop(
      sprintf(
        "Pweights for model '%s' must be finite and nonmissing.",
        model_name
      ),
      call. = FALSE
    )

  if (any(weights <= 0))
    stop(
      sprintf(
        "Pweights for model '%s' must be strictly positive.",
        model_name
      ),
      call. = FALSE
    )

  invisible(TRUE)
}

.suest_weights_equal <- function(x, y) {
  if (length(x) != length(y))
    return(FALSE)

  tolerance <- sqrt(.Machine$double.eps)
  scale <- pmax(abs(x), abs(y), .Machine$double.eps)
  all(abs(x - y) <= tolerance * scale)
}

.suest_validate_overlap_weights <- function(
    model_weights,
    sample_keys,
    model_names) {
  pairs <- utils::combn(seq_along(sample_keys), 2L, simplify = FALSE)

  for (pair in pairs) {
    overlap_keys <- intersect(sample_keys[[pair[1L]]], sample_keys[[pair[2L]]])
    if (!length(overlap_keys))
      next

    index1 <- match(overlap_keys, sample_keys[[pair[1L]]])
    index2 <- match(overlap_keys, sample_keys[[pair[2L]]])
    weights1 <- model_weights[[pair[1L]]][index1]
    weights2 <- model_weights[[pair[2L]]][index2]

    if (!.suest_weights_equal(weights1, weights2)) {
      difference <- abs(weights1 - weights2)
      tolerance <- sqrt(.Machine$double.eps) *
        pmax(abs(weights1), abs(weights2), .Machine$double.eps)
      mismatch <- which(difference > tolerance)[1L]

      stop(
        sprintf(
          paste0(
            "Pweights must agree for observations included in both models. ",
            "The first mismatch is between models '%s' and '%s' at ",
            "overlapping row '%s'."
          ),
          model_names[pair[1L]],
          model_names[pair[2L]],
          sub("^[^\r]*\r", "", overlap_keys[mismatch])
        ),
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

.suest_lm_pweight_components <- function(model) {
  mf <- stats::model.frame(model)
  X <- stats::model.matrix(model)
  y <- stats::model.response(mf)
  weights <- .suest_model_weights(mf)
  parameters <- stats::coef(model)
  offset <- stats::model.offset(mf)
  if (is.null(offset))
    offset <- 0
  residual <- as.numeric(y - X %*% parameters - offset)

  n <- nrow(X)
  k <- ncol(X)
  sum_weights <- sum(weights)
  sigma2 <- sum(weights * residual^2) / (sum_weights - k)
  if (!is.finite(sigma2) || sigma2 <= 0)
    stop(
      paste0(
        "Unable to construct a positive pweighted linear-model residual ",
        "variance. The sum of the weights must exceed the number of ",
        "estimated mean parameters."
      ),
      call. = FALSE
    )

  # suest2 refits a plain-pweight regression with iweights and lets official
  # suest append its Gaussian lnvar equation. The resulting score uses the
  # evaluated weights, while the ancillary bread uses Stata's integer-valued
  # effective iweight count.
  U_beta <- X * (weights * residual / sigma2)
  U_lnvar <- 0.5 * weights * (residual^2 / sigma2 - 1)
  U <- cbind(U_beta, lnvar = U_lnvar)

  parameters <- c(parameters, lnvar = log(sigma2))
  colnames(U) <- names(parameters)
  rownames(U) <- rownames(X)

  A <- crossprod(X, X * weights)
  effective_n <- floor(sum_weights)
  if (!is.finite(effective_n) || effective_n < 1)
    stop(
      "The pweighted linear model has no positive effective weight count.",
      call. = FALSE
    )

  B <- matrix(0, nrow = k + 1L, ncol = k + 1L)
  B[seq_len(k), seq_len(k)] <-
    n * sigma2 * .suest_information_inverse(A)
  B[k + 1L, k + 1L] <- 2 * n / effective_n
  dimnames(B) <- list(names(parameters), names(parameters))

  list(
    score = U,
    bread = B,
    parameters = parameters
  )
}


.suest_lm_components <- function(model) {
  mf <- stats::model.frame(model)
  X <- stats::model.matrix(model)
  residual <- as.numeric(stats::residuals(model))
  parameters <- stats::coef(model)
  n <- nrow(X)
  k <- ncol(X)

  # Stata's regress/suest contract treats the residual variance as an
  # ancillary parameter. The reported point estimate is e(rmse)^2, namely
  # RSS / df_residual rather than the Gaussian maximum-likelihood RSS / N.
  sigma2 <- sum(residual^2) / stats::df.residual(model)
  if (!is.finite(sigma2) || sigma2 <= 0)
    stop(
      "Unable to construct a positive linear-model residual variance.",
      call. = FALSE
    )

  U_beta <- X * as.numeric(residual / sigma2)
  U_lnvar <- 0.5 * (residual^2 / sigma2 - 1)
  U <- cbind(U_beta, lnvar = U_lnvar)

  parameters <- c(parameters, lnvar = log(sigma2))
  colnames(U) <- names(parameters)
  rownames(U) <- rownames(X)

  B <- matrix(0, nrow = k + 1L, ncol = k + 1L)
  B[seq_len(k), seq_len(k)] <-
    n * sigma2 * .suest_information_inverse(crossprod(X))
  B[k + 1L, k + 1L] <- 2
  dimnames(B) <- list(names(parameters), names(parameters))

  list(score = U, bread = B, parameters = parameters)
}

.suest_plm_fe_components <- function(model) {
  parameters <- .suest_plm_fe_parameters(model)
  beta_names <- setdiff(names(parameters), "(Intercept)")
  X <- stats::model.matrix(model, model = "within")
  missing <- setdiff(beta_names, colnames(X))
  if (length(missing))
    stop(
      "Unable to align the plm within-model matrix to its coefficients.",
      call. = FALSE
    )
  X <- X[, beta_names, drop = FALSE]
  residual <- as.numeric(stats::residuals(model))
  n <- nrow(X)
  p <- ncol(X)

  U <- cbind(`(Intercept)` = 0, X * residual)
  colnames(U) <- names(parameters)
  rownames(U) <- rownames(X)

  X_pooling <- stats::model.matrix(model, model = "pooling")
  means <- colMeans(X_pooling[, beta_names, drop = FALSE])
  mapping <- rbind(-means, diag(p))
  rownames(mapping) <- names(parameters)
  colnames(mapping) <- beta_names

  B <- matrix(0, nrow = p + 1L, ncol = p + 1L)
  B[, -1L] <- n * mapping %*%
    .suest_information_inverse(crossprod(X))
  dimnames(B) <- list(names(parameters), names(parameters))

  list(
    score = U,
    bread = B,
    parameters = parameters,
    cluster_df_k = p + 1L
  )
}

.suest_plm_be_components <- function(model) {
  parameters <- stats::coef(model)
  X_panel <- stats::model.matrix(model)
  residual_panel <- as.numeric(stats::residuals(model))
  panel_rows <- rownames(X_panel)
  panel <- as.character(plm::index(model)[[1L]])
  first <- match(panel_rows, panel)
  if (anyNA(first) || length(residual_panel) != nrow(X_panel))
    stop(
      "Unable to align the plm between-model panel means to its sample.",
      call. = FALSE
    )

  U <- matrix(0, nrow = length(panel), ncol = ncol(X_panel))
  U[first, ] <- X_panel * residual_panel
  colnames(U) <- names(parameters)
  rownames(U) <- rownames(stats::model.frame(model))

  n <- nrow(U)
  B <- n * .suest_information_inverse(crossprod(X_panel))
  dimnames(B) <- list(names(parameters), names(parameters))

  list(
    score = U,
    bread = B,
    parameters = parameters,
    cluster_df_k = length(parameters),
    cluster_df_n = nrow(X_panel)
  )
}

.suest_plm_re_components <- function(model) {
  parameters <- stats::coef(model)
  X <- stats::model.matrix(model)
  residual <- as.numeric(stats::residuals(model))
  U <- X * residual
  colnames(U) <- names(parameters)
  rownames(U) <- rownames(stats::model.frame(model))

  n <- nrow(X)
  B <- n * .suest_information_inverse(crossprod(X))
  dimnames(B) <- list(names(parameters), names(parameters))

  list(
    score = U,
    bread = B,
    parameters = parameters,
    cluster_df_k = length(parameters),
    cluster_df_n = n
  )
}

.suest_lme_ml_components <- function(model) {
  parameters <- .suest_lme_parameters(model)
  beta_names <- names(model$coefficients$fixed)
  beta <- parameters[beta_names]
  sigma_u <- unname(parameters["sigma_u"])
  sigma_e <- unname(parameters["sigma_e"])
  if (!is.finite(sigma_u) || !is.finite(sigma_e) ||
      sigma_u <= 0 || sigma_e <= 0)
    stop(
      "Random-intercept ML variance components must be positive and finite.",
      call. = FALSE
    )

  frame <- .suest_model_frame(model, "nlme::lme")
  terms <- stats::delete.response(stats::terms(model))
  X <- stats::model.matrix(
    terms,
    frame,
    contrasts.arg = model$contrasts
  )
  missing <- setdiff(beta_names, colnames(X))
  if (length(missing))
    stop(
      "Unable to align the nlme::lme fixed-effects design matrix.",
      call. = FALSE
    )
  X <- X[, beta_names, drop = FALSE]
  y <- as.numeric(stats::model.response(frame))
  panel <- as.character(.suest_panel_id(model))
  n <- nrow(X)
  p <- ncol(X)
  if (length(y) != n || length(panel) != n)
    stop(
      "Unable to align the nlme::lme response and panel identifiers.",
      call. = FALSE
    )

  residual <- y - as.numeric(X %*% beta)
  U <- matrix(
    0,
    nrow = n,
    ncol = p + 2L,
    dimnames = list(rownames(frame), names(parameters))
  )
  information <- matrix(
    0,
    nrow = p + 2L,
    ncol = p + 2L,
    dimnames = list(names(parameters), names(parameters))
  )
  sigma_e2 <- sigma_e^2
  sigma_u2 <- sigma_u^2

  for (rows in split(seq_len(n), panel)) {
    Xi <- X[rows, , drop = FALSE]
    ri <- residual[rows]
    Ti <- length(rows)
    residual_sum <- sum(ri)
    denominator <- sigma_e2 + Ti * sigma_u2
    q <- ri / sigma_e2 -
      sigma_u2 * residual_sum / (sigma_e2 * denominator)
    q_squared <- sum(q^2)

    U[rows, seq_len(p)] <- Xi * q
    last <- rows[Ti]
    U[last, p + 1L] <- sigma_u * (
      residual_sum^2 / denominator^2 - Ti / denominator
    )
    U[last, p + 2L] <- sigma_e * (
      q_squared - (Ti - 1) / sigma_e2 - 1 / denominator
    )

    column_sums <- colSums(Xi)
    inverse_covariance_X <- Xi / sigma_e2 -
      matrix(rep(column_sums, each = Ti), nrow = Ti) *
      sigma_u2 / (sigma_e2 * denominator)
    information[seq_len(p), seq_len(p)] <-
      information[seq_len(p), seq_len(p)] +
      crossprod(Xi, inverse_covariance_X)

    beta_sigma_u <-
      2 * sigma_u * residual_sum * column_sums / denominator^2
    inverse_covariance_q <- q / sigma_e2 -
      sigma_u2 * residual_sum / (sigma_e2 * denominator^2)
    beta_sigma_e <-
      2 * sigma_e * as.numeric(crossprod(Xi, inverse_covariance_q))
    information[seq_len(p), p + 1L] <-
      information[seq_len(p), p + 1L] + beta_sigma_u
    information[p + 1L, seq_len(p)] <-
      information[p + 1L, seq_len(p)] + beta_sigma_u
    information[seq_len(p), p + 2L] <-
      information[seq_len(p), p + 2L] + beta_sigma_e
    information[p + 2L, seq_len(p)] <-
      information[p + 2L, seq_len(p)] + beta_sigma_e

    sigma_u_gradient <- residual_sum^2 / denominator^2 -
      Ti / denominator
    sigma_u_second <- sigma_u_gradient + sigma_u * (
      -4 * Ti * sigma_u * residual_sum^2 / denominator^3 +
        2 * Ti^2 * sigma_u / denominator^2
    )
    sigma_cross_second <- sigma_u * (
      -4 * sigma_e * residual_sum^2 / denominator^3 +
        2 * Ti * sigma_e / denominator^2
    )
    inverse_q_quadratic <- q_squared / sigma_e2 -
      sigma_u2 * residual_sum^2 /
        (sigma_e2 * denominator^3)
    sigma_e_gradient <- q_squared -
      (Ti - 1) / sigma_e2 - 1 / denominator
    sigma_e_second <- sigma_e_gradient + sigma_e * (
      -4 * sigma_e * inverse_q_quadratic +
        2 * (Ti - 1) / sigma_e^3 +
        2 * sigma_e / denominator^2
    )

    information[p + 1L, p + 1L] <-
      information[p + 1L, p + 1L] - sigma_u_second
    information[p + 1L, p + 2L] <-
      information[p + 1L, p + 2L] - sigma_cross_second
    information[p + 2L, p + 1L] <-
      information[p + 2L, p + 1L] - sigma_cross_second
    information[p + 2L, p + 2L] <-
      information[p + 2L, p + 2L] - sigma_e_second
  }

  information <- (information + t(information)) / 2
  B <- n * .suest_information_inverse(information)
  dimnames(B) <- dimnames(information)

  if (any(!is.finite(U)) || any(!is.finite(B)))
    stop(
      "Unable to construct finite random-intercept ML scores.",
      call. = FALSE
    )

  list(score = U, bread = B, parameters = parameters)
}

.suest_biprobit_components <- function(model) {
  parameters <- .suest_extract_parameters(
    model,
    "biprobit",
    "mvProbit::mvProbit"
  )
  raw <- stats::coef(model)
  rho <- unname(raw[length(raw)])
  if (!is.finite(rho) || abs(rho) >= 1)
    stop(
      "The bivariate-probit correlation must be strictly between -1 and 1.",
      call. = FALSE
    )

  U <- sandwich::estfun(model)
  hessian <- model$hessian
  if (!is.matrix(U) || !is.matrix(hessian) ||
      ncol(U) != length(parameters) ||
      any(dim(hessian) != length(parameters)))
    stop(
      "The bivariate-probit scores and Hessian do not align.",
      call. = FALSE
    )

  correlation_derivative <- 1 - rho^2
  U[, ncol(U)] <- U[, ncol(U)] * correlation_derivative
  # Compute curvature from analytic bivariate scores in athrho coordinates.
  # The fitting engine's finite-difference Hessian can lose useful precision.
  frame <- .suest_model_frame(model, "mvProbit::mvProbit")
  X <- stats::model.matrix(stats::delete.response(stats::terms(stats::formula(model))), frame)
  signs <- 2*stats::model.response(frame) - 1
  k <- ncol(X)
  score_sum <- function(b) {
    eta1 <- as.numeric(X %*% b[seq_len(k)])
    eta2 <- as.numeric(X %*% b[k + seq_len(k)])
    rho <- tanh(b[length(b)])
    scale <- sqrt(1 - rho^2)
    probability <- vapply(seq_along(eta1), function(i) {
      r <- prod(signs[i, ])*rho
      as.numeric(mvtnorm::pmvnorm(upper = signs[i, ]*c(eta1[i], eta2[i]),
        corr = matrix(c(1, r, r, 1), 2), algorithm = mvtnorm::TVPACK()))
    }, numeric(1))
    first <- signs[, 1]*stats::dnorm(eta1)*
      stats::pnorm(signs[, 2]*(eta2 - rho*eta1)/scale)/probability
    second <- signs[, 2]*stats::dnorm(eta2)*
      stats::pnorm(signs[, 1]*(eta1 - rho*eta2)/scale)/probability
    density <- exp(-(eta1^2 - 2*rho*eta1*eta2 + eta2^2)/(2*scale^2))/(2*pi*scale)
    colSums(cbind(X*first, X*second,
      signs[, 1]*signs[, 2]*density*scale^2/probability))
  }
  information <- .suest_score_information(parameters, score_sum)
  B <- nrow(U)*.suest_information_inverse(information)
  colnames(U) <- names(parameters)
  dimnames(B) <- list(names(parameters), names(parameters))

  if (any(!is.finite(U)) || any(!is.finite(B)))
    stop(
      "Unable to construct finite bivariate-probit scores.",
      call. = FALSE
    )

  list(score = U, bread = B, parameters = parameters)
}

.suest_ivprobit_components <- function(model) {
  parameters <- .suest_extract_parameters(
    model,
    "ivprobit",
    "Rchoice::ivpml"
  )
  U <- sandwich::estfun(model)
  # Rchoice's analytic Hessian can disagree with derivatives of its scores
  # in overidentified IV-probit fits. Recompute curvature at the fitted point.
  X <- stats::model.matrix(model$formula, model$mf, rhs = 1)
  Z <- stats::model.matrix(model$formula, model$mf, rhs = 2)
  y <- stats::model.response(model$mf)
  endogenous <- X[, !colnames(X) %in% colnames(Z), drop = FALSE]
  score_sum <- function(b) {
    k <- ncol(X)
    p <- ncol(Z)
    sigma <- exp(b[k + p + 1L])
    rho <- tanh(b[k + p + 2L])
    scale <- sqrt(1 - rho^2)
    u <- (as.numeric(endogenous) - as.numeric(Z %*% b[k + seq_len(p)]))/sigma
    eta <- as.numeric(X %*% b[seq_len(k)])
    index <- (eta + rho*u)/scale
    sign <- 2*y - 1
    mills <- sign * exp(stats::dnorm(index, log = TRUE) -
      stats::pnorm(sign*index, log.p = TRUE))
    colSums(cbind(X * (mills/scale), Z * (u - mills*rho/scale)/sigma,
      u^2 - 1 - mills*rho*u/scale, mills*(rho*eta + u)/scale))
  }
  B <- nrow(U) * .suest_information_inverse(.suest_score_information(parameters, score_sum))
  if (!is.matrix(U) || !is.matrix(B) ||
      ncol(U) != length(parameters) ||
      any(dim(B) != length(parameters)))
    stop(
      "The instrumental-variable probit scores and bread do not align.",
      call. = FALSE
    )
  colnames(U) <- names(parameters)
  dimnames(B) <- list(names(parameters), names(parameters))
  if (any(!is.finite(U)) || any(!is.finite(B)))
    stop(
      "Unable to construct finite instrumental-variable probit scores.",
      call. = FALSE
    )
  list(score = U, bread = B, parameters = parameters)
}


.suest_poisson_pweight_components <- function(model) {
  mf <- stats::model.frame(model)
  X <- stats::model.matrix(model)
  y <- as.numeric(stats::model.response(mf))
  weights <- .suest_model_weights(mf)
  parameters <- stats::coef(model)
  offset <- stats::model.offset(mf)
  if (is.null(offset))
    offset <- 0
  eta <- as.numeric(X %*% parameters + offset)
  mu <- exp(eta)

  if (any(!is.finite(mu)) || any(mu <= 0))
    stop(
      "Unable to construct finite pweighted Poisson-model means.",
      call. = FALSE
    )

  U <- X * (weights * (y - mu))
  colnames(U) <- names(parameters)
  rownames(U) <- rownames(X)

  A <- crossprod(X, X * (weights * mu))
  A <- (A + t(A)) / 2
  dimnames(A) <- list(names(parameters), names(parameters))

  n <- nrow(X)
  B <- n * .suest_information_inverse(A)
  dimnames(B) <- dimnames(A)

  list(
    score = U,
    bread = B,
    parameters = parameters
  )
}


.suest_binary_pweight_components <- function(model, type) {
  mf <- stats::model.frame(model)
  X <- stats::model.matrix(model)
  y <- stats::model.response(mf)
  y <- if (is.factor(y)) as.integer(y) - 1L else as.numeric(y)
  weights <- .suest_model_weights(mf)
  parameters <- stats::coef(model)
  offset <- stats::model.offset(mf)
  if (is.null(offset))
    offset <- 0
  eta <- as.numeric(X %*% parameters + offset)

  if (type == "logit") {
    probability <- stats::plogis(eta)
    variance <- probability * (1 - probability)
    if (any(!is.finite(probability)) || any(variance <= 0))
      stop(
        paste0(
          "Pweighted binary models require fitted probabilities strictly ",
          "between zero and one."
        ),
        call. = FALSE
      )
    score_eta <- y - probability
    information_eta <- variance
  } else if (type == "probit") {
    probability <- stats::pnorm(eta)
    density <- stats::dnorm(eta)
    variance <- probability * (1 - probability)
    if (any(!is.finite(probability)) || any(variance <= 0))
      stop(
        paste0(
          "Pweighted binary models require fitted probabilities strictly ",
          "between zero and one."
        ),
        call. = FALSE
      )
    residual <- y - probability
    score_eta <- residual * density / variance
    information_eta <- density * (
      ((eta * residual + density) * variance +
         residual * density * (1 - 2 * probability)) /
        variance^2
    )
  } else {
    stop("Internal error: unsupported pweighted binary model.", call. = FALSE)
  }

  if (any(!is.finite(score_eta)) || any(!is.finite(information_eta)))
    stop(
      "Unable to construct finite pweighted binary-model scores.",
      call. = FALSE
    )

  U <- X * (weights * score_eta)
  colnames(U) <- names(parameters)
  rownames(U) <- rownames(X)

  A <- crossprod(X, X * (weights * information_eta))
  A <- (A + t(A)) / 2
  dimnames(A) <- list(names(parameters), names(parameters))

  n <- nrow(X)
  B <- n * .suest_information_inverse(A)
  dimnames(B) <- dimnames(A)

  list(
    score = U,
    bread = B,
    parameters = parameters
  )
}


.suest_information_inverse <- function(A) {
  A <- (A + t(A)) / 2
  scale <- sqrt(pmax(diag(A), .Machine$double.eps))
  standardized <- A / outer(scale, scale)
  inverse <- solve(standardized) / outer(scale, scale)
  (inverse + t(inverse)) / 2
}

.suest_polr_components <- function(model) {
  mf <- stats::model.frame(model)
  X <- stats::model.matrix(model)[, -1L, drop = FALSE]
  y <- as.integer(stats::model.response(mf))
  weights <- stats::model.weights(mf)
  if (is.null(weights))
    weights <- rep.int(1, nrow(X))

  beta <- model$coefficients
  zeta <- model$zeta
  offset <- stats::model.offset(mf)
  if (is.null(offset))
    offset <- 0
  eta <- as.vector(X %*% beta + offset)
  thresholds <- c(-Inf, zeta, Inf)

  lower <- thresholds[y] - eta
  upper <- thresholds[y + 1L] - eta

  if (identical(model$method, "logistic")) {
    F_lower <- stats::plogis(lower)
    F_upper <- stats::plogis(upper)
    f_lower <- stats::dlogis(lower)
    f_upper <- stats::dlogis(upper)
    fp_lower <- f_lower * (1 - 2 * F_lower)
    fp_upper <- f_upper * (1 - 2 * F_upper)
  } else if (identical(model$method, "probit")) {
    F_lower <- stats::pnorm(lower)
    F_upper <- stats::pnorm(upper)
    f_lower <- stats::dnorm(lower)
    f_upper <- stats::dnorm(upper)
    fp_lower <- ifelse(is.finite(lower), -lower * f_lower, 0)
    fp_upper <- ifelse(is.finite(upper), -upper * f_upper, 0)
  } else {
    stop("Unsupported ordered-response link.", call. = FALSE)
  }

  probability <- F_upper - F_lower
  if (any(!is.finite(probability)) || any(probability <= 0))
    stop(
      "Ordered-model fitted probabilities must be positive and finite.",
      call. = FALSE
    )

  d_eta <- f_lower - f_upper
  score_eta <- d_eta / probability
  U_beta <- X * (weights * score_eta)

  m <- length(zeta)
  U_zeta <- matrix(
    0,
    nrow = nrow(X),
    ncol = m,
    dimnames = list(rownames(X), names(zeta))
  )

  upper_exists <- y <= m
  lower_exists <- y > 1L

  if (any(upper_exists)) {
    U_zeta[
      cbind(which(upper_exists), y[upper_exists])
    ] <- weights[upper_exists] *
      f_upper[upper_exists] /
      probability[upper_exists]
  }

  if (any(lower_exists)) {
    U_zeta[
      cbind(which(lower_exists), y[lower_exists] - 1L)
    ] <- -weights[lower_exists] *
      f_lower[lower_exists] /
      probability[lower_exists]
  }

  U <- cbind(U_beta, U_zeta)
  colnames(U_beta) <- colnames(X)
  colnames(U) <- c(colnames(X), names(zeta))
  rownames(U) <- rownames(X)

  p_eta_eta <- fp_upper - fp_lower
  loglik_eta_eta <- p_eta_eta / probability - score_eta^2
  A_beta_beta <- crossprod(
    X,
    X * (-weights * loglik_eta_eta)
  )

  A_beta_zeta <- matrix(
    0,
    nrow = ncol(X),
    ncol = m,
    dimnames = list(colnames(X), names(zeta))
  )
  A_zeta_zeta <- matrix(
    0,
    nrow = m,
    ncol = m,
    dimnames = list(names(zeta), names(zeta))
  )

  for (threshold in seq_len(m)) {
    is_upper <- y == threshold
    is_lower <- y == threshold + 1L

    derivative <- numeric(nrow(X))

    if (any(is_upper)) {
      derivative[is_upper] <- (
        -fp_upper[is_upper] * probability[is_upper] -
          d_eta[is_upper] * f_upper[is_upper]
      ) / probability[is_upper]^2

      loglik_upper_upper <- (
        fp_upper[is_upper] / probability[is_upper] -
          (f_upper[is_upper] / probability[is_upper])^2
      )

      A_zeta_zeta[threshold, threshold] <-
        A_zeta_zeta[threshold, threshold] -
        sum(weights[is_upper] * loglik_upper_upper)
    }

    if (any(is_lower)) {
      derivative[is_lower] <- (
        fp_lower[is_lower] * probability[is_lower] +
          d_eta[is_lower] * f_lower[is_lower]
      ) / probability[is_lower]^2

      loglik_lower_lower <- (
        -fp_lower[is_lower] / probability[is_lower] -
          (f_lower[is_lower] / probability[is_lower])^2
      )

      A_zeta_zeta[threshold, threshold] <-
        A_zeta_zeta[threshold, threshold] -
        sum(weights[is_lower] * loglik_lower_lower)
    }

    A_beta_zeta[, threshold] <- -colSums(
      X * (weights * derivative)
    )
  }

  if (m >= 2L) {
    for (category in 2L:m) {
      observations <- y == category
      if (any(observations)) {
        cross_second <- (
          f_lower[observations] *
            f_upper[observations] /
            probability[observations]^2
        )
        value <- -sum(weights[observations] * cross_second)
        A_zeta_zeta[category - 1L, category] <- value
        A_zeta_zeta[category, category - 1L] <- value
      }
    }
  }

  A <- rbind(
    cbind(A_beta_beta, A_beta_zeta),
    cbind(t(A_beta_zeta), A_zeta_zeta)
  )
  dimnames(A) <- list(colnames(U), colnames(U))

  n <- nrow(X)
  B <- n * .suest_information_inverse(A)

  list(
    score = U,
    bread = B,
    parameters = .suest_extract_parameters(
      model,
      if (identical(model$method, "logistic")) "ologit" else "oprobit",
      "MASS::polr"
    )
  )
}

.suest_multinom_components <- function(model) {
  mf <- stats::model.frame(model)
  X <- stats::model.matrix(model)
  weights <- .suest_model_weights(mf, model, "nnet::multinom")

  probability <- stats::fitted(model)
  if (is.null(dim(probability)))
    probability <- cbind(1 - probability, probability)

  residual <- stats::residuals(model)
  if (is.null(dim(residual)))
    residual <- cbind(-residual, residual)

  categories <- model$lev
  nonreference <- seq.int(2L, length(categories))
  k <- ncol(X)
  q <- length(nonreference)
  n <- nrow(X)

  parameter_names <- as.vector(t(outer(
    categories[-1L],
    colnames(X),
    paste,
    sep = ":"
  )))

  U <- matrix(
    0,
    nrow = n,
    ncol = q * k,
    dimnames = list(rownames(X), parameter_names)
  )

  for (outcome in seq_len(q)) {
    columns <- (outcome - 1L) * k + seq_len(k)
    U[, columns] <- X *
      (weights * residual[, nonreference[outcome]])
  }

  A <- matrix(
    0,
    nrow = q * k,
    ncol = q * k,
    dimnames = list(parameter_names, parameter_names)
  )

  for (outcome1 in seq_len(q)) {
    p1 <- probability[, nonreference[outcome1]]
    columns1 <- (outcome1 - 1L) * k + seq_len(k)

    for (outcome2 in seq_len(q)) {
      p2 <- probability[, nonreference[outcome2]]
      columns2 <- (outcome2 - 1L) * k + seq_len(k)

      weight <- if (outcome1 == outcome2) {
        p1 * (1 - p1)
      } else {
        -p1 * p2
      }

      A[columns1, columns2] <- crossprod(
        X,
        X * (weights * weight)
      )
    }
  }

  B <- n * .suest_information_inverse(A)

  list(
    score = U,
    bread = B,
    parameters = .suest_extract_parameters(
      model,
      "multinom",
      "nnet::multinom"
    )
  )
}


.suest_negbin_components <- function(model) {
  mf <- stats::model.frame(model)
  X <- stats::model.matrix(model)
  y <- stats::model.response(mf)
  weights <- .suest_model_weights(mf)
  mu <- stats::fitted(model)
  theta <- as.numeric(model$theta)
  n <- nrow(X)

  beta_score_factor <- weights * theta * (y - mu) / (theta + mu)
  U_beta <- X * beta_score_factor

  g <- digamma(y + theta) - digamma(theta) +
    log(theta) + 1 - log(theta + mu) -
    (theta + y) / (theta + mu)
  U_theta <- weights * theta * g

  U <- cbind(U_beta, ln_theta = U_theta)
  colnames(U_beta) <- colnames(X)
  colnames(U) <- c(colnames(X), "ln_theta")
  rownames(U) <- rownames(X)

  weight_beta <- weights * theta * mu * (theta + y) / (theta + mu)^2
  A_bb <- crossprod(X, X * weight_beta)

  h_beta_theta <- weights * theta * mu * (y - mu) / (theta + mu)^2
  A_btheta <- -colSums(X * h_beta_theta)

  g_prime <- trigamma(y + theta) - trigamma(theta) +
    1 / theta - 1 / (theta + mu) +
    (y - mu) / (theta + mu)^2
  A_thetatheta <- -sum(weights * (theta * g + theta^2 * g_prime))

  A <- rbind(
    cbind(A_bb, A_btheta),
    c(A_btheta, A_thetatheta)
  )
  A <- (A + t(A)) / 2
  dimnames(A) <- list(colnames(U), colnames(U))

  B <- n * solve(A)
  B <- (B + t(B)) / 2

  list(
    score = U,
    bread = B,
    parameters = .suest_extract_parameters(
      model,
      "negbin",
      "MASS::glm.nb"
    )
  )
}

.suest_zinb_components <- function(model) {
  X <- model$x$count
  Z <- model$x$zero
  y <- as.numeric(model$y)
  weights <- as.numeric(model$weights)
  beta <- model$coefficients$count
  gamma <- model$coefficients$zero
  theta <- as.numeric(model$theta)

  count_offset <- if (is.null(model$offset$count)) 0 else model$offset$count
  zero_offset <- if (is.null(model$offset$zero)) 0 else model$offset$zero

  mu <- exp(as.numeric(X %*% beta) + count_offset)
  zero_eta <- as.numeric(Z %*% gamma) + zero_offset
  zero_link <- stats::make.link(model$link)
  probability_zero <- zero_link$linkinv(zero_eta)
  probability_derivative <- zero_link$mu.eta(zero_eta)
  probability_nb_zero <- stats::dnbinom(0, size = theta, mu = mu)

  is_zero <- y == 0
  mixture_probability <- probability_zero +
    (1 - probability_zero)*probability_nb_zero
  posterior_count <- rep.int(1, length(y))
  posterior_count[is_zero] <-
    (1 - probability_zero[is_zero])*probability_nb_zero[is_zero] /
    mixture_probability[is_zero]

  score_count_eta <- theta*(y - mu)/(theta + mu)
  score_count_eta[is_zero] <-
    posterior_count[is_zero]*score_count_eta[is_zero]

  score_zero_eta <- -probability_derivative/(1 - probability_zero)
  score_zero_eta[is_zero] <-
    probability_derivative[is_zero]*
    (1 - probability_nb_zero[is_zero]) /
    mixture_probability[is_zero]

  theta_score <- theta*(
    digamma(y + theta) - digamma(theta) + log(theta) + 1 -
      log(theta + mu) - (theta + y)/(theta + mu)
  )
  theta_score[is_zero] <-
    posterior_count[is_zero]*theta_score[is_zero]

  U <- cbind(
    X*as.numeric(weights*score_count_eta),
    Z*as.numeric(weights*score_zero_eta),
    ln_theta = weights*theta_score
  )
  parameters <- .suest_extract_parameters(
    model,
    "zinb",
    "pscl::zeroinfl"
  )
  colnames(U) <- names(parameters)

  H <- model$optim$hessian
  if (!is.matrix(H) || any(dim(H) != length(parameters)))
    stop(
      "Unable to retrieve the full observed-information matrix for the ZINB model.",
      call. = FALSE
    )

  information <- -(H + t(H))/2
  B <- nrow(U)*solve(information)
  B <- (B + t(B))/2
  dimnames(B) <- list(names(parameters), names(parameters))

  if (any(!is.finite(U)) || any(!is.finite(B)))
    stop(
      "Unable to construct finite zero-inflated negative-binomial scores.",
      call. = FALSE
    )

  list(score = U, bread = B, parameters = parameters)
}

.suest_truncreg_components <- function(model) {
  parameters <- .suest_extract_parameters(
    model,
    "truncreg",
    "truncreg::truncreg"
  )
  U <- sandwich::estfun(model)
  B <- sandwich::bread(model)

  if (ncol(U) != length(parameters))
    stop(
      "The truncated-regression score does not match its coefficient vector.",
      call. = FALSE
    )

  colnames(U) <- names(parameters)
  dimnames(B) <- list(names(parameters), names(parameters))
  list(score = U, bread = B, parameters = parameters)
}

.suest_model_components <- function(
    model,
    type,
    engine,
    weight_type = NULL) {
  engine <- .suest_engine_name(engine)
  if (identical(weight_type, "pweight") && type == "lm") {
    .suest_lm_pweight_components(model)
  } else if (type == "panel_fe") {
    .suest_plm_fe_components(model)
  } else if (type == "panel_be") {
    .suest_plm_be_components(model)
  } else if (type == "panel_re") {
    .suest_plm_re_components(model)
  } else if (type == "panel_ml") {
    .suest_lme_ml_components(model)
  } else if (type == "panel_gee") {
    .suest_gee_components(model)
  } else if (type == "biprobit") {
    .suest_biprobit_components(model)
  } else if (type == "ivprobit") {
    .suest_ivprobit_components(model)
  } else if (type == "lm") {
    .suest_lm_components(model)
  } else if (type == "probit" ||
             (identical(weight_type, "pweight") && type == "logit")) {
    .suest_binary_pweight_components(model, type)
  } else if (identical(weight_type, "pweight") && type == "poisson") {
    .suest_poisson_pweight_components(model)
  } else if (type == "negbin") {
    .suest_negbin_components(model)
  } else if (type == "zinb") {
    .suest_zinb_components(model)
  } else if (type == "truncreg") {
    .suest_truncreg_components(model)
  } else if (identical(engine, "ordinal::clm")) {
    list(
      score = sandwich::estfun(model),
      bread = sandwich::bread(model),
      parameters = .suest_extract_parameters(model, type, engine)
    )
  } else if (type %in% c("ologit", "oprobit")) {
    .suest_polr_components(model)
  } else if (type == "multinom") {
    .suest_multinom_components(model)
  } else {
    list(
      score = sandwich::estfun(model),
      bread = sandwich::bread(model),
      parameters = .suest_extract_parameters(model, type, engine)
    )
  }
}


.suest_block_diag <- function(...) {
  matrices <- list(...)
  if (length(matrices) == 1L && is.list(matrices[[1L]]))
    matrices <- matrices[[1L]]
  if (!length(matrices))
    return(matrix(numeric(0), 0L, 0L))

  dimensions <- vapply(matrices, nrow, integer(1))
  out <- matrix(0, sum(dimensions), sum(dimensions))
  starts <- cumsum(c(1L, dimensions[-length(dimensions)]))

  for (i in seq_along(matrices)) {
    index <- starts[i] - 1L + seq_len(dimensions[i])
    out[index, index] <- matrices[[i]]
  }
  out
}

.suest_align_scores <- function(U, keys, union_keys) {
  out <- matrix(
    0,
    nrow = length(union_keys),
    ncol = ncol(U),
    dimnames = list(union_keys, colnames(U))
  )
  out[match(keys, union_keys), ] <- U
  out
}

.suest_pairwise_overlap <- function(sample_keys, model_names) {
  n_models <- length(sample_keys)
  out <- matrix(
    0L,
    nrow = n_models,
    ncol = n_models,
    dimnames = list(model_names, model_names)
  )

  for (i in seq_len(n_models)) {
    for (j in i:n_models) {
      value <- if (i == j) {
        length(sample_keys[[i]])
      } else {
        length(intersect(sample_keys[[i]], sample_keys[[j]]))
      }
      out[i, j] <- value
      out[j, i] <- value
    }
  }
  out
}
