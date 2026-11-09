# Methods used internally by marginaleffects.

get_coef.suest_model <- function(model, ...) {
  model$coefficients
}

set_coef.suest_model <- function(model, coefs, ...) {
  out <- model

  if (is.null(names(coefs))) {
    if (length(coefs) != length(out$coefficients))
      stop("Incorrect number of coefficients.", call. = FALSE)
    names(coefs) <- names(out$coefficients)
  }

  if (!all(names(out$coefficients) %in% names(coefs)))
    stop(
      "The supplied coefficients do not match the SUEST parameters.",
      call. = FALSE
    )

  coefs <- coefs[names(out$coefficients)]
  out$coefficients <- coefs

  for (i in seq_along(out$models)) {
    b <- coefs[out$index[[i]]]
    names(b) <- out$local_names[[i]]
    type <- out$model_types[i]
    engine <- if (!is.null(out$model_engines)) {
      .suest_engine_name(out$model_engines[[i]])
    } else {
      .suest_model_engine(out$models[[i]])
    }

    out$models[[i]] <- .suest_set_parameters(
      out$models[[i]],
      b,
      type,
      engine
    )
  }

  out
}

get_vcov.suest_model <- function(model, vcov = NULL, ...) {
  model$vcov
}

# Data used when newdata is omitted and by datagrid().
get_data.suest_model <- function(x, ...) .suest_default_data(x)

# Each model's outcomes, fixed-effect predictors, and every variable it uses.
.suest_variable_roles <- function(m) {
  if (inherits(m, "suest_gjrm")) {
    variables <- lapply(m$formulas[1:2], all.vars)
    return(list(response = m$responses,
      predictors = unique(unlist(lapply(variables, `[`, -1L))),
      all = unique(unlist(variables))))
  }
  fit <- if (inherits(m, c("suest_lme4", "suest_vglm"))) m$fit else m
  f <- tryCatch(stats::formula(m), error = function(e) m$formula)
  two_sided <- inherits(f, "formula") && length(f) == 3L
  response <- tryCatch(insight::find_response(fit, combine = FALSE), error = function(e) NULL)
  predictors <- tryCatch(insight::find_predictors(fit, effects = "fixed", component = "all"),
    error = function(e) NULL)
  # Instruments, survey design variables, and fixed-effect ids are not predictors.
  predictors <- unique(unlist(predictors[setdiff(names(predictors),
    c("instruments", "endogenous", "design", "cluster"))], use.names = FALSE))
  if (is.null(response) && two_sided) response <- all.vars(f[[2L]])
  if (is.null(predictors) && two_sided) predictors <- all.vars(f[[3L]])
  list(response = response, predictors = predictors,
       all = unique(c(if (inherits(f, "formula")) all.vars(f), all.vars(m$call))))
}

# Predictors (default variables for slopes and comparisons) plus, as auxiliary,
# every other variable the models use, so counterfactual grids (e.g.
# predictions(variables = )) keep these columns and get a fresh rowid. NULL if
# any model's variables are unknown.
find_predictors.suest_model <- function(x, component = "all", flatten = FALSE, ...) {
  roles <- lapply(x$models, .suest_variable_roles)
  if (any(lengths(lapply(roles, `[[`, "all")) == 0L)) return(NULL)
  predictors <- unique(unlist(lapply(roles, `[[`, "predictors")))
  auxiliary <- setdiff(
    unique(c(unlist(lapply(roles, `[[`, "all")), ".suest_model", ".suest_rowid", ".suest_weight")),
    predictors
  )
  out <- list(conditional = predictors, auxiliary = auxiliary)
  if (!"all" %in% component) out <- out[intersect(names(out), component)]
  if (isTRUE(flatten)) unique(unlist(out, use.names = FALSE)) else out
}

# Outcomes that are not a predictor in another model (a mediator stays a
# predictor), so default slopes and comparisons skip them.
find_response.suest_model <- function(x, ...) {
  roles <- lapply(x$models, .suest_variable_roles)
  out <- setdiff(
    unique(unlist(lapply(roles, `[[`, "response"))),
    unlist(lapply(roles, `[[`, "predictors"))
  )
  if (length(out)) out else NULL
}

get_coef.suest_mi <- function(model, ...) model$coefficients

get_vcov.suest_mi <- function(model, vcov = NULL, ...) model$vcov

# marginaleffects estimates within each imputation and pools (suest_mi is a
# "mira"); this is reached only by functions without that support.
get_predict.suest_mi <- function(model, newdata, type = "response", ...) {
  stop(
    paste0(
      "This function does not pool across imputations. Use predictions(), ",
      "comparisons(), slopes(), hypotheses(), or their avg_ versions, which ",
      "estimate within each imputation and pool with Rubin's rules."
    ),
    call. = FALSE
  )
}

.suest_prepare_newdata <- function(model, newdata, model_frame) {
  out <- newdata
  xlevels <- model$xlevels

  if (!is.null(xlevels)) {
    for (nm in names(xlevels)) {
      if (nm %in% names(out)) {
        ordered_predictor <- nm %in% names(model_frame) &&
          is.ordered(model_frame[[nm]])

        values <- as.character(out[[nm]])
        out[[nm]] <- if (ordered_predictor) {
          ordered(values, levels = xlevels[[nm]])
        } else {
          factor(values, levels = xlevels[[nm]])
        }
      }
    }
  }

  out
}

get_predict.suest_model <- function(model, newdata, type = "response", ...) {
  if (is.null(type))
    type <- "response"

  categorical <- model$model_types %in% c("ologit", "oprobit", "multinom", "gologit", "re_ologit", "re_oprobit")
  zero_inflated <- model$model_types %in% c("zip", "zinb")
  zero_inflated_type <- type %in% .suest_zeroinfl_types

  if (zero_inflated_type && !all(zero_inflated))
    stop(
      sprintf(
        paste0(
          "type = \"%s\" is for zero-inflated models (pscl::zeroinfl), and ",
          "every model in the system must be one. Not zero-inflated: %s."
        ),
        type,
        paste(model$model_names[!zero_inflated], collapse = ", ")
      ),
      call. = FALSE
    )

  # A zip and a zinb model are different types, but they share every
  # zero-inflated prediction type.
  if (model$mixed_models && !all(zero_inflated) && !type %in% c("response", "probs"))
    stop(
      "Cross-type SUEST comparisons are supported on the response scale only.",
      call. = FALSE
    )

  # Count probabilities, Pr(y = 0), Pr(y = 1), ..., have one outcome per count,
  # so they are returned the way category probabilities are.
  if (identical(type, "prob"))
    categorical <- categorical | zero_inflated

  model_specific <- ".suest_model" %in% names(newdata)
  # With model-specific rows, marginaleffects can replicate the stacked data
  # many times when it constructs multi-level factor contrasts. For scalar
  # outcomes each input row still produces exactly one prediction, so the
  # returned predictions must remain in the same physical row order as the
  # counterfactual data. Concatenating model 1, then model 2, otherwise
  # misaligns predictions with contrast labels when samples are disjoint.
  preserve_input_order <- model_specific && !any(categorical)
  # Categorical comparisons recycle counterfactual metadata once per outcome.
  # Every outcome block must span the entire input grid; rows belonging to
  # another model are missing predictions, removed before averaging.
  pad_outcome_blocks <- model_specific && any(categorical)

  get_model_data <- function(i) {
    if (model_specific) {
      keep <- as.character(newdata$.suest_model) == model$model_names[i]
      out <- newdata[keep, , drop = FALSE]
    } else {
      out <- newdata
    }

    if (nrow(out) == 0L)
      stop(
        sprintf(
          paste0(
            "No newdata rows were supplied for model '%s'. The .suest_model ",
            "column says which model each row belongs to. If newdata came ",
            "from datagrid(), add .suest_model = c(%s) to the datagrid() call."
          ),
          model$model_names[i],
          paste0('"', model$model_names, '"', collapse = ", ")
        ),
        call. = FALSE
      )

    .suest_prepare_newdata(
      model$models[[i]],
      out,
      model$model_frames[[i]]
    )
  }

  get_rowid <- function(x) {
    # marginaleffects adds an integer rowid used internally for joins.
    # Prefer it over the SUEST bookkeeping column so returned predictions
    # retain the same rowid type as the incoming newdata.
    if ("rowid" %in% names(x)) {
      x$rowid
    } else if (".suest_rowid" %in% names(x)) {
      as.integer(x$.suest_rowid)
    } else {
      seq_len(nrow(x))
    }
  }

  if (any(categorical) && !zero_inflated_type && !type %in% c("response", "probs"))
    stop(
      "Systems containing categorical models require response-scale predictions.",
      call. = FALSE
    )
  if (!any(categorical) && !zero_inflated_type && !type %in% c("response", "link"))
    stop(
      if (all(zero_inflated)) {
        "'type' must be 'response', 'count', 'zero', or 'prob' for zero-inflated models."
      } else {
        "'type' must be 'response' or 'link'."
      },
      call. = FALSE
    )

  category_levels <- lapply(seq_along(model$models), function(i) {
    if (!categorical[i])
      return(NULL)
    if (zero_inflated[i])
      return(.suest_zeroinfl_count_levels(model$models[[i]]))
    if (!is.null(model$category_levels))
      return(model$category_levels[[i]])
    .suest_category_levels(
      model$models[[i]],
      .suest_model_engine(model$models[[i]])
    )
  })
  group_levels <- unlist(lapply(seq_along(model$models), function(i) {
    if (categorical[i]) {
      paste0(model$model_names[i], "::", category_levels[[i]])
    } else {
      model$model_names[i]
    }
  }), use.names = FALSE)

  model_predictions <- function(i) {
    nd <- get_model_data(i)
    rowid <- get_rowid(nd)
    input_order <- if (preserve_input_order) {
      which(as.character(newdata$.suest_model) == model$model_names[i])
    } else {
      NULL
    }
    engine <- if (!is.null(model$model_engines)) {
      .suest_engine_name(model$model_engines[[i]])
    } else {
      .suest_model_engine(model$models[[i]])
    }

    if (categorical[i]) {
      p <- if (zero_inflated[i]) {
        .suest_predict_values(model$models[[i]], nd, "prob", engine, model$model_types[i])
      } else {
        .suest_predict_probabilities(model$models[[i]], nd, engine)
      }
      if (pad_outcome_blocks) {
        keep <- which(as.character(newdata$.suest_model) == model$model_names[i])
        full <- matrix(NA_real_, nrow(newdata), ncol(p),
          dimnames = list(NULL, colnames(p)))
        full[keep, ] <- p
        p <- full
        rowid <- get_rowid(newdata)
      }
      group <- rep(
        paste0(model$model_names[i], "::", colnames(p)),
        each = nrow(p)
      )
      out <- data.frame(
        rowid = rep(rowid, times = ncol(p)),
        group = factor(group, levels = group_levels),
        estimate = as.vector(p),
        check.names = FALSE
      )
    } else {
      prediction_type <- if (identical(type, "probs")) "response" else type
      p <- .suest_predict_values(
        model$models[[i]],
        nd,
        prediction_type,
        engine,
        model$model_types[i]
      )
      if (pad_outcome_blocks) {
        keep <- which(as.character(newdata$.suest_model) == model$model_names[i])
        full <- rep(NA_real_, nrow(newdata))
        full[keep] <- as.numeric(p)
        p <- full
        rowid <- get_rowid(newdata)
      }
      out <- data.frame(
        rowid = rowid,
        group = factor(model$model_names[i], levels = group_levels),
        estimate = as.numeric(p),
        check.names = FALSE
      )
    }

    if (preserve_input_order)
      out$.suest_input_order <- input_order
    out
  }

  out <- do.call(
    rbind,
    lapply(seq_along(model$models), model_predictions)
  )

  if (preserve_input_order) {
    out <- out[order(out$.suest_input_order), , drop = FALSE]
    out$.suest_input_order <- NULL
    rownames(out) <- NULL
  }

  out
}
