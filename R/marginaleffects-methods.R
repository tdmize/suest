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

get_coef.suest_mi <- function(model, ...) model$coefficients

get_vcov.suest_mi <- function(model, vcov = NULL, ...) model$vcov

get_predict.suest_mi <- function(model, newdata, type = "response", ...) {
  stop(
    paste0(
      "Predictions, slopes, and comparisons are not yet implemented for ",
      "multiply imputed SUEST systems. Estimate the target estimand within ",
      "each imputation and pool it using Rubin's rules."
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

  categorical <- model$model_types %in% c("ologit", "oprobit", "multinom")

  if (model$mixed_models && !type %in% c("response", "probs"))
    stop(
      "Cross-type SUEST comparisons are supported on the response scale only.",
      call. = FALSE
    )

  model_specific <- ".suest_model" %in% names(newdata)

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
          "No newdata rows were supplied for model '%s'.",
          model$model_names[i]
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

  if (any(categorical) && !type %in% c("response", "probs"))
    stop(
      "Systems containing categorical models require response-scale predictions.",
      call. = FALSE
    )
  if (!any(categorical) && !type %in% c("response", "link"))
    stop("'type' must be 'response' or 'link'.", call. = FALSE)

  category_levels <- lapply(seq_along(model$models), function(i) {
    if (!categorical[i])
      return(NULL)
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
    engine <- if (!is.null(model$model_engines)) {
      .suest_engine_name(model$model_engines[[i]])
    } else {
      .suest_model_engine(model$models[[i]])
    }

    if (categorical[i]) {
      p <- .suest_predict_probabilities(model$models[[i]], nd, engine)
      group <- rep(
        paste0(model$model_names[i], "::", colnames(p)),
        each = nrow(p)
      )
      data.frame(
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
      data.frame(
        rowid = rowid,
        group = factor(model$model_names[i], levels = group_levels),
        estimate = as.numeric(p),
        check.names = FALSE
      )
    }
  }

  do.call(
    rbind,
    lapply(seq_along(model$models), model_predictions)
  )
}
