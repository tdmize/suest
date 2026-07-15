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

    if (type %in% c("ologit", "oprobit")) {
      k <- length(out$models[[i]]$coefficients)
      q <- length(out$models[[i]]$zeta)
      out$models[[i]]$coefficients <- b[seq_len(k)]
      out$models[[i]]$zeta <- b[k + seq_len(q)]
    } else if (type == "multinom") {
      m <- out$models[[i]]
      r <- length(m$vcoefnames)
      q <- length(m$lev) - 1L

      cf <- matrix(
        b,
        nrow = q,
        ncol = r,
        byrow = TRUE,
        dimnames = list(m$lev[-1L], m$vcoefnames)
      )

      W <- matrix(m$wts, nrow = m$n[3L], byrow = TRUE)
      W[-1L, 1L + seq_len(r)] <- cf
      m$wts <- as.vector(t(W))
      out$models[[i]] <- m
    } else if (type == "negbin") {
      m <- out$models[[i]]
      beta_names <- setdiff(names(b), "ln_theta")
      m$coefficients <- b[beta_names]
      m$theta <- exp(unname(b["ln_theta"]))
      m$family <- MASS::negative.binomial(
        theta = m$theta,
        link = m$family$link
      )
      out$models[[i]] <- m
    } else {
      out$models[[i]]$coefficients <- b
    }
  }

  out
}

get_vcov.suest_model <- function(model, vcov = NULL, ...) {
  model$vcov
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

  if (any(categorical) && !all(categorical))
    stop(
      "Internal error: categorical and scalar models cannot be predicted together.",
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

  if (all(categorical)) {
    if (!type %in% c("response", "probs"))
      stop(
        "Categorical models require type = 'response' or type = 'probs'.",
        call. = FALSE
      )

    category_levels <- model$models[[1]]$lev
    group_levels <- unlist(lapply(
      seq_along(model$models),
      function(i) paste0(model$model_names[i], "::", category_levels)
    ), use.names = FALSE)

    categorical_predictions <- function(i) {
      nd <- get_model_data(i)
      rowid <- get_rowid(nd)
      p <- stats::predict(
        model$models[[i]],
        newdata = nd,
        type = "probs"
      )

      if (is.null(dim(p))) {
        category <- names(p)
        p <- matrix(
          as.numeric(p),
          nrow = 1L,
          dimnames = list(NULL, category)
        )
      }

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
    }

    return(rbind(
      categorical_predictions(1L),
      categorical_predictions(2L)
    ))
  }

  if (!type %in% c("response", "link"))
    stop("'type' must be 'response' or 'link'.", call. = FALSE)

  scalar_predictions <- function(i) {
    nd <- get_model_data(i)
    rowid <- get_rowid(nd)
    predict_type <- if (model$model_types[i] == "lm") "response" else type

    p <- stats::predict(
      model$models[[i]],
      newdata = nd,
      type = predict_type
    )

    data.frame(
      rowid = rowid,
      group = model$model_names[i],
      estimate = as.numeric(p),
      check.names = FALSE
    )
  }

  rbind(
    scalar_predictions(1L),
    scalar_predictions(2L)
  )
}
