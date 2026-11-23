# Two-part models: a binary first part (is the outcome positive?) and a second
# part for the outcome among the positive outcomes, combined as
# E(y) = Pr(y > 0) * E(y | y > 0).

.suest_twopart_first_types <- c("logit", "probit", "cloglog", "fraclogit",
  "fracprobit", "fraccloglog")
.suest_twopart_second_types <- c("lm", "survey_lm", "glm", "poisson", "negbin")
.suest_twopart_types <- c("response", "positive", "amount")

.suest_twopart_setup <- function(fit, twopart, retransform) {
  if (is.null(twopart) || isFALSE(twopart)) {
    if (!is.null(retransform))
      stop("'retransform' applies to two-part models; specify 'twopart' too.", call. = FALSE)
    return(fit)
  }
  nm <- fit$model_names
  in_pairs <- function(x) {
    if (!length(x) || length(x) %% 2L)
      stop("'twopart' must name the models in pairs: each first part, then its second part.",
        call. = FALSE)
    unname(split(x, rep(seq_len(length(x) / 2L), each = 2L)))
  }
  pairs <- if (isTRUE(twopart)) {
    in_pairs(nm)
  } else if (is.character(twopart)) {
    in_pairs(twopart)
  } else if (is.list(twopart) &&
             all(vapply(twopart, function(p) is.character(p) && length(p) == 2L, logical(1)))) {
    twopart
  } else {
    stop(paste0("'twopart' must be TRUE, the model names in pairs, or a list of pairs, ",
      "each naming a first part and its second part."), call. = FALSE)
  }
  unknown <- setdiff(unlist(pairs), nm)
  if (length(unknown))
    stop(sprintf("'twopart' names models not in the system: %s. The models are %s.",
      paste(unknown, collapse = ", "), paste(nm, collapse = ", ")), call. = FALSE)
  firsts <- vapply(pairs, `[`, "", 1L)
  seconds <- vapply(pairs, `[`, "", 2L)
  if (length(intersect(firsts, seconds)))
    stop("A model cannot be both a first part and a second part.", call. = FALSE)
  labels <- names(pairs)
  if (is.null(labels)) labels <- rep("", length(pairs))
  labels[labels == ""] <- paste(firsts, seconds, sep = "+")[labels == ""]
  if (!is.null(retransform))
    retransform <- match.arg(retransform, c("duan", "normal"))

  types <- fit$model_types
  units <- lapply(seq_along(pairs), function(k) {
    f <- match(firsts[k], nm)
    s <- match(seconds[k], nm)
    if (!types[f] %in% .suest_twopart_first_types)
      stop(sprintf(paste0("The first part of '%s' must be a binary logit, probit, or ",
        "complementary log-log model."), labels[k]), call. = FALSE)
    if (!types[s] %in% .suest_twopart_second_types)
      stop(sprintf(paste0("The second part of '%s' must be a linear model, a GLM such as ",
        "a Gamma or (quasi-)Poisson model with a log link, or a negative binomial model."),
        labels[k]), call. = FALSE)
    y <- fit$models[[f]]$y
    if (is.null(y)) y <- stats::model.response(fit$model_frames[[f]])
    y <- as.numeric(y)
    keys <- fit$sample_keys
    if (!all(y %in% c(0, 1)) || length(y) != length(keys[[f]]) ||
        !setequal(keys[[s]], keys[[f]][y == 1]))
      stop(sprintf(paste0("The second part of '%s' must be fit to exactly the first part's ",
        "observations with a positive outcome (for example with subset = y > 0), with the ",
        "first part's outcome 1 for a positive outcome. If the parts use different data ",
        "objects, identify the observations with observation_id."), labels[k]), call. = FALSE)
    log_outcome <- !is.null(retransform) && types[s] %in% c("lm", "survey_lm")
    if (log_outcome && types[s] == "survey_lm")
      stop("'retransform' is not available for a survey linear second part.", call. = FALSE)
    list(kind = "pair", label = labels[k], first = f, second = s,
      retransform = if (log_outcome) retransform)
  })
  if (!is.null(retransform) && all(vapply(units, function(u) is.null(u$retransform), logical(1))))
    stop("'retransform' applies to lm() second parts, regressions of the log of the outcome.",
      call. = FALSE)
  others <- setdiff(seq_along(nm), c(match(firsts, nm), match(seconds, nm)))
  categorical <- types[others] %in% c("ologit", "oprobit", "multinom", "gologit",
    "re_ologit", "re_oprobit")
  if (any(categorical))
    stop(sprintf("Models beside two-part models must have a single outcome: %s does not.",
      paste(nm[others][categorical], collapse = ", ")), call. = FALSE)
  units <- c(units, lapply(others, function(i) list(kind = "model", label = nm[i], model = i)))
  if (anyDuplicated(vapply(units, `[[`, "", "label")))
    stop("Each two-part model needs a label that differs from the model names.", call. = FALSE)

  # The retransformation factor is held fixed at its estimate, as in suest2 and
  # twopm: Duan's mean of exp(residual), or exp(sigma^2/2) with sigma^2 the
  # mean squared residual.
  units <- lapply(units, function(u) {
    if (is.null(u$retransform)) return(u)
    r <- stats::residuals(fit$models[[u$second]])
    w <- if (is.null(fit$weight_type)) rep(1, length(r)) else fit$model_weights[[u$second]]
    u$factor <- if (u$retransform == "duan") stats::weighted.mean(exp(r), w) else
      exp(stats::weighted.mean(r^2, w) / 2)
    u
  })
  fit$twopart <- list(units = units, retransform = retransform)

  p <- get_predict.suest_model(fit, .suest_twopart_data(fit, stack = TRUE))
  if (anyNA(p$estimate))
    stop(paste0("Some observations in a first part's sample lack a predictor of the ",
      "second part. Fit both parts to the same observations."), call. = FALSE)
  fit
}

# Default data: each two-part model's first-part sample, with every variable of
# the first part's data, and each other model's own sample. When every unit has
# the same sample, one row per observation.
.suest_twopart_data <- function(x, stack = FALSE) {
  units <- x$twopart$units
  frames <- .suest_stacked_frames(x)
  # Columns missing from a block are taken from another block of the same people.
  fill <- function(out, other, keys, other_keys) {
    at <- match(keys, other_keys)
    for (nm in names(other))
      if (!nm %in% names(out) || (all(is.na(out[[nm]])) && !all(is.na(other[[nm]]))))
        out[[nm]] <- other[[nm]][at]
    out
  }
  model_block <- function(i) {
    out <- frames[as.character(frames$.suest_model) == x$model_names[i], , drop = FALSE]
    out$.suest_model <- out$.suest_rowid <- NULL
    out
  }
  blocks <- lapply(units, function(u) {
    if (u$kind == "model") return(model_block(u$model))
    out <- model_block(u$first)
    data <- .suest_twopart_source(x, u$first)
    if (!is.null(data)) {
      index <- match(rownames(x$model_frames[[u$first]]), rownames(data))
      if (!anyNA(index)) out[names(data)] <- as.data.frame(data)[index, , drop = FALSE]
    }
    fill(out, model_block(u$second), x$sample_keys[[u$first]], x$sample_keys[[u$second]])
  })
  keys <- lapply(units, function(u) x$sample_keys[[if (u$kind == "pair") u$first else u$model]])
  if (!stack && all(vapply(keys[-1L], setequal, logical(1), keys[[1L]]))) {
    out <- blocks[[1L]]
    for (k in seq_along(blocks)[-1L]) out <- fill(out, blocks[[k]], keys[[1L]], keys[[k]])
    rownames(out) <- NULL
    return(out)
  }
  labels <- vapply(units, `[[`, "", "label")
  out <- .suest_rbind_blocks(Map(function(b, l) {
    b$.suest_model <- l
    b
  }, blocks, labels))
  out$.suest_rowid <- seq_len(nrow(out))
  class(out) <- c("suest_newdata", "data.frame")
  out
}

.suest_twopart_source <- function(x, i) {
  m <- x$models[[i]]
  if (identical(x$weight_type, "survey") || is.null(m$call$data)) return(NULL)
  data <- tryCatch(eval(m$call$data, environment(stats::formula(m))), error = function(e) NULL)
  if (is.data.frame(data)) data else NULL
}

# Row-bind data frames whose factors may have different levels.
.suest_rbind_blocks <- function(blocks) {
  all_names <- unique(unlist(lapply(blocks, names)))
  specs <- list()
  for (nm in all_names) {
    present <- Filter(Negate(is.null), lapply(blocks, `[[`, nm))
    if (any(vapply(present, is.factor, logical(1)))) {
      ordered <- all(vapply(present, is.ordered, logical(1)))
      specs[[nm]] <- list(levels = unique(unlist(lapply(present, function(v)
        if (is.factor(v)) levels(v) else unique(as.character(v[!is.na(v)]))))),
        ordered = ordered)
    }
  }
  blocks <- lapply(blocks, function(b) {
    for (nm in names(specs)) if (nm %in% names(b)) b[[nm]] <- as.character(b[[nm]])
    for (nm in setdiff(all_names, names(b))) b[[nm]] <- NA
    b[all_names]
  })
  out <- do.call(rbind, blocks)
  rownames(out) <- NULL
  for (nm in names(specs))
    out[[nm]] <- factor(out[[nm]], levels = specs[[nm]]$levels, ordered = specs[[nm]]$ordered)
  out
}

.suest_twopart_values <- function(model, i, nd, type = "response") {
  engine <- .suest_engine_name(model$model_engines[[i]])
  nd <- .suest_prepare_newdata(model$models[[i]], nd, model$model_frames[[i]])
  as.numeric(.suest_predict_values(model$models[[i]], nd, type, engine, model$model_types[i]))
}

.suest_twopart_amount <- function(model, u, nd) {
  m <- .suest_twopart_values(model, u$second, nd)
  if (is.null(u$retransform)) m else exp(m) * u$factor
}

.suest_twopart_predict <- function(model, newdata, type = "response") {
  if (is.null(type)) type <- "response"
  if (!type %in% .suest_twopart_types)
    stop(paste0("For two-part models, 'type' must be \"response\" (the overall mean), ",
      "\"positive\" (the probability of a positive outcome), or \"amount\" (the mean ",
      "among positive outcomes)."), call. = FALSE)
  units <- model$twopart$units
  keep <- type == "response" | vapply(units, function(u) u$kind == "pair", logical(1))
  if (!all(keep))
    .suest_twopart_note(sprintf("type = \"%s\" applies to two-part models only; left out: %s.",
      type, paste(vapply(units[!keep], `[[`, "", "label"), collapse = ", ")))
  routed <- ".suest_model" %in% names(newdata)
  labels <- vapply(units, `[[`, "", "label")
  rowid <- if ("rowid" %in% names(newdata)) newdata$rowid else
    if (".suest_rowid" %in% names(newdata)) as.integer(newdata$.suest_rowid) else
      seq_len(nrow(newdata))
  # With one block of rows per model (.suest_model), every row needs a
  # prediction: rows of a model left out get NA under a kept label, so averages
  # drop them.
  pieces <- lapply(seq_along(units), function(k) {
    u <- units[[k]]
    rows <- if (routed) which(as.character(newdata$.suest_model) == u$label) else
      seq_len(nrow(newdata))
    if (!keep[k]) {
      if (!routed) return(NULL)
      return(data.frame(rowid = rowid[rows], group = labels[keep][1L], estimate = NA_real_,
        .order = rows))
    }
    if (!length(rows))
      stop(sprintf(paste0("No newdata rows were supplied for '%s'. The .suest_model column ",
        "says which model each row belongs to; with datagrid(), add .suest_model = c(%s)."),
        u$label, paste0('"', labels, '"', collapse = ", ")), call. = FALSE)
    nd <- newdata[rows, , drop = FALSE]
    estimate <- if (u$kind == "model") {
      .suest_twopart_values(model, u$model, nd)
    } else {
      p <- if (type != "amount") .suest_twopart_values(model, u$first, nd)
      m <- if (type != "positive") .suest_twopart_amount(model, u, nd)
      switch(type, response = p * m, positive = p, amount = m)
    }
    data.frame(rowid = rowid[rows], group = u$label, estimate = estimate, .order = rows)
  })
  out <- do.call(rbind, pieces)
  if (routed) out <- out[order(out$.order), , drop = FALSE]
  out$.order <- NULL
  rownames(out) <- NULL
  if (length(units) == 1L) {
    out$group <- NULL
  } else {
    out$group <- factor(out$group, levels = labels[keep])
  }
  out
}

# Show a note once per marginaleffects call, not on every prediction it makes.
.suest_twopart_note <- function(msg) {
  ns <- asNamespace("marginaleffects")
  for (frame in sys.frames()) {
    if (identical(topenv(frame), ns)) {
      if (exists(".suest_twopart_noted", envir = frame, inherits = FALSE)) return(invisible())
      assign(".suest_twopart_noted", TRUE, envir = frame)
      break
    }
  }
  message(msg)
}
