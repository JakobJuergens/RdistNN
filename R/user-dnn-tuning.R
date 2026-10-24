# Resolve tuning orders and candidate identities before grid computation.

.dnn_candidate_format_lambda <- function(lambda) {
  label <- vapply(
    lambda,
    function(x) format(x, trim = TRUE, scientific = FALSE),
    character(1L)
  )
  label <- gsub("\\.", "p", label)
  gsub("-", "m", label)
}

.dnn_candidate_unique_ids <- function(base_ids) {
  seen <- new.env(parent = emptyenv())
  ids <- character(length(base_ids))
  for (i in seq_along(base_ids)) {
    base <- base_ids[[i]]
    count <- if (exists(base, envir = seen, inherits = FALSE)) {
      get(base, envir = seen, inherits = FALSE) + 1L
    } else {
      1L
    }
    assign(base, count, envir = seen)
    ids[[i]] <- if (count == 1L) base else sprintf("%s__%d", base, count)
  }
  ids
}

.dnn_candidate_default_index <- function(candidates, default) {
  idx_s <- which(as.integer(candidates$s) == as.integer(default$s))
  if (length(idx_s) == 0L) {
    return(integer(0L))
  }

  default_lambda <- default$meta[["lambda", exact = TRUE]]
  default_lambda_num <- suppressWarnings(as.numeric(default_lambda))
  if (!is.null(default_lambda) && length(default_lambda_num) == 1L &&
      !is.na(default_lambda_num) && is.finite(default_lambda_num)) {
    lambda <- suppressWarnings(as.numeric(as.character(candidates$lambda)))
    tol <- sqrt(.Machine$double.eps) * pmax(1, abs(default_lambda_num))
    idx_lambda <- idx_s[!is.na(lambda[idx_s]) & abs(lambda[idx_s] - default_lambda_num) <= tol]
    if (length(idx_lambda) > 0L) {
      return(idx_lambda[1L])
    }
  }

  idx_s[1L]
}

.dnn_prepend_synthetic_default_candidate <- function(candidates, default) {
  synthetic <- candidates[NA_integer_, , drop = FALSE]
  rownames(synthetic) <- NULL
  synthetic$s <- as.integer(default$s)
  synthetic$lambda <- NA_real_
  synthetic$feasible <- !isFALSE(default$feasibility$k_le_s) &&
    !isFALSE(default$feasibility$s_lt_n)
  synthetic$reason <- if (isFALSE(default$feasibility$k_le_s)) "k > s" else
    if (isFALSE(default$feasibility$s_lt_n)) "jackknife requires s < n" else "ok"
  rbind(synthetic, candidates)
}

.dnn_resolve_candidate_schema <- function(grid, default) {
  if (!is.list(grid) || !is.data.frame(grid$candidates)) {
    .dnn_abort("internal_invariant",
      "`tuning$grid$candidates` must be a data.frame.",
      hint = "Report a reproducible example of this internal inconsistency."
    )
  }
  candidates <- as.data.frame(grid$candidates, stringsAsFactors = FALSE)
  if (nrow(candidates) < 1L) {
    .dnn_abort("internal_invariant",
      "`tuning$grid$candidates` must contain at least one row.",
      hint = "Report a reproducible example of this internal inconsistency."
    )
  }
  required <- c("s", "feasible", "reason")
  missing <- setdiff(required, names(candidates))
  if (length(missing) > 0L) {
    .dnn_abort("internal_invariant", sprintf(
      "`tuning$grid$candidates` is missing required column(s): %s.",
      paste(missing, collapse = ", ")
    ), hint = "Report a reproducible example of this internal inconsistency.")
  }

  candidates$s <- as.integer(candidates$s)
  candidates$feasible <- as.logical(candidates$feasible)
  candidates$reason <- as.character(candidates$reason)
  if (!("lambda" %in% names(candidates))) {
    candidates$lambda <- NA_real_
  }
  candidates$lambda <- suppressWarnings(as.numeric(as.character(candidates$lambda)))

  default_idx <- .dnn_candidate_default_index(candidates, default)
  if (length(default_idx) == 0L) {
    candidates <- .dnn_prepend_synthetic_default_candidate(candidates, default)
    default_idx <- 1L
  }

  candidates$is_default <- FALSE
  candidates$is_default[default_idx[1L]] <- TRUE
  candidates$status <- ifelse(candidates$feasible, "ok", "infeasible")
  candidates$computed_key <- sprintf("s_%d", as.integer(candidates$s))

  base_ids <- sprintf("s_%d", as.integer(candidates$s))
  has_lambda <- !is.na(candidates$lambda)
  base_ids[has_lambda] <- paste0(
    "lambda_",
    .dnn_candidate_format_lambda(candidates$lambda[has_lambda])
  )
  candidates$candidate_id <- .dnn_candidate_unique_ids(base_ids)

  schema_cols <- c(
    "candidate_id", "computed_key", "is_default", "feasible",
    "reason", "status", "lambda", "s"
  )
  candidates <- candidates[, c(schema_cols, setdiff(names(candidates), schema_cols)), drop = FALSE]
  rownames(candidates) <- NULL
  grid$candidates <- candidates
  grid
}

.dnn_resolve_s <- function(s, tuning, n, p, k, se, jackknife, ctx = NULL) {
  needs_jk <- isTRUE(se) && !identical(jackknife, "none")

  tuning_source <- "tuning"
  if (is.null(tuning) && is.null(s)) {
    # Neither `s` nor `tuning` supplied: default to the paper undersmoothing
    # rule s = floor(0.4 * n^((kx + 1)/(kx + 2))) via tuning_dnn_paper(). A fixed
    # constant default (for example s = 1) yields uniform weights and hence
    # no localization at all, so it must never be the silent default.
    tuning <- tuning_dnn_paper()
    tuning_source <- "paper_default"
  }

  if (is.null(tuning)) {
    s_out <- .validate_scalar_int(s, "s")
    k_ok <- as.integer(k) <= as.integer(s_out)
    s_ok <- if (needs_jk) as.integer(s_out) < as.integer(n) else TRUE
    if (!k_ok) {
      .abort_bad_arg("k", sprintf("must satisfy k <= s (got k=%d, s=%d)", as.integer(k), as.integer(s_out)))
    }
    if (!s_ok) {
      .abort_bad_arg(
        "s",
        "must satisfy s < n when delete-1 jackknife is requested",
        expected = sprintf("s < n=%d", as.integer(n)),
        actual = sprintf("s=%d (with `se`=%s, `jackknife`='%s')", as.integer(s_out), as.character(se), jackknife)
      )
    }
    default <- list(
      s = as.integer(s_out),
      meta = list(name = "manual", n = as.integer(n), kx = as.integer(p)),
      feasibility = list(k_le_s = k_ok, s_lt_n = s_ok)
    )
    grid <- .dnn_resolve_candidate_schema(
      list(
        candidates = data.frame(
          lambda = NA_real_,
          s = as.integer(s_out),
          feasible = TRUE,
          reason = "ok",
          stringsAsFactors = FALSE
        ),
        meta = list(name = "manual", n = as.integer(n), kx = as.integer(p))
      ),
      default = default
    )
    payload <- list(
      source = "manual",
      default = default,
      grid = grid
    )
    return(list(s = as.integer(s_out), tuning = payload, tuning_used = FALSE))
  }

  if (!is.null(s)) {
    .abort_bad_arg(
      "s/tuning",
      "cannot provide both explicit `s` and `tuning`; omit `s` when `tuning` is used",
      quoted = FALSE
    )
  }

  default <- tuning_eval(
    tuning = tuning,
    n = n,
    kx = p,
    k = k,
    se = se,
    jackknife = jackknife,
    ctx = ctx,
    mode = "default"
  )
  grid <- tuning_eval(
    tuning = tuning,
    n = n,
    kx = p,
    k = k,
    se = se,
    jackknife = jackknife,
    ctx = ctx,
    mode = "grid"
  )

  if (!is.list(default$feasibility)) {
    .dnn_abort("internal_invariant",
      "`tuning_eval(..., mode='default')` must return `feasibility`.",
      hint = "Report a reproducible example of this internal inconsistency."
    )
  }

  s_out <- .validate_scalar_int(default$s, "s")
  if (!isTRUE(default$feasibility$k_le_s)) {
    .abort_bad_arg("k", sprintf("must satisfy k <= s (got k=%d, s=%d)", as.integer(k), as.integer(s_out)))
  }
  if (needs_jk && !isTRUE(default$feasibility$s_lt_n)) {
    .abort_bad_arg(
      "s",
      "must satisfy s < n when delete-1 jackknife is requested",
      expected = sprintf("s < n=%d", as.integer(n)),
      actual = sprintf("s=%d (with `se`=%s, `jackknife`='%s')", as.integer(s_out), as.character(se), jackknife)
    )
  }

  payload <- list(
    source = tuning_source,
    default = default,
    grid = .dnn_resolve_candidate_schema(grid, default = default)
  )
  list(s = as.integer(s_out), tuning = payload, tuning_used = TRUE)
}

.dnn_validate_retained_tuning <- function(payload, fitted_s, field = "object$meta$tuning") {
  invalid <- function(field, reason) {
    .dnn_abort("state_invalid",
      sprintf("Retained `%s` %s. Refit with the current estimator.", field, reason), argument = field)
  }
  record <- function(value, fields, path) {
    if (!is.list(value)) invalid(path, "must retain the current tuning record")
    for (field in fields) {
      if (sum(names(value) == field, na.rm = TRUE) != 1L) {
        invalid(paste0(path, "$", field), "must be retained exactly once")
      }
    }
    value
  }
  integer <- function(value, path) {
    if (!is.null(dim(value))) invalid(path, "must retain an integer scalar")
    .dnn_validate_retained(.validate_scalar_int(value, path))
  }
  payload <- record(payload, c("source", "default", "grid"), field)
  source <- payload[["source", exact = TRUE]]
  if (!is.character(source) || !is.null(dim(source)) || length(source) != 1L ||
      is.na(source) || !source %in% c("manual", "tuning", "paper_default")) {
    invalid(paste0(field, "$source"), "must retain the tuning source")
  }
  default <- record(payload[["default", exact = TRUE]], c("s", "meta", "feasibility"),
    paste0(field, "$default"))
  if (integer(default[["s", exact = TRUE]], paste0(field, "$default$s")) !=
      fitted_s) {
    invalid(paste0(field, "$default$s"), "must match the fitted order")
  }
  grid <- record(payload[["grid", exact = TRUE]], c("candidates", "meta"), paste0(field, "$grid"))
  for (location in c("default", "grid")) {
    path <- paste0(field, "$", location, "$meta")
    descriptor <- record(payload[[location]][["meta", exact = TRUE]], c("name", "n", "kx"), path)
    name <- descriptor[["name", exact = TRUE]]
    if (!is.character(name) || !is.null(dim(name)) || length(name) != 1L || is.na(name) || !nzchar(name)) {
      invalid(paste0(path, "$name"), "must retain a non-empty name")
    }
    for (key in c("n", "kx")) integer(descriptor[[key, exact = TRUE]], paste0(path, "$", key))
  }
  feasibility_path <- paste0(field, "$default$feasibility")
  feasible <- record(default[["feasibility", exact = TRUE]], c("k_le_s", "s_lt_n"), feasibility_path)
  for (flag in c("k_le_s", "s_lt_n")) {
    if (!identical(feasible[[flag, exact = TRUE]], TRUE)) {
      invalid(paste0(feasibility_path, "$", flag), "must retain a feasible default")
    }
  }
  candidates <- grid[["candidates", exact = TRUE]]
  fields <- c("candidate_id", "computed_key", "is_default", "feasible", "reason", "status", "lambda", "s")
  if (!is.data.frame(candidates) || !all(fields %in% names(candidates)) || nrow(candidates) < 1L) {
    invalid(paste0(field, "$grid$candidates"), "must retain the current candidate schema")
  }
  payload
}
