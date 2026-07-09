# Shared internal helpers for the ggplot2-based plotting layer (Phase 1).
#
# ggplot2 is a suggested dependency: every public builder guards through
# .plot_require_ggplot2() and fails with an actionable message when it is
# not installed. All helpers here are internal.
# See markdown_notes/plotting_layer_implementation_plan.md.

# The `.data` pronoun inside ggplot2::aes() is supplied by the tidy-eval
# data mask at run time; this declaration only silences the R CMD check
# visible-binding NOTE without importing rlang.
utils::globalVariables(".data")

.plot_pkg_installed <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

.plot_require_ggplot2 <- function(fn) {
  if (!.plot_pkg_installed("ggplot2")) {
    .abort_bad_state(sprintf(
      "%s() requires the suggested package 'ggplot2'; install it with install.packages(\"ggplot2\").",
      fn
    ))
  }
  invisible(NULL)
}

# Okabe-Ito assignments per the plotting plan (3.0): colour-blind safe,
# hard-coded so no palette dependency is needed.
.plot_col_treated <- "#0072B2"
.plot_col_control <- "#D55E00"
.plot_col_single <- "#1B1B1B"

.rdistnn_theme <- function() {
  ggplot2::theme_minimal(base_size = 11)
}

# Profile detection (plan 2.5): if the targets vary in exactly one x0_*
# column, that column is the x axis; `x_var` overrides by coordinate name
# (from attr(df, "x0_names")), by x0_* column name, or by column index.
# Returns list(var = <column name or NULL>, label = <axis label>).
.plot_profile_axis <- function(df, x_var = NULL) {
  x0_cols <- grep("^x0_[0-9]+$", names(df), value = TRUE)
  x0_names <- attr(df, "x0_names")

  label_for <- function(col) {
    idx <- as.integer(sub("^x0_", "", col))
    if (!is.null(x0_names) && length(x0_names) >= idx && nzchar(x0_names[idx])) {
      x0_names[idx]
    } else {
      sprintf("x[%d]", idx)
    }
  }

  if (!is.null(x_var)) {
    if (length(x0_cols) == 0L) {
      .abort_bad_arg("x_var", "was supplied, but the fit stores no target coordinates")
    }
    if (is.numeric(x_var)) {
      col <- sprintf("x0_%d", as.integer(x_var))
    } else if (!is.null(x0_names) && x_var %in% x0_names) {
      col <- sprintf("x0_%d", match(x_var, x0_names))
    } else {
      col <- as.character(x_var)
    }
    if (!col %in% x0_cols) {
      .abort_bad_arg(
        "x_var",
        sprintf("does not identify a stored target coordinate (got '%s')", as.character(x_var)),
        quoted = FALSE
      )
    }
    return(list(var = col, label = label_for(col)))
  }

  if (length(x0_cols) == 0L) {
    return(list(var = NULL, label = "target"))
  }
  varying <- x0_cols[vapply(
    x0_cols,
    function(cl) length(unique(df[[cl]])) > 1L,
    logical(1L)
  )]
  if (length(varying) == 1L) {
    return(list(var = varying, label = label_for(varying)))
  }
  list(var = NULL, label = "target")
}

# Honesty caption (plan 3.0): median n_eff with range, plus s and its
# provenance. Degrades gracefully when diagnostics are unavailable (for
# example when W was not stored).
.plot_honesty_caption <- function(object) {
  meta <- object$meta
  s_txt <- if (!is.null(meta$s)) {
    src <- if (!is.null(meta$tuning$source)) meta$tuning$source else "manual"
    sprintf("s = %d (%s)", as.integer(meta$s), src)
  } else if (!is.null(meta$s1)) {
    sprintf("s1 = %d, s2 = %d", as.integer(meta$s1), as.integer(meta$s2))
  } else {
    NULL
  }
  neff_txt <- tryCatch(
    {
      dd <- plot_data_diagnostics(object)
      sprintf(
        "median n_eff = %.1f (%.1f-%.1f)",
        stats::median(dd$n_eff), min(dd$n_eff), max(dd$n_eff)
      )
    },
    error = function(e) NULL
  )
  paste(c(neff_txt, s_txt), collapse = "; ")
}

# Estimate display mode (plan 3.0 with the C-2 forest rule): profile when a
# profile axis exists and m >= 4; forest for small non-profile multi-target
# sets; pointrange otherwise.
.plot_estimates_mode <- function(style, has_profile, m) {
  if (style != "auto") {
    return(style)
  }
  if (has_profile && m >= 4L) {
    return("profile")
  }
  if (!has_profile && m > 1L && m <= 12L) {
    return("forest")
  }
  "pointrange"
}

# Default target subsetting for weight-domain and density panels: with more
# than eight targets, show the first, the median-n_eff, and the last.
.plot_select_targets <- function(all_targets, neff_by_target, targets = NULL) {
  m <- length(all_targets)
  if (!is.null(targets)) {
    targets <- unique(as.integer(targets))
    if (!all(targets %in% all_targets)) {
      .abort_bad_arg(
        "targets",
        sprintf("must be target ids in 1..%d", max(all_targets)),
        quoted = FALSE
      )
    }
    return(targets)
  }
  if (m <= 8L) {
    return(all_targets)
  }
  med <- all_targets[which.min(abs(neff_by_target - stats::median(neff_by_target)))]
  unique(c(all_targets[1L], med, all_targets[m]))
}

# Diagnostics frame augmented with stored target coordinates so overlap and
# diagnostics panels can reuse the profile-axis detection.
.plot_diag_frame <- function(object) {
  dd <- plot_data_diagnostics(object)
  x0 <- object$x0
  if (!is.null(x0) && nrow(as.matrix(x0)) == nrow(dd)) {
    x0 <- as.matrix(x0)
    for (j in seq_len(ncol(x0))) {
      dd[[sprintf("x0_%d", j)]] <- as.numeric(x0[, j])
    }
    attr(dd, "x0_names") <- colnames(x0)
  }
  dd
}

# Dashboard annotation string (plan C-6): metadata contract consumed by the
# dashboard tests. Fields absent from meta are omitted.
.plot_dashboard_annotation <- function(object) {
  meta <- object$meta
  s_txt <- if (!is.null(meta$s)) {
    src <- if (!is.null(meta$tuning$source)) meta$tuning$source else "manual"
    sprintf("s = %d (%s)", as.integer(meta$s), src)
  } else if (!is.null(meta$s1)) {
    sprintf("s1 = %d, s2 = %d", as.integer(meta$s1), as.integer(meta$s2))
  } else {
    NULL
  }
  parts <- c(
    if (!is.null(meta$n)) sprintf("n = %d", as.integer(meta$n)),
    if (!is.null(meta$m)) sprintf("m = %d", as.integer(meta$m)),
    if (!is.null(meta$k)) sprintf("k = %d", as.integer(meta$k)),
    s_txt,
    if (!is.null(meta$distance)) sprintf("distance = %s", meta$distance),
    if (!is.null(meta$jackknife)) sprintf("jackknife = %s", meta$jackknife)
  )
  paste(parts, collapse = ", ")
}

# Compact target labels for forest layouts, built from stored coordinates:
# "X1 = 0, X2 = 1" (or "x[j] = v" without names); plain "target j" when no
# coordinates are stored.
.plot_target_labels <- function(df) {
  x0_cols <- grep("^x0_[0-9]+$", names(df), value = TRUE)
  ids <- sort(unique(df$target_id))
  if (length(x0_cols) == 0L) {
    return(stats::setNames(sprintf("target %d", ids), ids))
  }
  x0_names <- attr(df, "x0_names")
  labs <- vapply(ids, function(j) {
    row <- df[match(j, df$target_id), x0_cols, drop = FALSE]
    parts <- vapply(seq_along(x0_cols), function(i) {
      nm <- if (!is.null(x0_names) && nzchar(x0_names[i])) {
        x0_names[i]
      } else {
        sprintf("x[%d]", i)
      }
      sprintf("%s = %s", nm, format(row[[i]], digits = 3))
    }, character(1L))
    paste(parts, collapse = ", ")
  }, character(1L))
  stats::setNames(labs, ids)
}

.plot_candidates_is_default <- function(candidates) {
  is.character(candidates) &&
    length(candidates) == 1L &&
    identical(candidates[[1L]], "default")
}

.plot_require_grid <- function(object, fn) {
  grid <- object$tuning_grid
  if (is.null(grid)) {
    .abort_bad_arg(
      "candidates",
      sprintf(
        "%s() can compare non-default candidates only for grid fits; refit with `tuning_mode = \"grid\"`.",
        fn
      ),
      quoted = FALSE
    )
  }
  if (!is.list(grid) || !is.data.frame(grid$candidates) ||
      !is.data.frame(grid$estimates) || !is.data.frame(grid$diagnostics)) {
    .abort_bad_state("`object$tuning_grid` is malformed.")
  }
  grid
}

.plot_format_candidate_value <- function(x) {
  format(as.numeric(x), trim = TRUE, scientific = FALSE)
}

.plot_candidate_label <- function(candidate_row) {
  s <- as.integer(candidate_row$s[[1L]])
  lambda <- if ("lambda" %in% names(candidate_row)) {
    suppressWarnings(as.numeric(candidate_row$lambda[[1L]]))
  } else {
    NA_real_
  }
  if (!is.na(lambda) && is.finite(lambda)) {
    sprintf("lambda = %s (s = %d)", .plot_format_candidate_value(lambda), s)
  } else {
    sprintf("s = %d", s)
  }
}

.plot_add_candidate_labels <- function(candidates) {
  candidates$candidate_label <- vapply(
    seq_len(nrow(candidates)),
    function(i) .plot_candidate_label(candidates[i, , drop = FALSE]),
    character(1L)
  )
  candidates
}

.plot_order_candidates <- function(candidates) {
  candidates <- .plot_add_candidate_labels(candidates)
  candidates[order(as.integer(candidates$s), candidates$candidate_id), , drop = FALSE]
}

.plot_resolve_candidates <- function(object, candidates, fn) {
  grid <- .plot_require_grid(object, fn = fn)
  rows <- grid$candidates
  if (.plot_candidates_is_default(candidates)) {
    rows <- rows[as.logical(rows$is_default), , drop = FALSE]
  } else if (is.character(candidates) &&
             length(candidates) == 1L &&
             identical(candidates[[1L]], "all")) {
    rows <- rows
  } else if (is.character(candidates) &&
             length(candidates) > 0L &&
             !anyNA(candidates) &&
             all(nzchar(candidates))) {
    idx <- match(candidates, rows$candidate_id)
    if (anyNA(idx)) {
      .abort_bad_arg(
        "candidates",
        sprintf(
          "Unknown candidate id(s): %s. Available candidate ids: %s",
          paste(candidates[is.na(idx)], collapse = ", "),
          paste(rows$candidate_id, collapse = ", ")
        ),
        quoted = FALSE
      )
    }
    rows <- rows[idx, , drop = FALSE]
  } else {
    .abort_bad_arg(
      "candidates",
      "must be \"default\", \"all\", or a character vector of candidate ids",
      quoted = FALSE
    )
  }
  .plot_order_candidates(rows)
}

.plot_candidate_drop_caption <- function(candidates) {
  dropped <- candidates[as.character(candidates$status) != "ok", , drop = FALSE]
  if (nrow(dropped) < 1L) {
    return(NULL)
  }
  details <- paste(
    sprintf(
      "%s (%s: %s)",
      dropped$candidate_id,
      dropped$status,
      dropped$reason
    ),
    collapse = ", "
  )
  sprintf(
    "Dropped %d unavailable candidate(s): %s.",
    nrow(dropped),
    details
  )
}

.plot_ok_candidates <- function(candidates) {
  candidates[as.character(candidates$status) == "ok", , drop = FALSE]
}

.plot_grid_estimate_coef <- function(object, component) {
  coef <- as.character(component)
  coef[is.na(component) | !nzchar(coef)] <- if (inherits(object, "np_result")) {
    "mu"
  } else if (inherits(object, "cate_result")) {
    "tau"
  } else {
    "estimate"
  }
  coef
}

.plot_attach_x0 <- function(out, object) {
  if (is.null(object$x0)) {
    return(out)
  }
  x0 <- as.matrix(object$x0)
  if (!"target_id" %in% names(out)) {
    return(out)
  }
  target_id <- as.integer(out$target_id)
  if (anyNA(target_id) || any(target_id < 1L | target_id > nrow(x0))) {
    return(out)
  }
  x0_df <- as.data.frame(x0[target_id, , drop = FALSE])
  names(x0_df) <- paste0("x0_", seq_len(ncol(x0)))
  out <- cbind(out, x0_df)
  attr(out, "x0_names") <- colnames(x0)
  out
}

.plot_candidate_levels <- function(candidates) {
  as.character(.plot_order_candidates(.plot_ok_candidates(candidates))$candidate_label)
}

.plot_order_candidate_frame <- function(df) {
  if (nrow(df) < 1L) {
    return(df)
  }
  df[order(
    as.logical(df$is_default),
    as.integer(df$target_id),
    as.character(df$coef),
    as.integer(df$s),
    as.character(df$candidate_id)
  ), , drop = FALSE]
}

.plot_candidate_estimates_frame <- function(object, candidates, level, fn) {
  grid <- .plot_require_grid(object, fn = fn)
  selected <- .plot_resolve_candidates(object, candidates = candidates, fn = fn)
  ok <- .plot_ok_candidates(selected)
  if (nrow(ok) < 1L) {
    .abort_bad_arg(
      "candidates",
      "did not select any available candidates to display",
      quoted = FALSE
    )
  }

  rows <- grid$estimates[grid$estimates$candidate_id %in% ok$candidate_id, , drop = FALSE]
  rows$.candidate_order <- match(rows$candidate_id, ok$candidate_id)
  rows <- rows[order(
    rows$.candidate_order,
    as.integer(rows$target_id),
    as.character(rows$component)
  ), , drop = FALSE]
  rows$.candidate_order <- NULL
  idx <- match(rows$candidate_id, ok$candidate_id)
  z <- stats::qnorm(1 - (1 - .validate_scalar_prob(level, arg = "level")) / 2)
  out <- data.frame(
    candidate_id = rows$candidate_id,
    candidate_label = ok$candidate_label[idx],
    computed_key = rows$computed_key,
    lambda = rows$lambda,
    s = as.integer(rows$s),
    is_default = as.logical(ok$is_default[idx]),
    target_id = as.integer(rows$target_id),
    coef = .plot_grid_estimate_coef(object, rows$component),
    estimate = as.numeric(rows$estimate),
    se = as.numeric(rows$se),
    status = rows$status,
    reason = rows$reason,
    stringsAsFactors = FALSE
  )
  out$ci_lower <- out$estimate - z * out$se
  out$ci_upper <- out$estimate + z * out$se
  out <- .plot_attach_x0(out, object)
  out$candidate_label <- factor(
    out$candidate_label,
    levels = .plot_candidate_levels(selected)
  )
  out$default_linewidth <- factor(
    ifelse(out$is_default, "default", "candidate"),
    levels = c("candidate", "default")
  )
  attr(out, "level") <- level
  list(
    data = .plot_order_candidate_frame(out),
    drop_caption = .plot_candidate_drop_caption(selected),
    selected = selected,
    ok = ok
  )
}

.plot_candidate_weight_frame <- function(object, candidates, fn) {
  selected <- .plot_resolve_candidates(object, candidates = candidates, fn = fn)
  ok <- .plot_ok_candidates(selected)
  if (nrow(ok) < 1L) {
    .abort_bad_arg(
      "candidates",
      "did not select any available candidates to display",
      quoted = FALSE
    )
  }
  ranks <- .candidate_cached_ranks(object)
  n <- nrow(ranks)
  frames <- vector("list", nrow(ok))
  for (i in seq_len(nrow(ok))) {
    row <- ok[i, , drop = FALSE]
    spec <- .candidate_rank_weight_spec(object, row)
    tab <- .rank_weight_table(spec, n)
    if (is.null(tab)) {
      ctx <- if (is.list(object$.model)) object$.model$ctx else NULL
      tab <- rank_weight_eval(spec, seq_len(n), ctx = ctx)
    }
    total <- sum(tab)
    frames[[i]] <- data.frame(
      candidate_id = row$candidate_id,
      candidate_label = row$candidate_label,
      s = as.integer(row$s),
      is_default = as.logical(row$is_default),
      rank = seq_len(n),
      weight = as.numeric(tab),
      cum_share = if (is.finite(total) && abs(total) > 0) cumsum(tab) / total else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, frames)
  out$candidate_label <- factor(
    out$candidate_label,
    levels = .plot_candidate_levels(selected)
  )
  out$default_linewidth <- factor(
    ifelse(out$is_default, "default", "candidate"),
    levels = c("candidate", "default")
  )
  out <- out[order(as.logical(out$is_default), out$rank, out$s), , drop = FALSE]
  list(
    data = out,
    drop_caption = .plot_candidate_drop_caption(selected),
    selected = selected,
    ok = ok
  )
}
