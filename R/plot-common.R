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
