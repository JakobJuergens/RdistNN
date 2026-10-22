# Weight-geometry panels (plan panels NP-2, NP-3): weight decay by rank,
# cumulative concentration, and the weight-vs-distance "effective kernel"
# view. The data layer is plot_data_weights().

.dnn_profile_support_edge <- function(object, n) {
  meta <- object$meta
  method <- if (!is.null(meta$method)) {
    as.character(meta$method[[1L]])
  } else {
    NA_character_
  }
  if (!identical(method, "exact")) {
    return(NA_integer_)
  }

  if (inherits(object, "tdnn_np_result")) {
    if (is.null(meta$s1)) {
      return(NA_integer_)
    }
    edge <- as.integer(n) - as.integer(meta$s1[[1L]]) + 1L
    return(max(1L, min(as.integer(n), edge)))
  }

  if (is.null(meta$s) || is.null(meta$k)) {
    return(NA_integer_)
  }
  edge <- as.integer(n) - as.integer(meta$s[[1L]]) + as.integer(meta$k[[1L]])
  max(1L, min(as.integer(n), edge))
}

#' Plot DNN weight profiles, concentration, and effective kernels
#'
#' Three views of the localization weights stored on a fitted object.
#' `type = "profile"` draws weight against rank with a guide at the exact
#' nominal unaveraged rank-table support edge `n - s + k` for DNN fits (`n - s + 1` in
#' the paper-default `k = 1` case); `type = "cumulative"` draws the
#' cumulative weight share with guides at 0.5 and 0.9 and markers at the
#' top-1/5/10 shares; `type = "kernel"` draws weight against stored distance
#' — the effective data-driven kernel — with a guide at the rank-`s` radius.
#' The rank-`s` radius is a scale diagnostic, not the DNN support radius.
#'
#' When more than eight targets are stored, panels default to a subset
#' (first, median effective sample size, last); override with `targets =`.
#' Two-scale (TDNN) fits have signed weights: a dashed zero line is added
#' automatically whenever negative weights are shown.
#'
#' @param object A fitted result object holding the weight matrix `W`
#'   (stored by default; `type = "kernel"` additionally requires stored
#'   distances via `return = c("W", "dx")`).
#' @param type One of `"profile"`, `"cumulative"`, or `"kernel"`.
#' @param targets Optional integer vector of target ids to display, including
#'   target-specific profiles when comparing candidates.
#' @param log_weights Logical (default `FALSE`); for rank-domain views,
#'   drop non-positive weights and use a log10 y scale (the number of
#'   dropped weights is reported in the caption).
#' @param candidates Candidate rows to compare for grid fits. The default
#'   `"default"` preserves the stored weight-matrix display. Use `"all"` or
#'   a character vector of candidate ids with `type = "profile"` to compare
#'   weight profiles without refitting. Duplicate-covariate block averaging
#'   can make these profiles vary by target; the selected target profiles are
#'   shown in separate panels. Untied default profiles coincide across targets
#'   and are shown once unless `targets` is explicitly supplied.
#' @param ... Ignored; reserved for future options.
#'
#' @return A ggplot object.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 80
#'   x_mat <- matrix(rnorm(n * 2), n, 2)
#'   y_vec <- x_mat[, 1] + rnorm(n, sd = 0.3)
#'   fit <- np_dnn(
#'     X = x_mat, x0 = c(0, 0), Y = y_vec,
#'     k = 1L, s = 12L, se = FALSE, return = c("W", "dx")
#'   )
#'   plot_weights(fit, type = "kernel")
#' }
#' @export
plot_weights <- function(object,
                         type = c("profile", "cumulative", "kernel"),
                         targets = NULL,
                         log_weights = FALSE,
                         candidates = "default",
                         ...) {
  .dnn_reuse_scope(object)
  .plot_require_ggplot2("plot_weights")
  type <- match.arg(type)
  if (!.plot_candidates_is_default(candidates)) {
    return(.plot_weights_candidates(
      object = object,
      type = type,
      targets = targets,
      log_weights = log_weights,
      candidates = candidates
    ))
  }

  dw <- plot_data_weights(object)
  n <- max(dw$rank)
  all_targets <- sort(unique(dw$target_id))
  m <- length(all_targets)

  if (type == "kernel" && all(is.na(dw$distance))) {
    .abort_bad_arg(
      "object",
      "must store distances for type = \"kernel\"; refit with `return = c(\"W\", \"dx\")`",
      quoted = FALSE
    )
  }

  neff_by_target <- vapply(
    all_targets,
    function(j) 1 / sum(dw$weight[dw$target_id == j]^2),
    numeric(1L)
  )
  targets <- .plot_select_targets(all_targets, neff_by_target, targets)
  dw <- dw[dw$target_id %in% targets, , drop = FALSE]
  dw$target <- factor(dw$target_id)
  multi <- length(targets) > 1L

  captions <- character(0)
  signed <- any(dw$weight < 0)
  if (signed) {
    captions <- c(captions, "two-scale weights are signed")
  }

  base_aes <- function(xcol, ycol) {
    if (multi) {
      ggplot2::aes(
        x = .data[[xcol]], y = .data[[ycol]], colour = .data$target
      )
    } else {
      ggplot2::aes(x = .data[[xcol]], y = .data[[ycol]])
    }
  }
  col_args <- if (multi) list() else list(colour = .plot_col_single)

  s_meta <- object$meta$s
  support_edge <- .dnn_profile_support_edge(object, n = n)

  if (type == "profile") {
    if (!is.na(support_edge) && !is.null(object$.model$x_mat) &&
        anyDuplicated(.dnn_duplicate_groups(object$.model$x_mat))) {
      captions <- c(captions, "dotted guide: nominal rank-table edge before duplicate-block averaging")
    }
    if (isTRUE(log_weights)) {
      dropped <- sum(dw$weight <= 0)
      dw <- dw[dw$weight > 0, , drop = FALSE]
      captions <- c(
        captions,
        sprintf("log scale: %d non-positive weights dropped", dropped)
      )
    }
    p <- ggplot2::ggplot(dw, base_aes("rank", "weight")) +
      do.call(ggplot2::geom_step, col_args)
    if (signed && !isTRUE(log_weights)) {
      p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed")
    }
    if (!is.na(support_edge)) {
      p <- p + ggplot2::geom_vline(
        xintercept = support_edge,
        linetype = "dotted"
      )
    }
    if (isTRUE(log_weights)) {
      p <- p + ggplot2::scale_y_log10()
    }
    p <- p + ggplot2::labs(x = "Rank", y = "DNN weight")
  } else if (type == "cumulative") {
    p <- ggplot2::ggplot(dw, base_aes("rank", "cum_share")) +
      do.call(ggplot2::geom_line, col_args) +
      ggplot2::geom_hline(yintercept = c(0.5, 0.9), linetype = "dotted") +
      ggplot2::geom_point(
        data = dw[dw$rank %in% c(1L, 5L, 10L), , drop = FALSE]
      ) +
      ggplot2::labs(x = "Rank", y = "Cumulative weight share")
  } else {
    p <- ggplot2::ggplot(dw, base_aes("distance", "weight")) +
      do.call(ggplot2::geom_step, col_args)
    if (signed) {
      p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed")
    }
    if (!is.null(s_meta)) {
      radius <- vapply(
        targets,
        function(j) dw$distance[dw$target_id == j][as.integer(s_meta)],
        numeric(1L)
      )
      p <- p + ggplot2::geom_vline(
        xintercept = radius, linetype = "dotted"
      )
    }
    p <- p + ggplot2::labs(x = "Distance", y = "DNN weight")
  }

  if (multi) {
    p <- p + ggplot2::labs(colour = "target")
  }
  honesty <- .plot_honesty_caption(object)
  if (nzchar(honesty)) captions <- c(captions, honesty)
  if (length(captions) > 0L) {
    p <- p + ggplot2::labs(caption = paste(captions, collapse = "\n"))
  }
  p + .rdistnn_theme()
}

.plot_weights_candidates <- function(object, type, log_weights, candidates, targets = NULL) {
  if (!identical(type, "profile")) {
    .abort_bad_arg(
      "type",
      paste0(
        "candidate comparisons are available only for type = \"profile\"; ",
        "use candidate_fit() for target-specific kernel or cumulative views"
      ),
      quoted = FALSE
    )
  }
  payload <- .plot_candidate_weight_frame(
    object = object,
    candidates = candidates,
    fn = "plot_weights",
    targets = targets
  )
  dw <- payload$data
  captions <- character(0)
  signed <- any(dw$weight < 0)
  if (signed) {
    captions <- c(captions, "two-scale weights are signed")
  }
  if (isTRUE(log_weights)) {
    dropped <- sum(dw$weight <= 0)
    dw <- dw[dw$weight > 0, , drop = FALSE]
    captions <- c(
      captions,
      sprintf("log scale: %d non-positive weights dropped", dropped)
    )
  }
  if (!is.null(payload$drop_caption)) {
    captions <- c(captions, payload$drop_caption)
  }

  if ("target_id" %in% names(dw)) {
    dw$profile_group <- interaction(dw$candidate_id, dw$target_id)
  } else {
    dw$profile_group <- dw$candidate_id
  }
  p <- ggplot2::ggplot(
    dw,
    ggplot2::aes(
      x = .data$rank,
      y = .data$weight,
      colour = .data$candidate_label,
      linetype = .data$candidate_label,
      linewidth = .data$default_linewidth,
      group = .data$profile_group
    )
  ) +
    ggplot2::geom_step()
  if ("target_id" %in% names(dw)) p <- p + ggplot2::facet_wrap(~target_id)
  if (signed && !isTRUE(log_weights)) {
    p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed")
  }
  if (isTRUE(log_weights)) {
    p <- p + ggplot2::scale_y_log10()
  }
  p <- p +
    ggplot2::scale_linewidth_manual(
      values = c(candidate = 0.5, default = 0.9),
      guide = "none"
    ) +
    ggplot2::labs(
      x = "Rank",
      y = "DNN weight",
      colour = "candidate",
      linetype = "candidate"
    )
  if (length(captions) > 0L) {
    p <- p + ggplot2::labs(caption = paste(captions, collapse = "\n"))
  }
  p + .rdistnn_theme()
}
