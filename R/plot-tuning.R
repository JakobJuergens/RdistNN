# Kernel-order sensitivity panel (plan panel NP-6): estimates with
# confidence bounds across the candidate kernel orders of a tuning grid,
# with the fitted object's own s marked. Answers "how sensitive is this
# estimate to s?" by refitting at each candidate.

.plot_tuning_candidates <- function(object, tuning) {
  meta <- object$meta
  grid <- tuning_eval(
    tuning = tuning, n = meta$n, kx = meta$p, k = meta$k,
    se = TRUE, jackknife = "fast", mode = "grid"
  )
  s_cand <- sort(unique(c(grid$candidates$s, as.integer(meta$s))))
  s_cand[s_cand >= max(2L, as.integer(meta$k)) & s_cand < meta$n]
}

.plot_tuning_panel <- function(object, refit_fn, tuning, targets, level) {
  s_cand <- .plot_tuning_candidates(object, tuning)
  fits <- lapply(s_cand, refit_fn)

  frames <- Map(function(fit, s) {
    d <- plot_data_estimates(fit, level = level)
    d$s <- s
    d
  }, fits, s_cand)
  df <- do.call(rbind, frames)

  all_targets <- sort(unique(df$target_id))
  if (is.null(targets) && length(all_targets) > 3L) {
    targets <- unique(c(
      all_targets[1L],
      all_targets[ceiling(length(all_targets) / 2)],
      all_targets[length(all_targets)]
    ))
  } else if (is.null(targets)) {
    targets <- all_targets
  } else {
    targets <- unique(as.integer(targets))
    if (!all(targets %in% all_targets)) {
      .abort_bad_arg(
        "targets",
        sprintf("must be target ids in 1..%d", max(all_targets)),
        quoted = FALSE
      )
    }
  }
  df <- df[df$target_id %in% targets, , drop = FALSE]
  df$target <- factor(sprintf("target %d", df$target_id),
                      levels = sprintf("target %d", targets))
  s_selected <- as.integer(object$meta$s)
  df$selected <- df$s == s_selected

  has_se <- !all(is.na(df$se))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$s, y = .data$estimate)) +
    ggplot2::geom_vline(xintercept = s_selected, linetype = "dotted")
  if (has_se) {
    p <- p + ggplot2::geom_pointrange(
      ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper,
                   shape = .data$selected),
      colour = .plot_col_single
    )
  } else {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(shape = .data$selected),
      colour = .plot_col_single
    )
  }
  p <- p +
    ggplot2::scale_shape_manual(
      values = c(`TRUE` = 16, `FALSE` = 1), guide = "none"
    ) +
    ggplot2::scale_x_continuous(breaks = sort(unique(df$s)))
  if (length(targets) > 1L) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$target),
                                 scales = "free_y")
  }
  src <- if (!is.null(object$meta$tuning$source)) {
    object$meta$tuning$source
  } else {
    "manual"
  }
  p + ggplot2::labs(
    x = "kernel order s", y = "Estimate",
    caption = sprintf(
      "filled point / dotted line: fitted s = %d (%s); n = %d",
      s_selected, src, as.integer(object$meta$n)
    )
  ) + .rdistnn_theme()
}

#' Plot estimate sensitivity to the kernel order
#'
#' Refits the estimator at each candidate kernel order of a tuning grid
#' (default: the paper rule's kappa grid via [tuning_dnn_paper()], always
#' including the fitted `s`) and plots the estimates with confidence
#' bounds against `s`, marking the fitted kernel order with a filled point
#' and a dotted guide. Each candidate refit recomputes the delete-one
#' jackknife, so the panel costs roughly one fit per candidate.
#'
#' Refits use the training data stored on the fit; fits from package
#' versions predating that storage must be re-estimated first.
#'
#' @param object A fitted `np_result` or `cate_result`.
#' @param tuning A tuning specification for the candidate grid (default
#'   [tuning_dnn_paper()]).
#' @param targets Optional integer vector of target ids (defaults to at
#'   most three: first, middle, last).
#' @param level Confidence level for the bounds (default `0.95`).
#' @param se Logical; compute jackknife bounds for each candidate
#'   (default `TRUE`).
#' @param ... Ignored; reserved for future options.
#'
#' @return A ggplot object.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 100
#'   x_mat <- matrix(rnorm(n * 2), n, 2)
#'   y_vec <- x_mat[, 1] + rnorm(n, sd = 0.3)
#'   fit <- np_dnn(X = x_mat, x0 = c(0, 0), Y = y_vec, k = 1L)
#'   plot_tuning(fit)
#' }
#' }
#' @export
plot_tuning <- function(object, tuning = tuning_dnn_paper(),
                        targets = NULL, level = 0.95, se = TRUE, ...) {
  UseMethod("plot_tuning")
}

#' @rdname plot_tuning
#' @export
plot_tuning.default <- function(object, tuning = tuning_dnn_paper(),
                                targets = NULL, level = 0.95, se = TRUE,
                                ...) {
  .abort_bad_arg(
    "object",
    sprintf(
      "plot_tuning() is not implemented for class '%s'",
      class(object)[1L]
    ),
    quoted = FALSE
  )
}

#' @rdname plot_tuning
#' @export
plot_tuning.np_result <- function(object, tuning = tuning_dnn_paper(),
                                  targets = NULL, level = 0.95, se = TRUE,
                                  ...) {
  .plot_require_ggplot2("plot_tuning")
  model <- object$.model
  if (!is.list(model) || is.null(model$x_mat) || is.null(model$y_vec)) {
    .abort_bad_arg(
      "object",
      "does not store its training data (fit predates this feature); re-estimate before calling plot_tuning()",
      quoted = FALSE
    )
  }
  refit_fn <- function(s) {
    np_dnn(
      X = model$x_mat, X0 = object$x0, Y = model$y_vec,
      distance = model$distance, method = model$method,
      k = object$meta$k, s = s, normalize = model$normalize,
      ctx = model$ctx, se = se, return = character(0)
    )
  }
  .plot_tuning_panel(object, refit_fn, tuning, targets, level)
}

#' @rdname plot_tuning
#' @export
plot_tuning.cate_result <- function(object, tuning = tuning_dnn_paper(),
                                    targets = NULL, level = 0.95,
                                    se = TRUE, ...) {
  .plot_require_ggplot2("plot_tuning")
  model <- object$.model
  if (!is.list(model) || is.null(model$x_mat) || is.null(model$psi)) {
    .abort_bad_arg(
      "object",
      "does not store its training data (fit predates this feature); re-estimate before calling plot_tuning()",
      quoted = FALSE
    )
  }
  refit_fn <- function(s) {
    cate_dnn(
      X = model$x_mat, X0 = object$x0, pseudo_outcome = model$psi,
      distance = model$distance, method = model$method,
      k = object$meta$k, s = s, normalize = model$normalize,
      ctx = model$ctx, se = se, return = character(0)
    )
  }
  .plot_tuning_panel(object, refit_fn, tuning, targets, level)
}
