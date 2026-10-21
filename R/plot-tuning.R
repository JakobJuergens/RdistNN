# Kernel-order sensitivity panel (plan panel NP-6): estimates with
# confidence bounds across the candidate kernel orders of a tuning grid,
# with the fitted object's own s marked. Answers "how sensitive is this
# estimate to s?" from stored or internally evaluated candidate grids.

.plot_tuning_jackknife <- function(se) {
  if (isTRUE(se)) "fast" else "none"
}

.plot_tuning_payload <- function(object, tuning, se) {
  payload <- if (!is.null(object$meta$tuning)) object$meta$tuning else object$tuning
  if (is.null(tuning) && is.list(payload) &&
      !identical(payload$source, "manual") &&
      !is.null(payload$grid$candidates) &&
      is.data.frame(payload$grid$candidates)) {
    payload$grid <- .tuning_validate_grid_output(
      payload$grid$candidates, name = "retained fit", n = object$meta$n,
      kx = object$meta$p, k = object$meta$k, se = se,
      jackknife = .plot_tuning_jackknife(se))
    payload$grid <- .dnn_resolve_candidate_schema(payload$grid, default = payload$default)
    return(payload)
  }
  if (is.null(tuning)) tuning <- tuning_dnn_paper()
  meta <- object$meta
  model <- object$.model
  ctx <- if (is.list(model)) model$ctx else NULL
  grid <- tuning_eval(
    tuning = tuning,
    n = meta$n,
    kx = meta$p,
    k = meta$k,
    se = se,
    jackknife = .plot_tuning_jackknife(se),
    ctx = ctx,
    mode = "grid"
  )
  default <- list(
    s = as.integer(meta$s),
    meta = list(name = "fit", n = as.integer(meta$n), kx = as.integer(meta$p)),
    feasibility = list(
      k_le_s = as.integer(meta$k) <= as.integer(meta$s),
      s_lt_n = if (isTRUE(se)) as.integer(meta$s) < as.integer(meta$n) else TRUE
    )
  )
  list(
    source = "fit",
    default = default,
    grid = .dnn_resolve_candidate_schema(grid, default = default)
  )
}

.plot_tuning_ranks <- function(object) {
  ranking <- .dnn_retained_ranking(object)
  if (!is.null(ranking)) return(ranking$ranks)
  model <- object$.model
  if (!is.list(model) || is.null(model$x_mat)) {
    .abort_bad_arg(
      "object",
      "does not contain the required fitted training data; refit before calling plot_tuning()",
      quoted = FALSE
    )
  }
  .abort_bad_state("This fit does not retain the geometry needed to reconstruct ranks.")
}

.plot_tuning_grid_np <- function(object, tuning, se) {
  if (!is.null(object$tuning_grid)) {
    return(object$tuning_grid)
  }
  model <- object$.model
  if (!is.list(model) || is.null(model$x_mat) || is.null(model$y_vec)) {
    .abort_bad_arg(
      "object",
      "does not contain the required fitted training data; refit before calling plot_tuning()",
      quoted = FALSE
    )
  }
  payload <- .plot_tuning_payload(object, tuning = tuning, se = se)
  suppressWarnings(.np_dnn_tuning_grid(
    candidates = payload$grid$candidates,
    ranks = .plot_tuning_ranks(object),
    mu_hat = object$mu,
    se_hat = object$se,
    x_mat = model$x_mat,
    x0_mat = as.matrix(object$x0),
    y_vec = model$y_vec,
    distance = model$distance,
    method = model$method,
    k = object$meta$k,
    normalize = model$normalize,
    ctx = model$ctx,
    se = se,
    jackknife = .plot_tuning_jackknife(se),
    tuning_payload = payload,
    reuse_default = FALSE
  ))
}

.plot_tuning_grid_cate <- function(object, tuning, se) {
  if (!is.null(object$tuning_grid)) {
    return(object$tuning_grid)
  }
  model <- object$.model
  if (!is.list(model) || is.null(model$x_mat) || is.null(model$psi)) {
    .abort_bad_arg(
      "object",
      "does not contain the required fitted training data; refit before calling plot_tuning()",
      quoted = FALSE
    )
  }
  payload <- .plot_tuning_payload(object, tuning = tuning, se = se)
  suppressWarnings(.cate_dnn_tuning_grid(
    candidates = payload$grid$candidates,
    ranks = .plot_tuning_ranks(object),
    tau_hat = object$tau,
    se_hat = object$se,
    x_mat = model$x_mat,
    x0_mat = as.matrix(object$x0),
    psi = model$psi,
    distance = model$distance,
    method = model$method,
    k = object$meta$k,
    normalize = model$normalize,
    ctx = model$ctx,
    se = se,
    jackknife = .plot_tuning_jackknife(se),
    tuning_payload = payload,
    reuse_default = FALSE
  ))
}

.plot_tuning_plqr_score <- function(object) {
  if (!is.null(object$score)) {
    return(as.matrix(object$score))
  }
  model <- object$.model
  if (!is.list(model) || is.null(model$d_resid_mat) ||
      is.null(model$below_threshold)) {
    .abort_bad_arg(
      "object",
      "does not store PLQR score components; re-estimate before calling plot_tuning()",
      quoted = FALSE
    )
  }
  as.matrix(model$d_resid_mat) *
    as.numeric(object$tau - model$below_threshold)
}

.plot_tuning_grid_plqr <- function(object, tuning, se) {
  if (!is.null(object$tuning_grid)) {
    return(object$tuning_grid)
  }
  model <- object$.model
  if (!is.list(model) || is.null(model$x_mat)) {
    .abort_bad_arg(
      "object",
      "does not contain the required fitted training data; refit before calling plot_tuning()",
      quoted = FALSE
    )
  }
  score_mat <- .plot_tuning_plqr_score(object)
  payload <- .plot_tuning_payload(object, tuning = tuning, se = se)
  suppressWarnings(.plqr_dnn_tuning_grid(
    candidates = payload$grid$candidates,
    ranks = .plot_tuning_ranks(object),
    moment_hat = object$moment,
    se_hat = object$moment_se,
    x_mat = model$x_mat,
    x0_mat = as.matrix(object$x0),
    score_mat = score_mat,
    distance = model$distance,
    method = model$method,
    k = object$meta$k,
    normalize = model$normalize,
    ctx = model$ctx,
    se = se,
    jackknife = .plot_tuning_jackknife(se),
    tuning_payload = payload,
    reuse_default = FALSE
  ))
}

.plot_tuning_panel <- function(object, grid, targets, level, candidates) {
  object_grid <- object
  object_grid$tuning_grid <- grid
  payload <- .plot_candidate_estimates_frame(
    object = object_grid,
    candidates = candidates,
    level = level,
    fn = "plot_tuning"
  )
  df <- payload$data

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
  df$selected <- as.logical(df$is_default)
  df$selected_shape <- factor(
    ifelse(df$selected, "selected", "candidate"),
    levels = c("candidate", "selected")
  )
  if (length(unique(df$coef)) > 1L && length(targets) > 1L) {
    df$facet <- factor(
      paste(df$target, df$coef, sep = " / "),
      levels = unique(paste(df$target, df$coef, sep = " / "))
    )
  } else if (length(unique(df$coef)) > 1L) {
    df$facet <- factor(df$coef, levels = unique(df$coef))
  } else {
    df$facet <- df$target
  }

  has_se <- !all(is.na(df$se))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$s, y = .data$estimate)) +
    ggplot2::geom_vline(xintercept = as.integer(object$meta$s), linetype = "dotted")
  if (has_se) {
    p <- p + ggplot2::geom_pointrange(
      ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper,
                   shape = .data$selected_shape,
                   colour = .data$candidate_label), na.rm = TRUE
    )
  } else {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(shape = .data$selected_shape,
                   colour = .data$candidate_label), na.rm = TRUE
    )
  }
  p <- p +
    ggplot2::scale_shape_manual(
      values = c(candidate = 1, selected = 16), guide = "none"
    ) +
    ggplot2::scale_x_continuous(breaks = sort(unique(df$s))) +
    ggplot2::labs(colour = "candidate")
  if (length(unique(df$facet)) > 1L) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$facet),
                                 scales = "free_y")
  }
  src <- if (!is.null(object$meta$tuning$source)) {
    object$meta$tuning$source
  } else {
    "manual"
  }
  p + ggplot2::labs(
    x = "kernel order s", y = "Estimate",
    caption = paste(
      c(
        sprintf(
          "filled point / dotted line: fitted s = %d (%s); n = %d",
          as.integer(object$meta$s), src, as.integer(object$meta$n)
        ),
        payload$drop_caption,
        .plot_estimate_availability_caption(object, df)
      ),
      collapse = "\n"
    )
  ) + .rdistnn_theme()
}

#' Plot estimate sensitivity to the kernel order
#'
#' Evaluates or reuses the candidate kernel orders of a tuning grid
#' (default: the paper rule's lambda grid via [tuning_dnn_paper()], always
#' including the fitted `s`) and plots the estimates with confidence
#' bounds against `s`, marking the fitted kernel order with a filled point
#' and a dotted guide. For grid fits, the stored `tuning_grid` payload is
#' consumed directly. For scalar fits, the internal grid engine is used on
#' the stored model state; public estimators are not refit.
#'
#' Scalar grid evaluation uses the training data stored on the fit; fits
#' without the required fitted training state must be refitted first.
#'
#' @section Experimental methods:
#' Methods for `plqr_result` objects are experimental, as is all PLQR
#' functionality in this release; see [plqr_dnn()].
#'
#' @param object A fitted `np_result`, `cate_result`, `affine_ratio_result`
#'   (including `clate_result`), or `plqr_result`.
#' @param tuning A tuning specification for the candidate grid. Scalar fits
#'   use their retained tuning specification when omitted. Fits with a manually
#'   chosen `s` use [tuning_dnn_paper()], including that fitted order; optional
#'   return fields do not change the sweep. An explicit specification overrides
#'   that scalar grid. Stored grid fits reuse their existing candidates.
#' @param targets Optional integer vector of target ids (defaults to at
#'   most three: first, middle, last).
#' @param level Confidence level for the bounds (default `0.95`).
#' @param se Logical; compute jackknife bounds for each candidate
#'   when `object` is a scalar fit (default `TRUE`). Stored grid fits reuse
#'   the standard errors already present in `object$tuning_grid`.
#' @param candidates Candidate rows to display. The default `"all"` shows
#'   the tuning sweep; a character vector of candidate ids or `"default"`
#'   can be used to restrict the display.
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
                        targets = NULL, level = 0.95, se = TRUE,
                        candidates = "all", ...) {
  UseMethod("plot_tuning")
}

#' @rdname plot_tuning
#' @export
plot_tuning.default <- function(object, tuning = tuning_dnn_paper(),
                                targets = NULL, level = 0.95, se = TRUE,
                                candidates = "all",
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
                                  candidates = "all",
                                  ...) {
  .plot_require_ggplot2("plot_tuning")
  se <- .validate_scalar_logical(se, "se")
  if (missing(tuning)) tuning <- NULL
  grid <- .plot_tuning_grid_np(object, tuning = tuning, se = se)
  .plot_tuning_panel(object, grid = grid, targets = targets,
                     level = level, candidates = candidates)
}

#' @rdname plot_tuning
#' @export
plot_tuning.cate_result <- function(object, tuning = tuning_dnn_paper(),
                                    targets = NULL, level = 0.95,
                                    se = TRUE, candidates = "all", ...) {
  .plot_require_ggplot2("plot_tuning")
  se <- .validate_scalar_logical(se, "se")
  if (missing(tuning)) tuning <- NULL
  grid <- .plot_tuning_grid_cate(object, tuning = tuning, se = se)
  .plot_tuning_panel(object, grid = grid, targets = targets,
                     level = level, candidates = candidates)
}

#' @rdname plot_tuning
#' @export
plot_tuning.affine_ratio_result <- function(object, tuning = tuning_dnn_paper(),
    targets = NULL, level = .95, se = TRUE, candidates = "all", ...) {
  .plot_require_ggplot2("plot_tuning")
  if (missing(tuning)) tuning <- NULL
  se <- .validate_scalar_logical(se, "se")
  grid <- object$tuning_grid
  if (is.null(grid)) {
    model <- .affine_prediction_state(object)
    payload <- .plot_tuning_payload(object, tuning, se)
    engine <- .plot_tuning_jackknife(se)
    default <- object
    default$meta$jackknife <- engine
    grid <- suppressWarnings(.affine_ratio_tuning_grid(default, model$signal,
      model$geometry$ranks, model, se, engine, model$denominator_tol, payload,
      reuse_default = FALSE))
  }
  .plot_tuning_panel(object, grid, targets, level, candidates)
}

#' @rdname plot_tuning
#' @export
plot_tuning.plqr_result <- function(object, tuning = tuning_dnn_paper(),
                                    targets = NULL, level = 0.95,
                                    se = TRUE, candidates = "all", ...) {
  .plot_require_ggplot2("plot_tuning")
  if (missing(tuning)) tuning <- NULL
  grid <- .plot_tuning_grid_plqr(object, tuning = tuning, se = se)
  .plot_tuning_panel(object, grid = grid, targets = targets,
                     level = level, candidates = candidates)
}
