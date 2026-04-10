# Diagnostics methods for user-level DNN estimator results.

.dnn_diag_top_share <- function(w, k) {
  k_use <- max(1L, min(length(w), as.integer(k)))
  sum(sort(as.numeric(w), decreasing = TRUE)[seq_len(k_use)])
}

.dnn_weighted_quantile <- function(x, w, prob = 0.5) {
  if (!is.numeric(x) || !is.numeric(w) || length(x) != length(w)) {
    return(NA_real_)
  }
  if (.has_na_only(x) || .has_nonfinite_non_na(x)) {
    return(NA_real_)
  }
  if (.has_na_only(w) || .has_nonfinite_non_na(w)) {
    return(NA_real_)
  }

  w_sum <- sum(w)
  if (!is.finite(w_sum) || w_sum <= 0) {
    return(NA_real_)
  }

  ord <- order(x)
  x_ord <- as.numeric(x[ord])
  w_ord <- as.numeric(w[ord])
  cum_w <- cumsum(w_ord) / sum(w_ord)
  idx <- which(cum_w >= as.numeric(prob))[1L]
  if (is.na(idx)) {
    return(x_ord[length(x_ord)])
  }
  x_ord[idx]
}

.dnn_tuning_diag_fields <- function(tuning_payload, m) {
  if (is.null(tuning_payload) || !is.list(tuning_payload)) {
    return(list(
      tuning_source = rep(NA_character_, m),
      tuning_name = rep(NA_character_, m),
      tuning_default_s = rep(NA_integer_, m),
      tuning_grid_available = rep(NA, m)
    ))
  }

  source <- if (!is.null(tuning_payload$source)) {
    as.character(tuning_payload$source)
  } else {
    NA_character_
  }
  name <- if (!is.null(tuning_payload$default$meta$name)) {
    as.character(tuning_payload$default$meta$name)
  } else if (!is.null(source)) {
    source
  } else {
    NA_character_
  }
  default_s <- if (!is.null(tuning_payload$default$s)) {
    as.integer(tuning_payload$default$s)
  } else {
    NA_integer_
  }
  has_grid <- !is.null(tuning_payload$grid$candidates) &&
    is.data.frame(tuning_payload$grid$candidates)

  list(
    tuning_source = rep(source, m),
    tuning_name = rep(name, m),
    tuning_default_s = rep(default_s, m),
    tuning_grid_available = rep(has_grid, m)
  )
}

.dnn_weight_diagnostics <- function(W, meta = NULL, dx = NULL, tuning_payload = NULL) {
  w_mat <- as.matrix(W)
  if (!is.numeric(w_mat)) {
    .abort_bad_arg("W", "must be a numeric matrix", quoted = FALSE)
  }
  if (nrow(w_mat) < 1L || ncol(w_mat) < 1L) {
    .abort_bad_arg("W", "must have positive nrow and ncol", quoted = FALSE)
  }
  .abort_if_na_or_nonfinite(w_mat, arg = "W", quoted = FALSE)

  n <- nrow(w_mat)
  m <- ncol(w_mat)

  meta_n <- if (!is.null(meta$n)) as.integer(meta$n) else as.integer(n)
  meta_p <- if (!is.null(meta$p)) as.integer(meta$p) else NA_integer_
  meta_k <- if (!is.null(meta$k)) as.integer(meta$k) else NA_integer_
  meta_s <- if (!is.null(meta$s)) as.integer(meta$s) else NA_integer_
  meta_method <- if (!is.null(meta$method)) as.character(meta$method) else NA_character_
  meta_distance <- if (!is.null(meta$distance)) as.character(meta$distance) else NA_character_

  n_eff <- 1 / colSums(w_mat^2)
  max_weight <- apply(w_mat, 2L, max)
  top1 <- apply(w_mat, 2L, .dnn_diag_top_share, k = 1L)
  top5 <- apply(w_mat, 2L, .dnn_diag_top_share, k = 5L)
  top10 <- apply(w_mat, 2L, .dnn_diag_top_share, k = 10L)
  n_nonzero <- colSums(abs(w_mat) > sqrt(.Machine$double.eps))

  radius_rank_k <- rep(NA_real_, m)
  radius_rank_s <- rep(NA_real_, m)
  weighted_mean_distance <- rep(NA_real_, m)
  weighted_median_distance <- rep(NA_real_, m)
  distance_available <- rep(FALSE, m)

  if (!is.null(dx)) {
    d_mat <- as.matrix(dx)
    if (!is.numeric(d_mat)) {
      .abort_bad_arg("dx", "must be a numeric matrix", quoted = FALSE)
    }
    if (nrow(d_mat) != nrow(w_mat) || ncol(d_mat) != ncol(w_mat)) {
      .abort_bad_arg(
        "dx",
        "must have the same dimensions as `W`",
        expected = sprintf("nrow=%d, ncol=%d", nrow(w_mat), ncol(w_mat)),
        actual = sprintf("nrow=%d, ncol=%d", nrow(d_mat), ncol(d_mat)),
        quoted = FALSE
      )
    }
    .abort_if_na_or_nonfinite(d_mat, arg = "dx", quoted = FALSE)

    if (!is.na(meta_k) && meta_k >= 1L) {
      k_use <- max(1L, min(n, meta_k))
      radius_rank_k <- apply(d_mat, 2L, function(d) sort(d, partial = k_use)[k_use])
    }
    if (!is.na(meta_s) && meta_s >= 1L) {
      s_use <- max(1L, min(n, meta_s))
      radius_rank_s <- apply(d_mat, 2L, function(d) sort(d, partial = s_use)[s_use])
    }

    denom <- colSums(w_mat)
    good <- is.finite(denom) & denom > 0
    weighted_mean_distance[good] <- colSums(d_mat * w_mat)[good] / denom[good]

    weighted_median_distance <- vapply(
      seq_len(m),
      function(j) .dnn_weighted_quantile(d_mat[, j], w_mat[, j], prob = 0.5),
      numeric(1L)
    )
    distance_available <- rep(TRUE, m)
  }

  tuning <- .dnn_tuning_diag_fields(tuning_payload = tuning_payload, m = m)
  out <- data.frame(
    target_id = paste0("target_", seq_len(m)),
    n = rep(meta_n, m),
    p = rep(meta_p, m),
    k = rep(meta_k, m),
    s = rep(meta_s, m),
    method = rep(meta_method, m),
    distance = rep(meta_distance, m),
    n_eff = as.numeric(n_eff),
    max_weight = as.numeric(max_weight),
    top1_weight_share = as.numeric(top1),
    top5_weight_share = as.numeric(top5),
    top10_weight_share = as.numeric(top10),
    n_nonzero_weight = as.integer(n_nonzero),
    distance_diagnostics_available = as.logical(distance_available),
    radius_rank_k = as.numeric(radius_rank_k),
    radius_rank_s = as.numeric(radius_rank_s),
    weighted_mean_distance = as.numeric(weighted_mean_distance),
    weighted_median_distance = as.numeric(weighted_median_distance),
    tuning_source = tuning$tuning_source,
    tuning_name = tuning$tuning_name,
    tuning_default_s = tuning$tuning_default_s,
    tuning_grid_available = tuning$tuning_grid_available,
    stringsAsFactors = FALSE
  )
  rownames(out) <- out$target_id
  class(out) <- c("dnn_diagnostics", "data.frame")
  out
}

#' Compute diagnostics for fitted localized estimators
#'
#' Returns a compact diagnostics table for supported fit objects.
#'
#' @param object A fitted estimator object.
#' @param ... Not used.
#' @return An object of class `dnn_diagnostics`.
#' @examples
#' set.seed(1)
#' n <- 30
#' x_mat <- matrix(rnorm(n * 2), n, 2)
#' y_vec <- x_mat[, 1] + rnorm(n, sd = 0.2)
#' fit <- np_dnn(X = x_mat, x0 = c(0, 0), Y = y_vec, se = FALSE)
#' diagnostics(fit)
#' @export
diagnostics <- function(object, ...) {
  UseMethod("diagnostics")
}

#' @rdname diagnostics
#' @export
diagnostics.default <- function(object, ...) {
  .abort_bad_arg(
    "object",
    sprintf("diagnostics() is not implemented for class '%s'", class(object)[1L]),
    quoted = FALSE
  )
}

#' @rdname diagnostics
#' @param object An object returned by [np_dnn()].
#' @export
diagnostics.np_result <- function(object, ...) {
  if (is.null(object$W)) {
    .abort_bad_arg(
      "object$W",
      "is required for diagnostics; refit with `return` including \"W\"",
      quoted = FALSE
    )
  }
  tuning_payload <- if (!is.null(object$meta$tuning)) {
    object$meta$tuning
  } else {
    object$tuning
  }
  .dnn_weight_diagnostics(
    W = object$W,
    meta = object$meta,
    dx = object$dx,
    tuning_payload = tuning_payload
  )
}

#' @rdname diagnostics
#' @param x A `dnn_diagnostics` object.
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.dnn_diagnostics <- function(x, ...) {
  cat("<dnn_diagnostics>\n", sep = "")
  cat("  targets: ", nrow(x), "\n", sep = "")
  cat(
    "  n_eff  : min=", round(min(x$n_eff), 4), ", median=", round(stats::median(x$n_eff), 4),
    ", max=", round(max(x$n_eff), 4), "\n",
    sep = ""
  )
  cat(
    "  weights: max(top5 share)=", round(max(x$top5_weight_share), 4),
    ", max(weight)=", round(max(x$max_weight), 4), "\n",
    sep = ""
  )
  invisible(x)
}
