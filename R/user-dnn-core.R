# Common DNN weights, signal evaluation and additive result fields.

.dnn_distance_name <- function(distance) {
  if (is.character(distance)) {
    return(distance)
  }
  if (is.list(distance) && !is.null(distance$name)) {
    return(distance$name)
  }
  "user_distance"
}

.dnn_warn_degenerate_weights <- function(W) {
  n <- nrow(W)
  if (is.null(n) || n < 2L) {
    return(invisible(NULL))
  }
  ess_max <- max(1 / colSums(W^2))
  if (ess_max >= 0.9 * n) {
    warning(
      sprintf(
        paste0(
          "DNN weights are (near-)uniform at one or more target points ",
          "(max effective sample size %.1f of n = %d), so the fit provides ",
          "little or no localization. This typically means `s` is too small ",
          "relative to `k` (for example s = 1, or k = s). Increase `s`, or ",
          "omit both `s` and `tuning` to use the default paper tuning rule."
        ),
        ess_max, n
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

.dnn_weights_with_extras <- function(
  x_mat, x0_mat, distance, method, k, s, normalize, ctx, return
) {
  keep_extras <- intersect(return, c("ranks", "ord", "dx"))
  if ("neighbors" %in% return) keep_extras <- unique(c(keep_extras, "ranks"))
  wt <- weights_dnn(
    X = x_mat,
    X0 = x0_mat,
    distance = distance,
    k = as.integer(k),
    s = as.integer(s),
    method = method,
    normalize = normalize,
    ctx = ctx,
    return = unique(c("W", keep_extras))
  )
  .dnn_warn_degenerate_weights(wt$W)
  list(
    wt = wt,
    W = wt$W,
    m = ncol(wt$W)
  )
}

.dnn_signal_matrix <- function(signal, n, arg = "signal") {
  if (is.null(dim(signal))) {
    signal <- .validate_vector_nrow(signal, n, arg = arg, ref = "X", finite = TRUE)
    return(matrix(as.numeric(signal), ncol = 1L))
  }

  sig <- as.matrix(signal)
  if (nrow(sig) != n) {
    .abort_bad_state(sprintf("nrow(%s) must equal nrow(X).", arg))
  }
  if (!is.numeric(sig)) {
    .abort_bad_arg(arg, "must be numeric")
  }
  .abort_if_na_or_nonfinite(sig, arg = arg)
  sig
}

.dnn_estimate_from_signal <- function(signal, W, n, arg = "signal") {
  sig <- .dnn_signal_matrix(signal, n = n, arg = arg)
  .compensated_crossprod(sig, W)
}

.dnn_attach_common_extras <- function(out, return, wt, W, k, tuning_payload = NULL) {
  if ("W" %in% return) out$W <- W
  if ("ranks" %in% return && !is.null(wt$ranks)) out$ranks <- wt$ranks
  if ("ord" %in% return && !is.null(wt$ord)) out$ord <- wt$ord
  if ("dx" %in% return && !is.null(wt$dx)) out$dx <- wt$dx
  if ("tuning" %in% return && !is.null(tuning_payload)) out$tuning <- tuning_payload
  if ("ess" %in% return) out$ess <- 1 / colSums(W^2)
  if ("neighbors" %in% return) {
    k_eff <- max(1L, min(nrow(W), as.integer(k)))
    if (!is.null(wt$ord)) {
      out$neighbors <- wt$ord[seq_len(k_eff), , drop = FALSE]
    } else if (!is.null(wt$ranks)) {
      out$neighbors <- apply(wt$ranks, 2L, function(r) order(r)[seq_len(k_eff)])
      if (is.null(dim(out$neighbors))) {
        out$neighbors <- matrix(out$neighbors, nrow = k_eff)
      }
    }
  }
  out
}

.dnn_rank_weight_spec <- function(method, k, s, normalize = TRUE) {
  if (identical(method, "exact")) {
    rank_weight_dnn_exact(k = as.integer(k), s = as.integer(s))
  } else {
    rank_weight_dnn_asymp(k = as.integer(k), s = as.integer(s), normalize = normalize)
  }
}
