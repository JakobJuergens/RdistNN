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

.dnn_warn_degenerate_weights <- function(W, method = "DNN", concentration = NULL) {
  n <- nrow(W)
  ess <- if (is.null(concentration)) 1 / colSums(W^2) else concentration
  if (is.null(n) || n < 2L) return(invisible(ess))
  affected <- which(ess >= 0.9 * n)
  if (length(affected)) {
    message <- if (identical(method, "TDNN")) {
      sprintf("TDNN signed weights have high effective sample size at %d target(s) (maximum %.1f of n = %d).",
        length(affected), max(ess), n)
    } else {
      sprintf(paste0("DNN weights are (near-)uniform at %d target(s) ",
        "(maximum effective sample size %.1f of n = %d), providing little or no localization."),
        length(affected), max(ess), n)
    }
    hint <- if (identical(method, "TDNN")) {
      paste0("Inspect signed weights and distance diagnostics; ",
        "effective sample size alone does not establish geometric localization.")
    } else {
      paste0("Increase `s` if it is too small relative to `k`, or omit `s` and `tuning` for default tuning. ",
        "A larger `s` cannot distinguish observations at equal computed distance.")
    }
    .dnn_warn("weak_localization", message, affected = length(affected), examples = affected, hint = hint)
  }
  invisible(ess)
}

# Descriptive facts, not warning-history flags. Group preparation/validation
# already supplies target ties; weight concentration is computed while W exists.
.dnn_attach_notice_facts <- function(meta, groups, concentration, facts = NULL) {
  if (is.null(facts)) facts <- .dnn_group_facts(groups)
  m <- length(concentration)
  tied <- if (is.null(groups)) rep(FALSE, m) else unname(facts$tied_targets)
  meta$target_notices <- data.frame(target_id = seq_len(m), distance_ties = tied,
    weight_concentration = unname(concentration),
    weak_localization = unname(meta$n >= 2L & concentration >= .9 * meta$n))
  meta$n_tied_targets <- sum(tied)
  meta$n_weak_targets <- sum(meta$target_notices$weak_localization)
  meta
}

.dnn_weights_with_extras <- function(
  x_mat, x0_mat, distance, method, k, s, normalize, ctx, return
) {
  keep_extras <- unique(c("ranks", "dx", intersect(return, "ord")))
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
  .dnn_warn_degenerate_weights(wt$W,
    concentration = wt$meta$target_notices$weight_concentration)
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
    .abort_bad_arg(arg, "must have nrow(X) rows", expected = n, actual = nrow(sig))
  }
  if (!is.numeric(sig)) {
    .abort_bad_arg(arg, "must be numeric")
  }
  .abort_if_na_or_nonfinite(sig, arg = arg)
  sig
}

.dnn_estimate_from_signal <- function(signal, W, n, arg = "signal") {
  sig <- .dnn_signal_matrix(signal, n = n, arg = arg)
  estimate <- .compensated_crossprod(sig, W)
  if (any(!is.finite(estimate))) {
    targets <- which(colSums(!is.finite(estimate)) > 0L)
    .dnn_abort("numerical_range", sprintf("`%s` produces nonfinite weighted estimates at target(s) %s.",
      .dnn_brief(arg), .dnn_message_examples(targets)$text), argument = arg,
      hint = "Rescale the signal before fitting or predicting.")
  }
  estimate
}

.dnn_attach_common_extras <- function(out, return, wt, W, k, tuning_payload = NULL) {
  concentration <- wt$meta$target_notices$weight_concentration
  if (!is.numeric(concentration) || length(concentration) != ncol(W) ||
      anyNA(concentration) || any(concentration < 0)) {
    .dnn_abort("internal_invariant", "Prepared weights lack current target concentration facts.",
      hint = "Report a reproducible example of this internal inconsistency.")
  }
  out$meta <- .dnn_attach_notice_facts(out$meta, wt$groups, concentration)
  if ("W" %in% return) out$W <- W
  if ("ranks" %in% return && !is.null(wt$ranks)) out$ranks <- wt$ranks
  if ("ord" %in% return && !is.null(wt$ord)) out$ord <- wt$ord
  if ("dx" %in% return && !is.null(wt$dx)) out$dx <- wt$dx
  if ("tuning" %in% return && !is.null(tuning_payload)) out$tuning <- tuning_payload
  if ("ess" %in% return) out$ess <- concentration
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
