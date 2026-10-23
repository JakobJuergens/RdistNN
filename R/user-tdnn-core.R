# Shared internal helpers for user-level TDNN estimators.
#
# The TDNN wrapper reuses the same public-boundary helpers as the DNN wrapper
# wherever the contract is identical (for example target resolution and
# jackknife variance assembly).  TDNN-specific code starts only where the
# rank-weight construction or jackknife engine truly differs.

.tdnn_weights_with_extras <- function(
  x_mat, x0_mat, distance, method, s1, s2, dim, normalize, ctx, return
) {
  keep_extras <- unique(c("ranks", "dx", intersect(return, "ord")))
  if ("neighbors" %in% return) keep_extras <- unique(c(keep_extras, "ranks"))
  wt <- weights_tdnn(
    X         = x_mat,
    X0        = x0_mat,
    distance  = distance,
    s1        = as.integer(s1),
    s2        = as.integer(s2),
    dim       = as.integer(dim),
    method    = method,
    normalize = normalize,
    ctx       = ctx,
    return    = unique(c("W", keep_extras))
  )
  .dnn_warn_degenerate_weights(wt$W, method = "TDNN")
  list(
    wt = wt,
    W  = wt$W,
    m  = ncol(wt$W)
  )
}

#' Delete-one jackknife SE for TDNN rank-weighted estimators (internal)
#'
#' Mirrors `.dnn_jackknife_se_loo_matrix()` with the DNN rank-weight spec
#' replaced by the TDNN rank-weight spec.
#'
#' **Fast path** (`jackknife = "fast"`): computes full-sample ranks and the
#' deletion weight table once, then uses the shared centered prefix/suffix
#' engine. No distances are recomputed.
#'
#' **LOO fallback** (`jackknife = "loo"`): for each `i`, recomputes weights
#' from scratch on the `(n-1)` reduced dataset.  Strictly slower; useful for
#' validation.
#'
#' @param x_mat    Numeric matrix (n x p).
#' @param x0_mat   Numeric matrix (m x p) of target points.
#' @param signal   Numeric vector/matrix of per-observation signal values
#'   (n x q).
#' @param distance Distance spec or registered name.
#' @param method   `"exact"` or `"asymp"`.
#' @param s1       Integer; smaller U-statistic order.
#' @param s2       Integer; larger U-statistic order.
#' @param dim      Integer; localization dimension.
#' @param normalize Logical; passed to the asymptotic rank-weight constructor.
#' @param ctx      Optional context.
#' @param signal_arg Name of signal argument for error messages.
#' @param jackknife One of `"fast"` (default) or `"loo"`.
#' @param precomputed_ranks Optional precomputed n x m integer rank matrix
#'   (only used by the fast path; if `NULL` ranks are computed internally).
#'
#' @return Numeric matrix (q x m) of jackknife standard errors.
#' @keywords internal
.tdnn_jackknife_se_loo_matrix <- function(
  x_mat, x0_mat, signal, distance, method, s1, s2, dim, normalize, ctx,
  signal_arg         = "signal",
  jackknife          = c("fast", "loo"),
  precomputed_ranks  = NULL, precomputed_groups = NULL
) {
  jackknife <- match.arg(jackknife)
  n <- nrow(x_mat)
  if (as.integer(s2) >= n) {
    stop(
      "Delete-1 jackknife requires `s2 < n` so the same scales can be used after deletion.",
      call. = FALSE
    )
  }

  sig <- .dnn_signal_matrix(signal, n = n, arg = signal_arg)
  q   <- ncol(sig)
  m   <- nrow(x0_mat)

  # ---- fast path: rank-update loop ----------------------------------------
  if (identical(jackknife, "fast")) {
    ranks_full <- precomputed_ranks
    if (is.null(ranks_full) || is.null(precomputed_groups)) {
      ranking_full <- weights_tdnn(
        X         = x_mat,
        X0        = x0_mat,
        distance  = distance,
        s1        = s1,
        s2        = s2,
        dim       = dim,
        method    = method,
        normalize = normalize,
        ctx       = ctx,
        return    = "ranks"
      )
      if (is.null(ranks_full)) ranks_full <- ranking_full$ranks
      precomputed_groups <- ranking_full$groups
    }

    ranks_full <- as.matrix(ranks_full)
    if (nrow(ranks_full) != n || ncol(ranks_full) != m) {
      .abort_bad_state("`precomputed_ranks` must be an n x m ranks matrix.")
    }
    if (!is.integer(ranks_full)) storage.mode(ranks_full) <- "integer"

    rw <- if (identical(method, "exact")) {
      rank_weight_tdnn_exact(s1 = as.integer(s1), s2 = as.integer(s2),
                              dim = as.integer(dim))
    } else {
      rank_weight_tdnn_asymp(s1 = as.integer(s1), s2 = as.integer(s2),
                              dim = as.integer(dim), normalize = normalize)
    }

    loo <- .dnn_loo_centered(sig, ranks_full, .rank_weight_table(rw, n - 1L),
      tab_full = .rank_weight_table(rw, n), groups = precomputed_groups)$loo

    return(.dnn_jackknife_se_from_loo(loo = loo, q = q, n = n, m = m))
  }

  # Recompute each deletion's weights, including its sample-dependent metric.
  weight_fun <- function(i) {
    weights_tdnn(
      X = x_mat[-i, , drop = FALSE], X0 = x0_mat, distance = distance,
      s1 = s1, s2 = s2, dim = dim, method = method,
      normalize = normalize, ctx = ctx, return = "W"
    )$W
  }
  loo <- .jackknife_recomputed_centered(sig, weight_fun)
  .dnn_jackknife_se_from_loo(loo = loo, q = q, n = n, m = m)
}

# ---- Fast path: prefix/suffix jackknife via jackknife_rank_full --------------

#' Delete-one jackknife SE for TDNN via prefix/suffix rank-order engine
#'
#' Computes the jackknife standard error for a TDNN rank-weighted estimator
#' using the shared centered prefix/suffix engine. It precomputes the TDNN
#' rank tables at sizes \eqn{n} and \eqn{n-1} and passes them to
#' `.dnn_loo_centered()`, as does the `fast` route of
#' `.tdnn_jackknife_se_loo_matrix()`.
#'
#' The C++ engine is used automatically when the package has been compiled;
#' otherwise the same centered calculation uses its R fallback.
#'
#' @param x_mat    Numeric matrix (n x p).
#' @param x0_mat   Numeric matrix (m x p) of target points.
#' @param signal   Numeric vector or matrix (n x q) of per-observation signals.
#' @param distance Distance spec or registered name.
#' @param method   `"exact"` or `"asymp"`.
#' @param s1       Integer; smaller U-statistic order.
#' @param s2       Integer; larger U-statistic order.
#' @param dim      Integer; localization dimension.
#' @param normalize Logical; passed to the asymptotic rank-table constructor.
#' @param ctx      Optional context.
#' @param signal_arg Name of the signal argument used in error messages.
#' @param precomputed_ranks Optional precomputed n x m integer rank matrix.
#'   When `NULL` (default) ranks are computed internally via [weights_tdnn()].
#'
#' @return Numeric matrix (q x m) of jackknife standard errors.
#' @keywords internal
.tdnn_jackknife_se_loo_full <- function(
  x_mat, x0_mat, signal, distance, method, s1, s2, dim, normalize, ctx,
  signal_arg        = "signal",
  precomputed_ranks = NULL, precomputed_groups = NULL
) {
  n <- nrow(x_mat)
  if (as.integer(s2) >= n) {
    stop(
      "Delete-1 jackknife requires `s2 < n` so the same scales can be used after deletion.",
      call. = FALSE
    )
  }

  sig <- .dnn_signal_matrix(signal, n = n, arg = signal_arg)
  q   <- ncol(sig)
  m   <- nrow(x0_mat)

  # Ranks
  ranks_full <- precomputed_ranks
  if (is.null(ranks_full) || is.null(precomputed_groups)) {
    ranking_full <- weights_tdnn(
      X         = x_mat,
      X0        = x0_mat,
      distance  = distance,
      s1        = s1,
      s2        = s2,
      dim       = dim,
      method    = method,
      normalize = normalize,
      ctx       = ctx,
      return    = "ranks"
    )
    if (is.null(ranks_full)) ranks_full <- ranking_full$ranks
    precomputed_groups <- ranking_full$groups
  }
  if (!is.integer(ranks_full)) storage.mode(ranks_full) <- "integer"

  # Precompute rank tables at n and n-1
  s1i  <- as.integer(s1)
  s2i  <- as.integer(s2)
  dimi <- as.integer(dim)

  if (identical(method, "exact")) {
    wn   <- tdnn_rank_table_exact(n,       s1i, s2i, dimi)
    wnm1 <- tdnn_rank_table_exact(n - 1L,  s1i, s2i, dimi)
  } else {
    wn   <- tdnn_rank_table_asymp(n,       s1i, s2i, dimi, normalize = normalize)
    wnm1 <- tdnn_rank_table_asymp(n - 1L,  s1i, s2i, dimi, normalize = normalize)
  }

  loo <- .dnn_loo_centered(sig, ranks_full, wnm1, tab_full = wn,
                           groups = precomputed_groups)$loo

  .dnn_jackknife_se_from_loo(loo = loo, q = q, n = n, m = m)
}
