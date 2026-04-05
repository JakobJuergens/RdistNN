# Shared internal helpers for user-level TDNN estimators.
#
# Mirrors the structure of user-dnn-core.R, with DNN's (k, s, method)
# replaced by TDNN's (s1, s2, dim, method).  The jackknife logic is
# structurally identical: precompute full-sample ranks once, then for each
# deleted observation decrement the affected ranks and re-evaluate the TDNN
# rank-weight spec.  No distances are recomputed inside the loop.

.tdnn_resolve_targets <- function(x_mat, X0 = NULL, x0 = NULL) {
  if (!is.null(X0) && !is.null(x0)) {
    stop("Provide either `X0` or `x0`, not both.", call. = FALSE)
  }
  if (is.null(X0)) {
    if (is.null(x0)) {
      stop("Provide `X0` (m x p) or a single `x0` (length p).", call. = FALSE)
    }
    x0_vec <- .validate_target_vector(x0, ncol(x_mat), x0_arg = "x0", x_arg = "X")
    return(matrix(x0_vec, nrow = 1L))
  }
  .validate_ncol_relation(x_mat, X0, y_arg = "X0", x_arg = "X", relation = "equal")
}

.tdnn_weights_with_extras <- function(
  x_mat, x0_mat, distance, method, s1, s2, dim, normalize, ctx, return
) {
  keep_extras <- intersect(return, c("ranks", "ord", "dx"))
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
  list(
    wt = wt,
    W  = wt$W,
    m  = ncol(wt$W)
  )
}

.tdnn_jackknife_se_from_loo <- function(loo, q, n, m) {
  # Identical to .dnn_jackknife_se_from_loo; reproduced here so this file is
  # self-contained and does not depend on the DNN internals.
  if (q == 1L) {
    loo_bar  <- colMeans(loo)
    centered <- sweep(loo, 2L, loo_bar, FUN = "-")
    var_jk   <- ((n - 1) / n) * colSums(centered^2)
    return(matrix(sqrt(pmax(var_jk, 0)), nrow = 1L, ncol = m))
  }

  se <- matrix(NA_real_, nrow = q, ncol = m)
  for (j in seq_len(q)) {
    mj      <- matrix(loo[, j, ], nrow = n, ncol = m)
    mj_bar  <- colMeans(mj)
    centered <- sweep(mj, 2L, mj_bar, FUN = "-")
    var_jk  <- ((n - 1) / n) * colSums(centered^2)
    se[j, ] <- sqrt(pmax(var_jk, 0))
  }
  se
}

#' Delete-one jackknife SE for TDNN rank-weighted estimators (internal)
#'
#' Mirrors `.dnn_jackknife_se_loo_matrix()` with the DNN rank-weight spec
#' replaced by the TDNN rank-weight spec.
#'
#' **Fast path** (`jackknife = "fast"`): computes full-sample ranks once, then
#' for each deleted observation `i` adjusts ranks by decrementing entries that
#' were ranked below `i` and re-evaluates the TDNN rank-weight spec on the
#' resulting `(n-1) x m` rank matrix.  No distances are recomputed.
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
  precomputed_ranks  = NULL
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
    if (is.null(ranks_full)) {
      ranks_full <- weights_tdnn(
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
      )$ranks
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

    if (q == 1L) {
      loo <- matrix(NA_real_, nrow = n, ncol = m)
    } else {
      loo <- array(NA_real_, dim = c(n, q, m))
    }

    idx_all <- seq_len(n)
    for (i in idx_all) {
      keep <- (idx_all != i)

      # Decrement ranks of observations that were ranked below observation i
      r_i          <- matrix(ranks_full[i, ], nrow = n - 1L, ncol = m, byrow = TRUE)
      rank_matrix_i <- ranks_full[keep, , drop = FALSE]
      rank_matrix_i <- rank_matrix_i - (rank_matrix_i > r_i)
      if (!is.integer(rank_matrix_i)) storage.mode(rank_matrix_i) <- "integer"

      W_i   <- rank_weight_eval_batch(rw, rank_matrix_i)
      est_i <- crossprod(sig[keep, , drop = FALSE], W_i)
      if (q == 1L) {
        loo[i, ] <- as.numeric(est_i[1L, ])
      } else {
        loo[i, , ] <- est_i
      }
    }

    return(.tdnn_jackknife_se_from_loo(loo = loo, q = q, n = n, m = m))
  }

  # ---- loo fallback: recompute weights from scratch for each deletion ------
  if (q == 1L) {
    loo <- matrix(NA_real_, nrow = n, ncol = m)
    for (i in seq_len(n)) {
      keep  <- (seq_len(n) != i)
      wt_i  <- weights_tdnn(
        X         = x_mat[keep, , drop = FALSE],
        X0        = x0_mat,
        distance  = distance,
        s1        = s1,
        s2        = s2,
        dim       = dim,
        method    = method,
        normalize = normalize,
        ctx       = ctx,
        return    = "W"
      )
      loo[i, ] <- as.numeric(crossprod(sig[keep, 1L], wt_i$W))
    }
    return(.tdnn_jackknife_se_from_loo(loo = loo, q = q, n = n, m = m))
  }

  loo <- array(NA_real_, dim = c(n, q, m))
  for (i in seq_len(n)) {
    keep  <- (seq_len(n) != i)
    wt_i  <- weights_tdnn(
      X         = x_mat[keep, , drop = FALSE],
      X0        = x0_mat,
      distance  = distance,
      s1        = s1,
      s2        = s2,
      dim       = dim,
      method    = method,
      normalize = normalize,
      ctx       = ctx,
      return    = "W"
    )
    loo[i, , ] <- crossprod(sig[keep, , drop = FALSE], wt_i$W)
  }
  .tdnn_jackknife_se_from_loo(loo = loo, q = q, n = n, m = m)
}

# ---- Fast path: prefix/suffix jackknife via jackknife_rank_full --------------

#' Delete-one jackknife SE for TDNN via prefix/suffix rank-order engine
#'
#' Computes the jackknife standard error for a TDNN rank-weighted estimator
#' using the generic \eqn{O(n \cdot m \cdot q)} prefix/suffix engine
#' ([jackknife_rank_full()]).  Unlike `.tdnn_jackknife_se_loo_matrix()`, this
#' path does **not** re-evaluate the rank-weight spec for each deletion; it
#' precomputes the TDNN rank tables at sizes \eqn{n} and \eqn{n-1} once and
#' passes them directly to the engine.
#'
#' The C++ engine is used automatically when the package has been compiled;
#' otherwise [jackknife_rank_full()] falls back to its pure-R implementation.
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
  precomputed_ranks = NULL
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
  if (is.null(ranks_full)) {
    ranks_full <- weights_tdnn(
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
    )$ranks
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

  res <- jackknife_rank_full(ranks_full, sig, wn, wnm1)

  # theta_loo has dim c(n, q, m); squeeze q-dimension for q == 1
  loo <- res$theta_loo
  if (q == 1L) loo <- matrix(loo, n, m)

  .tdnn_jackknife_se_from_loo(loo = loo, q = q, n = n, m = m)
}
