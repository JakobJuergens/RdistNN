# Shared internal helpers for user-level DNN estimators.

.dnn_resolve_data_cols <- function(cols, data, arg, allow_multiple = TRUE) {
  if (is.null(cols)) {
    .abort_bad_arg(arg, "must be provided")
  }

  if (is.character(cols)) {
    if (length(cols) == 0L || anyNA(cols) || any(!nzchar(cols))) {
      .abort_bad_arg(arg, "must be a non-empty character vector of column names")
    }
    missing <- setdiff(cols, names(data))
    if (length(missing) > 0L) {
      stop(sprintf(
        "Column(s) %s (specified by `%s`) not found in `data`.\nAvailable columns: %s",
        paste(sprintf("'%s'", missing), collapse = ", "),
        arg,
        paste(names(data), collapse = ", ")
      ), call. = FALSE)
    }
    resolved <- cols
  } else if (is.numeric(cols)) {
    if (length(cols) == 0L) {
      .abort_bad_arg(arg, "must be a non-empty integer vector of column indices")
    }
    if (.has_na_only(cols)) {
      .abort_bad_arg(arg, "cannot contain `NA` indices")
    }
    if (.has_nonfinite_non_na(cols)) {
      .abort_bad_arg(arg, "cannot contain non-finite indices (`Inf`, `-Inf`, or `NaN`)")
    }
    if (any(cols != as.integer(cols))) {
      .abort_bad_arg(arg, "must be a non-empty integer vector of column indices")
    }
    idx <- as.integer(cols)
    p <- ncol(data)
    if (any(idx < 1L | idx > p)) {
      .abort_bad_arg(arg, sprintf("contains out-of-range indices; valid range is 1..%d", p))
    }
    resolved <- names(data)[idx]
  } else {
    .abort_bad_arg(arg, "must be a character vector of column names or an integer vector of indices")
  }

  if (!isTRUE(allow_multiple) && length(resolved) != 1L) {
    .abort_bad_arg(arg, "must reference exactly one column")
  }
  resolved
}

# Extract and validate a single data-frame column by name or index.
.dnn_col_vec <- function(col, data, n, arg, ref = "data") {
  col_name <- .dnn_resolve_data_cols(col, data = data, arg = arg, allow_multiple = FALSE)
  .validate_vector_nrow(data[[col_name]], n, arg = arg, ref = ref, strict = TRUE)
}

# Extract and validate multiple data-frame columns by names/indices.
.dnn_col_mat <- function(cols, data, arg, out_arg = "X") {
  col_names <- .dnn_resolve_data_cols(cols, data = data, arg = arg, allow_multiple = TRUE)
  .validate_matrix_nonempty(data[, col_names, drop = FALSE], arg = out_arg)
}

.dnn_resolve_targets <- function(x_mat, X0 = NULL, x0 = NULL) {
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

.dnn_distance_name <- function(distance) {
  if (is.character(distance)) {
    return(distance)
  }
  if (is.list(distance) && !is.null(distance$name)) {
    return(distance$name)
  }
  "user_distance"
}

.dnn_weights_with_extras <- function(
  x_mat, x0_mat, distance, method, k, s, normalize, ctx, return
) {
  keep_extras <- intersect(return, c("ranks", "ord", "dx"))
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
  crossprod(sig, W)
}

.dnn_attach_common_extras <- function(out, return, wt, W, k) {
  if ("W" %in% return) out$W <- W
  if ("ranks" %in% return && !is.null(wt$ranks)) out$ranks <- wt$ranks
  if ("ord" %in% return && !is.null(wt$ord)) out$ord <- wt$ord
  if ("dx" %in% return && !is.null(wt$dx)) out$dx <- wt$dx
  if ("ess" %in% return) out$ess <- 1 / colSums(W^2)
  if ("neighbors" %in% return) {
    k_eff <- max(1L, min(nrow(W), as.integer(k)))
    if (!is.null(wt$ord)) {
      out$neighbors <- wt$ord[seq_len(k_eff), , drop = FALSE]
    } else if (!is.null(wt$ranks)) {
      out$neighbors <- apply(wt$ranks, 2L, function(r) order(r)[seq_len(k_eff)])
    }
  }
  out
}

.dnn_jackknife_se_from_loo <- function(loo, q, n, m) {
  if (q == 1L) {
    loo_bar <- colMeans(loo)
    centered <- sweep(loo, 2L, loo_bar, FUN = "-")
    var_jk <- ((n - 1) / n) * colSums(centered^2)
    return(matrix(sqrt(pmax(var_jk, 0)), nrow = 1L, ncol = m))
  }

  se <- matrix(NA_real_, nrow = q, ncol = m)
  for (j in seq_len(q)) {
    mj <- matrix(loo[, j, ], nrow = n, ncol = m)
    mj_bar <- colMeans(mj)
    centered <- sweep(mj, 2L, mj_bar, FUN = "-")
    var_jk <- ((n - 1) / n) * colSums(centered^2)
    se[j, ] <- sqrt(pmax(var_jk, 0))
  }
  se
}

.dnn_jackknife_se_loo_matrix <- function(
  x_mat, x0_mat, signal, distance, method, k, s, normalize, ctx,
  signal_arg = "signal",
  jackknife = c("fast", "loo"),
  precomputed_ranks = NULL
) {
  jackknife <- match.arg(jackknife)
  n <- nrow(x_mat)
  if (as.integer(s) >= n) {
    stop(sprintf(
      "Delete-1 jackknife requires s < n (got s=%d, n=%d). Use `se = FALSE` or reduce `s`.",
      as.integer(s), n
    ), call. = FALSE)
  }

  sig <- .dnn_signal_matrix(signal, n = n, arg = signal_arg)
  q <- ncol(sig)
  m <- nrow(x0_mat)

  if (identical(jackknife, "fast")) {
    ranks_full <- precomputed_ranks
    if (is.null(ranks_full)) {
      ranks_full <- weights_dnn(
        X = x_mat,
        X0 = x0_mat,
        distance = distance,
        k = k,
        s = s,
        method = method,
        normalize = normalize,
        ctx = ctx,
        return = "ranks"
      )$ranks
    }

    ranks_full <- as.matrix(ranks_full)
    if (nrow(ranks_full) != n || ncol(ranks_full) != m) {
      .abort_bad_state("`precomputed_ranks` must be an n x m ranks matrix.")
    }
    if (!is.integer(ranks_full)) storage.mode(ranks_full) <- "integer"

    rw <- if (identical(method, "exact")) {
      rank_weight_dnn_exact(k = as.integer(k), s = as.integer(s))
    } else {
      rank_weight_dnn_asymp(k = as.integer(k), s = as.integer(s), normalize = normalize)
    }

    if (q == 1L) {
      loo <- matrix(NA_real_, nrow = n, ncol = m)
    } else {
      loo <- array(NA_real_, dim = c(n, q, m))
    }

    idx_all <- seq_len(n)
    for (i in idx_all) {
      keep <- (idx_all != i)

      r_i <- matrix(ranks_full[i, ], nrow = n - 1L, ncol = m, byrow = TRUE)
      rank_matrix_i <- ranks_full[keep, , drop = FALSE]
      rank_matrix_i <- rank_matrix_i - (rank_matrix_i > r_i)
      if (!is.integer(rank_matrix_i)) storage.mode(rank_matrix_i) <- "integer"

      W_i <- rank_weight_eval_batch(rw, rank_matrix_i)
      est_i <- crossprod(sig[keep, , drop = FALSE], W_i)
      if (q == 1L) {
        loo[i, ] <- as.numeric(est_i[1L, ])
      } else {
        loo[i, , ] <- est_i
      }
    }

    return(.dnn_jackknife_se_from_loo(loo = loo, q = q, n = n, m = m))
  }

  if (q == 1L) {
    loo <- matrix(NA_real_, nrow = n, ncol = m)
    for (i in seq_len(n)) {
      keep <- (seq_len(n) != i)
      wt_i <- weights_dnn(
        X = x_mat[keep, , drop = FALSE],
        X0 = x0_mat,
        distance = distance,
        k = k,
        s = s,
        method = method,
        normalize = normalize,
        ctx = ctx,
        return = "W"
      )
      loo[i, ] <- as.numeric(crossprod(sig[keep, 1L], wt_i$W))
    }
    return(.dnn_jackknife_se_from_loo(loo = loo, q = q, n = n, m = m))
  }

  loo <- array(NA_real_, dim = c(n, q, m))
  for (i in seq_len(n)) {
    keep <- (seq_len(n) != i)
    wt_i <- weights_dnn(
      X = x_mat[keep, , drop = FALSE],
      X0 = x0_mat,
      distance = distance,
      k = k,
      s = s,
      method = method,
      normalize = normalize,
      ctx = ctx,
      return = "W"
    )
    loo[i, , ] <- crossprod(sig[keep, , drop = FALSE], wt_i$W)
  }
  .dnn_jackknife_se_from_loo(loo = loo, q = q, n = n, m = m)
}
