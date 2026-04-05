# High-level API: from {X, X0/x0} + RANKING + rank_weight -> W
# Modular, user-friendly, ranker-agnostic. Works with distances, custom rankers,
# or precomputed ranks/orders/distances. DNN wrappers provided below.
# -----------------------------------------------------------------------------

# tiny internal helper
`%||%` <- function(a, b) if (is.null(a)) b else a

# deterministic tie-breaker (distance, then index)
.order_with_ties <- function(dx) order(dx, seq_along(dx))

# ---- resolvers ---------------------------------------------------------------

# Use your existing distance stack
.resolve_distance <- function(distance) {
  if (is_distance(distance)) {
    return(distance)
  }
  if (is.character(distance)) {
    return(get_distance(distance))
  }
  if (is.function(distance)) {
    return(as_distance(distance, name = "user_distance"))
  }
  .abort_bad_arg(
    "distance",
    "must be a dnn_distance, a registered name, or a function(X, x0, ctx)->numeric(n)"
  )
}

# Accept a variety of "ranking" providers and return ranks (and optionally ord, dx)
# ranking can be:
#   • a distance spec/name -> compute dx (n x m), ord, ranks
#   • a function(X, X0, ctx) -> list with $ranks / $ord / $dx OR a ranks matrix
#   • a list with $ranks (n x m), and/or $ord (n x m), and/or $dx (n x m)
#   • an integer matrix ranks (n x m) directly
# Notes: n = nrow(X), m = nrow(X0). Each column is one target.
.resolve_ranking <- function(ranking, X, X0, ctx = NULL) {
  x_mat <- as.matrix(X)
  x0_mat <- as.matrix(X0)
  n <- nrow(x_mat)
  m <- nrow(x0_mat)

  make_from_dx <- function(dist_matrix) {
    order_matrix <- apply(dist_matrix, 2L, .order_with_ties) # n x m
    rank_matrix <- apply(order_matrix, 2L, ord_to_ranks) # n x m
    list(ranks = rank_matrix, ord = order_matrix, dx = dist_matrix)
  }
  make_from_ord <- function(order_matrix) {
    if (!is.integer(order_matrix)) storage.mode(order_matrix) <- "integer"
    if (nrow(order_matrix) != n || ncol(order_matrix) != m) {
      .abort_bad_arg("ord", "must be n x m")
    }
    for (j in seq_len(m)) {
      o <- order_matrix[, j]
      if (!setequal(o, seq_len(n))) .abort_bad_state("Each column of `ord` must be a permutation of 1..n.")
    }
    rank_matrix <- apply(order_matrix, 2L, ord_to_ranks)
    list(ranks = rank_matrix, ord = order_matrix, dx = NULL)
  }
  make_from_ranks <- function(rank_matrix) {
    if (!is.integer(rank_matrix)) storage.mode(rank_matrix) <- "integer"
    if (nrow(rank_matrix) != n || ncol(rank_matrix) != m) {
      .abort_bad_arg("ranks", "must be n x m")
    }
    for (j in seq_len(m)) {
      r <- rank_matrix[, j]
      if (!setequal(r, seq_len(n))) .abort_bad_state("Each column of `ranks` must be a permutation of 1..n.")
    }
    list(ranks = rank_matrix, ord = NULL, dx = NULL)
  }

  # 1) Distance spec / name -> compute distances
  if (is_distance(ranking) || is.character(ranking)) {
    dist_spec <- .resolve_distance(ranking)
    dist_matrix <- if (!is.null(dist_spec$eval_batch)) {
      distance_eval_batch(dist_spec, x_mat, x0_mat, ctx) # n x m
    } else {
      vapply(seq_len(m), function(j) distance_eval(dist_spec, x_mat, x0_mat[j, ], ctx), numeric(n))
    }
    return(make_from_dx(dist_matrix))
  }

  # 2) Function -> run it and recurse on its output
  if (is.function(ranking)) {
    obj <- ranking(x_mat, x0_mat, ctx)
    return(.resolve_ranking(obj, x_mat, x0_mat, ctx))
  }

  # 3) List: respect whichever of ranks/ord/dx it provides
  if (is.list(ranking)) {
    if (!is.null(ranking$ranks)) {
      return(make_from_ranks(ranking$ranks))
    }
    if (!is.null(ranking$ord)) {
      return(make_from_ord(ranking$ord))
    }
    if (!is.null(ranking$dx)) {
      return(make_from_dx(ranking$dx))
    }
    .abort_bad_state("`ranking` list must contain one of: $ranks, $ord, or $dx.")
  }

  # 4) Matrix -> assume it is a ranks matrix
  if (is.matrix(ranking)) {
    return(make_from_ranks(ranking))
  }

  .abort_bad_state(
    paste(
      "Unsupported `ranking` type. Provide a distance spec/name, a",
      "ranks/ord/dx matrix/list, or a function(X, X0, ctx)->{ranks|ord|dx}."
    )
  )
}

# Accept a rank_weight (or a zero-arg factory returning one)
.resolve_weight <- function(weight) {
  # Accept a rank_weight spec
  if (is_rank_weight(weight)) {
    return(list(kind = "rank", spec = weight))
  }
  # Accept a zero-arg factory that returns a rank_weight
  if (is.function(weight)) {
    w <- weight()
    if (!is_rank_weight(w)) {
      stop("Zero-arg weight() factory must return a rank_weight.", call. = FALSE)
    }
    return(list(kind = "rank", spec = w))
  }

  # Clear error messages / deprecation for any dist-weight-like inputs
  if (inherits(weight, "dist_weight") ||
    (is.list(weight) && identical(weight$.__kind, "dist_weight")) ||
    (is.list(weight) && (is.function(weight$eval) || is.function(weight$eval_batch)))) {
    stop("dist_weight support has been removed from RdistNN. Please use a rank_weight.", call. = FALSE)
  }

  stop("`weight` must be a rank_weight (or a zero-arg factory returning one).", call. = FALSE)
}


# ---- main general API --------------------------------------------------------

#' Compute weight matrix W via an arbitrary ranking + weight mapping
#'
#' This is the general, ranker-agnostic top-layer. It consumes a *ranking provider*
#' (distance, custom function, or precomputed ranks/orders/distances) and a weighting
#' rule (`rank_weight`), then returns the weight matrix `W` (n x m).
#'
#' @param X  numeric matrix (n x p), reference data.
#' @param X0 numeric matrix (m x p), target points (rows). If missing, use `x0`.
#' @param x0 optional numeric length-p single target (convenience for m=1).
#' @param ranking one of:
#'   - a `dnn_distance` spec, a registered distance name (e.g., "l2","l1"),
#'   - a function `(X, X0, ctx) ->` list with `$ranks`/`$ord`/`$dx`, or an `n x m` ranks matrix,
#'   - precomputed object: list with `$ranks`, `$ord`, or `$dx`, or a ranks matrix.
#' @param weight a `rank_weight` spec (or zero-arg factory returning one).
#'   Back-compat: you may also pass `rank_weight=`; it will be used if `weight` is NULL.
#' @param rank_weight (deprecated alias) rank_weight spec or zero-arg factory.
#' @param ctx optional shared context passed to ranking and weight evaluators.
#'   On the ranking side, if `ctx` is `NULL` and a distance spec defines
#'   `setup(X)`, that setup path is used automatically.
#' @param return which extras to include besides `W`: any of c("ord","ranks","dx").
#' @return An object of class `weights_result` with components:
#'   - `W`   : numeric matrix n x m (columns typically sum to 1),
#'   - optional: `ord`, `ranks`, `dx`,
#'   - `meta`: list with `n, m, p, weight_kind, weight_name`.
#' @seealso [weights_dnn()], [rank_weight_eval_batch()], [distance_eval_batch()]
#' @examples
#' set.seed(1)
#' x_mat <- matrix(rnorm(20), 10, 2)
#' x0_mat <- matrix(rnorm(6), 3, 2)
#' rw <- as_rank_weight(function(r) {
#'   w <- 1 / as.numeric(r)
#'   w / sum(w)
#' }, name = "inv_rank")
#' out <- weights_matrix(x_mat, x0_mat, ranking = "l2", weight = rw, return = c("W", "ranks"))
#' dim(out$W)
#' @export
weights_matrix <- function(
  X, X0 = NULL, x0 = NULL,
  ranking,
  weight = NULL,
  rank_weight = NULL,
  ctx = NULL,
  return = c("W")
) {
  # shapes
  x_mat <- .validate_matrix_nonempty(X, arg = "X")
  if (!is.null(X0) && !is.null(x0)) .abort_bad_state("Provide either `X0` or `x0`, not both.")
  if (is.null(X0)) {
    if (is.null(x0)) .abort_bad_state("Provide `X0` (m x p) or a single `x0` (length p).")
    x0_vec <- .validate_target_vector(x0, ncol(x_mat), x0_arg = "x0", x_arg = "X")
    x0_mat <- matrix(x0_vec, nrow = 1L)
  } else {
    x0_mat <- as.matrix(X0)
    if (ncol(x0_mat) == 0L) .abort_bad_state("`X` and `X0` must have at least one column.")
    x0_mat <- .validate_ncol_relation(x_mat, x0_mat, y_arg = "X0", x_arg = "X", relation = "equal")
  }
  if (ncol(x0_mat) == 0L) .abort_bad_state("`X` and `X0` must have at least one column.")

  n <- nrow(x_mat)
  p <- ncol(x_mat)
  m <- nrow(x0_mat)

  # ranking -> ranks (and ord/dx if computable)
  ranking_result <- .resolve_ranking(ranking, x_mat, x0_mat, ctx)

  # resolve weight spec (rank- or distance-based)
  wt_in <- weight %||% rank_weight
  wt_res <- .resolve_weight(wt_in)

  # evaluate weights
  weight_matrix <- {
    # rank-only path
    if (is.null(ranking_result$ranks)) {
      stop("rank_weight provided but no ranks available from `ranking`.", call. = FALSE)
    }
    rank_weight_eval_batch(wt_res$spec, ranking_result$ranks, ctx = ctx)
  }

  out <- list(
    W = weight_matrix,
    meta = list(
      n = n, m = m, p = p,
      weight_kind = wt_res$kind,
      weight_name = wt_res$spec$name
    )
  )
  if ("ord" %in% return && !is.null(ranking_result$ord)) out$ord <- ranking_result$ord
  if ("ranks" %in% return && !is.null(ranking_result$ranks)) out$ranks <- ranking_result$ranks
  if ("dx" %in% return && !is.null(ranking_result$dx)) out$dx <- ranking_result$dx
  class(out) <- "weights_result"
  out
}

#' Print a `weights_result` object
#'
#' @param x A `weights_result` object created by [weights_matrix()].
#' @param ... Not used.
#' @return `x`, invisibly.
#' @examples
#' set.seed(1)
#' x_mat <- matrix(rnorm(20), 10, 2)
#' out <- weights_dnn_exact(x_mat, x0 = c(0, 0), k = 1L, s = 1L)
#' print(out)
#' @export
print.weights_result <- function(x, ...) {
  m <- x$meta$m
  n <- x$meta$n
  p <- x$meta$p
  cat("<weights_result>\n", sep = "")
  cat("  data  : n=", n, ", p=", p, ", m=", m, "\n", sep = "")
  cat("  weight: ", x$meta$weight_kind, " (", x$meta$weight_name, ")\n", sep = "")
  invisible(x)
}

# ---- DNN convenience wrappers ------------------------------------------------
# These are thin wrappers around weights_matrix() using your DNN rank_weight
# constructors and distance ranker. Users who want "just DNN" can call these.

#' DNN weights (exact full-U), distance-based ranking
#'
#' Exact DNN weights are derived from the full rank permutation for each target.
#' They are not constructed by keeping only the first `s` observations. For the
#' exact finite-population rule, the support in rank order is `j <= n - s + k`;
#' in particular, for the paper-faithful `k = 1` DNN this becomes `j <= n - s + 1`.
#'
#' @param X numeric matrix (n x p)
#' @param X0 numeric matrix (m x p); mutually exclusive with \code{x0}
#' @param x0 numeric length-p single target (convenience for m=1)
#' @param distance distance spec or registered name (e.g., "l2")
#' @param k integer; neighbors used in each s-subset (\code{1 <= k <= s})
#' @param s integer; U-statistic order / subset size
#' @param ctx optional context
#' @param return character vector of extras to include; any of \code{c("W","ord","ranks","dx")}
#' @return a \code{weights_result} with \code{W} and metadata
#' @examples
#' set.seed(1)
#' x_mat <- matrix(rnorm(20), 10, 2)
#' x0_mat <- matrix(rnorm(6), 3, 2)
#' out <- weights_dnn_exact(x_mat, x0_mat, distance = "l2", k = 2L, s = 3L, return = "W")
#' dim(out$W)
#' @export
weights_dnn_exact <- function(
  X, X0 = NULL, x0 = NULL, distance = "l2", k, s,
  ctx = NULL, return = c("W")
) {
  weights_matrix(
    X = X, X0 = X0, x0 = x0,
    ranking = distance,
    weight = rank_weight_dnn_exact(k = as.integer(k), s = as.integer(s)),
    ctx = ctx,
    return = return
  )
}

#' DNN weights (asymptotic full-U), distance-based ranking
#'
#' The asymptotic DNN rule is evaluated on the full rank grid and is therefore
#' typically treated as a full-support approximation rather than as a hard
#' finite-support cutoff.
#'
#' @inheritParams weights_dnn_exact
#' @param normalize logical; renormalize columns to sum 1
#' @return A `weights_result` object (same structure as [weights_dnn_exact()]).
#' @examples
#' set.seed(1)
#' x_mat <- matrix(rnorm(20), 10, 2)
#' x0_mat <- matrix(rnorm(4), 2, 2)
#' out <- weights_dnn_asymp(x_mat, x0_mat, k = 2L, s = 3L, normalize = TRUE)
#' dim(out$W)
#' @export
weights_dnn_asymp <- function(
  X, X0 = NULL, x0 = NULL, distance = "l2", k, s,
  normalize = TRUE, ctx = NULL, return = c("W")
) {
  weights_matrix(
    X = X, X0 = X0, x0 = x0,
    ranking = distance,
    weight = rank_weight_dnn_asymp(k = as.integer(k), s = as.integer(s), normalize = normalize),
    ctx = ctx,
    return = return
  )
}

#' DNN weights (choose method), distance-based ranking
#'
#' Use `method = "exact"` for the exact finite-population rank-weight rule and
#' `method = "asymp"` for the asymptotic full-rank approximation.
#'
#' @inheritParams weights_dnn_asymp
#' @param method one of \code{"exact"} or \code{"asymp"}
#' @return A `weights_result` object.
#' @examples
#' set.seed(1)
#' x_mat <- matrix(rnorm(20), 10, 2)
#' x0 <- c(0, 0)
#' out <- weights_dnn(x_mat, x0 = x0, k = 2L, s = 3L, method = "exact")
#' dim(out$W)
#' @export
weights_dnn <- function(
  X, X0 = NULL, x0 = NULL, distance = "l2", k, s,
  method = c("exact", "asymp"), normalize = TRUE, ctx = NULL,
  return = c("W")
) {
  method <- match.arg(method)
  if (method == "exact") {
    weights_dnn_exact(
      X = X, X0 = X0, x0 = x0,
      distance = distance,
      k = k, s = s,
      ctx = ctx,
      return = return
    )
  } else {
    weights_dnn_asymp(
      X = X, X0 = X0, x0 = x0,
      distance = distance,
      k = k, s = s,
      normalize = normalize,
      ctx = ctx,
      return = return
    )
  }
}

# ---- TDNN convenience wrappers -----------------------------------------------
# Thin wrappers around weights_matrix() using the TDNN rank_weight constructors.
# k is fixed to 1 (paper-faithful Demirkaya construction).

#' TDNN weights (exact, k = 1), distance-based ranking
#'
#' Computes Demirkaya two-scale DNN (TDNN) weights using the exact rank-weight
#' table ([tdnn_rank_table_exact()]).  TDNN is a signed linear combination of
#' two single-scale exact DNN tables at orders `s1` and `s2`; the combination
#' coefficients are chosen to cancel the leading bias term in the DNN MSE
#' expansion.
#'
#' **Signed weights.** The returned weight matrix `W` can contain negative
#' entries; columns still sum to 1.  This is by design and reflects the
#' bias-correction structure of TDNN.
#'
#' @param X Numeric matrix (n x p).
#' @param X0 Numeric matrix (m x p); mutually exclusive with \code{x0}.
#' @param x0 Numeric length-p single target (convenience for m = 1).
#' @param distance Distance spec or registered name (e.g., \code{"l2"}).
#' @param s1 Positive integer; smaller U-statistic order (\code{s1 < s2 <= n}).
#' @param s2 Positive integer; larger U-statistic order.
#' @param dim Positive integer; localization dimension.  Defaults to
#'   \code{ncol(X)} when \code{NULL}.
#' @param ctx Optional context passed through to distance and rank-weight
#'   evaluators.
#' @param return Character vector of extras to include; any of
#'   \code{c("W", "ord", "ranks", "dx")}.
#'
#' @return A \code{weights_result} object (same structure as
#'   [weights_dnn_exact()]).
#'
#' @seealso [weights_tdnn_asymp()], [weights_tdnn()], [rank_weight_tdnn_exact()]
#'
#' @examples
#' set.seed(1)
#' X  <- matrix(rnorm(30), 15, 2)
#' x0 <- c(0, 0)
#' out <- weights_tdnn_exact(X, x0 = x0, s1 = 3L, s2 = 6L)
#' colSums(out$W)  # 1
#' any(out$W < 0)  # TRUE
#' @export
weights_tdnn_exact <- function(
  X, X0 = NULL, x0 = NULL, distance = "l2", s1, s2, dim = NULL,
  ctx = NULL, return = c("W")
) {
  if (is.null(dim)) dim <- ncol(X)
  weights_matrix(
    X = X, X0 = X0, x0 = x0,
    ranking = distance,
    weight  = rank_weight_tdnn_exact(s1 = as.integer(s1),
                                      s2 = as.integer(s2),
                                      dim = as.integer(dim)),
    ctx    = ctx,
    return = return
  )
}

#' TDNN weights (asymptotic, k = 1), distance-based ranking
#'
#' Computes Demirkaya two-scale DNN (TDNN) weights using the asymptotic
#' rank-weight table ([tdnn_rank_table_asymp()]).  Full rank support: all
#' `n` observations receive non-zero weight (though some are negative).
#'
#' @inheritParams weights_tdnn_exact
#' @param normalize Logical; renormalize rank table columns to sum 1
#'   (default \code{TRUE}).
#'
#' @return A \code{weights_result} object (same structure as
#'   [weights_dnn_asymp()]).
#'
#' @seealso [weights_tdnn_exact()], [weights_tdnn()], [rank_weight_tdnn_asymp()]
#'
#' @examples
#' set.seed(1)
#' X  <- matrix(rnorm(30), 15, 2)
#' x0 <- c(0, 0)
#' out <- weights_tdnn_asymp(X, x0 = x0, s1 = 3L, s2 = 6L)
#' colSums(out$W)  # 1
#' @export
weights_tdnn_asymp <- function(
  X, X0 = NULL, x0 = NULL, distance = "l2", s1, s2, dim = NULL,
  normalize = TRUE, ctx = NULL, return = c("W")
) {
  if (is.null(dim)) dim <- ncol(X)
  weights_matrix(
    X = X, X0 = X0, x0 = x0,
    ranking = distance,
    weight  = rank_weight_tdnn_asymp(s1 = as.integer(s1),
                                      s2 = as.integer(s2),
                                      dim = as.integer(dim),
                                      normalize = normalize),
    ctx    = ctx,
    return = return
  )
}

#' TDNN weights (choose method, k = 1), distance-based ranking
#'
#' Unified entry point for TDNN weights.  Use \code{method = "exact"} for
#' the exact finite-population rank-weight rule ([weights_tdnn_exact()]) and
#' \code{method = "asymp"} for the asymptotic full-rank approximation
#' ([weights_tdnn_asymp()]).
#'
#' @inheritParams weights_tdnn_asymp
#' @param method One of \code{"exact"} or \code{"asymp"}.
#'
#' @return A \code{weights_result} object.
#'
#' @seealso [weights_tdnn_exact()], [weights_tdnn_asymp()]
#'
#' @examples
#' set.seed(1)
#' X  <- matrix(rnorm(30), 15, 2)
#' x0 <- c(0, 0)
#' out <- weights_tdnn(X, x0 = x0, s1 = 3L, s2 = 6L, method = "exact")
#' colSums(out$W)  # 1
#' @export
weights_tdnn <- function(
  X, X0 = NULL, x0 = NULL, distance = "l2", s1, s2, dim = NULL,
  method = c("exact", "asymp"), normalize = TRUE, ctx = NULL,
  return = c("W")
) {
  method <- match.arg(method)
  if (method == "exact") {
    weights_tdnn_exact(
      X = X, X0 = X0, x0 = x0, distance = distance,
      s1 = s1, s2 = s2, dim = dim,
      ctx = ctx, return = return
    )
  } else {
    weights_tdnn_asymp(
      X = X, X0 = X0, x0 = x0, distance = distance,
      s1 = s1, s2 = s2, dim = dim,
      normalize = normalize, ctx = ctx, return = return
    )
  }
}
