# R/api_weights.R
# -----------------------------------------------------------------------------
# High-level API: from {X, X0/x0} + RANKING + (rank_weight | dist_weight) -> W
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
  stop("`distance` must be a dnn_distance, a registered name, or a function(X, x0, ctx)->numeric(n).")
}

# Accept a variety of "ranking" providers and return ranks (and optionally ord, dx)
# ranking can be:
#   • a distance spec/name -> compute dx (n x m), ord, ranks
#   • a function(X, X0, ctx) -> list with $ranks / $ord / $dx OR a ranks matrix
#   • a list with $ranks (n x m), and/or $ord (n x m), and/or $dx (n x m)
#   • an integer matrix ranks (n x m) directly
# Notes: n = nrow(X), m = nrow(X0). Each column is one target.
.resolve_ranking <- function(ranking, X, X0, ctx = NULL) {
  X <- as.matrix(X)
  X0 <- as.matrix(X0)
  n <- nrow(X)
  m <- nrow(X0)

  make_from_dx <- function(DX) {
    ORD <- apply(DX, 2L, .order_with_ties) # n x m
    RKS <- apply(ORD, 2L, ord_to_ranks) # n x m
    list(ranks = RKS, ord = ORD, dx = DX)
  }
  make_from_ord <- function(ORD) {
    if (!is.integer(ORD)) storage.mode(ORD) <- "integer"
    if (nrow(ORD) != n || ncol(ORD) != m) stop("`ord` must be n x m.")
    for (j in seq_len(m)) {
      o <- ORD[, j]
      if (!setequal(o, seq_len(n))) stop("Each column of `ord` must be a permutation of 1..n.")
    }
    RKS <- apply(ORD, 2L, ord_to_ranks)
    list(ranks = RKS, ord = ORD, dx = NULL)
  }
  make_from_ranks <- function(RKS) {
    if (!is.integer(RKS)) storage.mode(RKS) <- "integer"
    if (nrow(RKS) != n || ncol(RKS) != m) stop("`ranks` must be n x m.")
    for (j in seq_len(m)) {
      r <- RKS[, j]
      if (!setequal(r, seq_len(n))) stop("Each column of `ranks` must be a permutation of 1..n.")
    }
    list(ranks = RKS, ord = NULL, dx = NULL)
  }

  # 1) Distance spec / name -> compute distances
  if (is_distance(ranking) || is.character(ranking)) {
    dist_spec <- .resolve_distance(ranking)
    DX <- if (!is.null(dist_spec$eval_batch)) {
      distance_eval_batch(dist_spec, X, X0, ctx) # n x m
    } else {
      vapply(seq_len(m), function(j) distance_eval(dist_spec, X, X0[j, ], ctx), numeric(n))
    }
    return(make_from_dx(DX))
  }

  # 2) Function -> run it and recurse on its output
  if (is.function(ranking)) {
    obj <- ranking(X, X0, ctx)
    return(.resolve_ranking(obj, X, X0, ctx))
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
    stop("`ranking` list must contain one of: $ranks, $ord, or $dx.")
  }

  # 4) Matrix -> assume it is a ranks matrix
  if (is.matrix(ranking)) {
    return(make_from_ranks(ranking))
  }

  stop("Unsupported `ranking` type. Provide a distance spec/name, a ranks/ord/dx matrix/list, or a function(X, X0, ctx)->{ranks|ord|dx}.")
}

# Accept a rank_weight or a dist_weight (or a zero-arg factory returning one)
.resolve_weight <- function(weight) {
  if (is_rank_weight(weight)) {
    return(list(kind = "rank", spec = weight))
  }
  if (is_dist_weight(weight)) {
    return(list(kind = "dist", spec = weight))
  }
  if (is.function(weight)) {
    w <- weight()
    return(.resolve_weight(w))
  }
  stop("`weight` must be a rank_weight or dist_weight (or a zero-arg factory returning one).")
}

# ---- main general API --------------------------------------------------------

#' Compute weight matrix W via an arbitrary ranking + weight mapping
#'
#' This is the general, ranker-agnostic top-layer. It consumes a *ranking provider*
#' (distance, custom function, or precomputed ranks/orders/distances) and a weighting
#' rule (either a `rank_weight` or a `dist_weight`), then returns the weight matrix
#' `W` (n x m).
#'
#' @param X  numeric matrix (n x p), reference data.
#' @param X0 numeric matrix (m x p), target points (rows). If missing, use `x0`.
#' @param x0 optional numeric length-p single target (convenience for m=1).
#' @param ranking one of:
#'   - a `dnn_distance` spec, a registered distance name (e.g., "l2","l1"),
#'   - a function `(X, X0, ctx) ->` list with `$ranks`/`$ord`/`$dx`, or an `n x m` ranks matrix,
#'   - precomputed object: list with `$ranks`, `$ord`, or `$dx`, or a ranks matrix.
#' @param weight a `rank_weight` **or** `dist_weight` spec (or zero-arg factory).
#'   Back-compat: you may also pass `rank_weight=`; it will be used if `weight` is NULL.
#' @param rank_weight (deprecated alias) rank_weight spec or zero-arg factory.
#' @param ctx optional context for rankers that support `setup`.
#' @param return which extras to include besides `W`: any of c("ord","ranks","dx").
#' @return An object of class `weights_result` with components:
#'   - `W`   : numeric matrix n x m (columns typically sum to 1),
#'   - optional: `ord`, `ranks`, `dx`,
#'   - `meta`: list with `n, m, p, weight_kind, weight_name`.
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
  X <- as.matrix(X)
  if (!is.null(X0) && !is.null(x0)) stop("Provide either `X0` or `x0`, not both.")
  if (is.null(X0)) {
    if (is.null(x0)) stop("Provide `X0` (m x p) or a single `x0` (length p).")
    X0 <- matrix(as.numeric(x0), nrow = 1L)
  } else {
    X0 <- as.matrix(X0)
  }
  if (ncol(X) == 0L || ncol(X0) == 0L) stop("`X` and `X0` must have at least one column.")
  if (ncol(X) != ncol(X0)) stop("ncol(X0) must equal ncol(X).")

  n <- nrow(X)
  p <- ncol(X)
  m <- nrow(X0)

  # ranking -> ranks (and ord/dx if computable)
  rk <- .resolve_ranking(ranking, X, X0, ctx)

  # resolve weight spec (rank- or distance-based)
  wt_in <- weight %||% rank_weight
  wt_res <- .resolve_weight(wt_in)

  # evaluate weights
  W <- switch(wt_res$kind,
    "rank" = {
      if (is.null(rk$ranks)) {
        stop("rank_weight provided but no ranks available from `ranking`.")
      }
      rank_weight_eval_batch(wt_res$spec, rk$ranks)
    },
    "dist" = {
      if (is.null(rk$dx)) {
        stop("dist_weight provided but `ranking` did not supply distances `dx`.")
      }
      dist_weight_eval_batch(wt_res$spec, rk$dx)
    },
    stop("Internal: unknown weight kind.")
  )

  out <- list(
    W = W,
    meta = list(
      n = n, m = m, p = p,
      weight_kind = wt_res$kind,
      weight_name = wt_res$spec$name
    )
  )
  if ("ord" %in% return && !is.null(rk$ord)) out$ord <- rk$ord
  if ("ranks" %in% return && !is.null(rk$ranks)) out$ranks <- rk$ranks
  if ("dx" %in% return && !is.null(rk$dx)) out$dx <- rk$dx
  class(out) <- "weights_result"
  out
}

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
#' @param X numeric matrix (n x p)
#' @param X0 numeric matrix (m x p); mutually exclusive with \code{x0}
#' @param x0 numeric length-p single target (convenience for m=1)
#' @param distance distance spec or registered name (e.g., "l2")
#' @param k integer; neighbors used in each s-subset (\code{1 <= k <= s})
#' @param s integer; U-statistic order / subset size
#' @param ctx optional context
#' @param return character vector of extras to include; any of \code{c("W","ord","ranks","dx")}
#' @return a \code{weights_result} with \code{W} and metadata
#' @export
weights_dnn_exact <- function(X, X0 = NULL, x0 = NULL, distance = "l2", k, s, ctx = NULL, return = c("W")) {
  weights_matrix(
    X, X0, x0,
    ranking = distance,
    weight = rank_weight_dnn_exact(k = as.integer(k), s = as.integer(s)),
    ctx = ctx,
    return = return
  )
}

#' DNN weights (asymptotic full-U), distance-based ranking
#' @inheritParams weights_dnn_exact
#' @param normalize logical; renormalize columns to sum 1
#' @export
weights_dnn_asymp <- function(X, X0 = NULL, x0 = NULL, distance = "l2", k, s, normalize = TRUE, ctx = NULL, return = c("W")) {
  weights_matrix(
    X, X0, x0,
    ranking = distance,
    weight = rank_weight_dnn_asymp(k = as.integer(k), s = as.integer(s), normalize = normalize),
    ctx = ctx,
    return = return
  )
}

#' DNN weights (choose method), distance-based ranking
#' @inheritParams weights_dnn_asymp
#' @param method one of \code{"exact"} or \code{"asymp"}
#' @export
weights_dnn <- function(X, X0 = NULL, x0 = NULL, distance = "l2", k, s, method = c("exact", "asymp"), normalize = TRUE, ctx = NULL, return = c("W")) {
  method <- match.arg(method)
  if (method == "exact") {
    weights_dnn_exact(X, X0, x0, distance, k, s, ctx, return)
  } else {
    weights_dnn_asymp(X, X0, x0, distance, k, s, normalize, ctx, return)
  }
}
