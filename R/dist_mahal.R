# Mahalanobis geometry is fitted once from the centered training design.

# The SVD acts on the design rather than its Gram matrix, whose condition
# number is squared. A binary scale prevents covariance-range failures.
.mahal_fit <- function(X) {
  if (ncol(X) < 1L)
    stop("distance_mahalanobis() requires at least one covariate.", call. = FALSE)
  if (nrow(X) < 2L)
    stop("distance_mahalanobis() requires at least two training observations.", call. = FALSE)
  # Preserve small represented coordinates beside a distant first row. This
  # origin is retained for prediction; the centered SVD arithmetic stays fixed.
  origin <- .numeric_signal_level(X)
  Z <- matrix(0, nrow(X), ncol(X))
  column_power <- numeric(ncol(X))
  spread_power <- rep(-Inf, ncol(X))
  # First center safely on each column's scale. Choosing a global scale from
  # raw levels would let a huge constant erase a tiny but varying coordinate.
  for (j in seq_len(ncol(X))) {
    amplitude <- max(abs(X[, j]))
    power <- if (amplitude == 0) 0 else
      max(-1074, min(1022, floor(log(amplitude) / log(2))))
    z <- X[, j] / 2^power
    z <- z - z[1L]
    z <- z - mean(z)
    Z[, j] <- z
    column_power[j] <- power
    spread <- max(abs(z))
    if (spread > 0) spread_power[j] <- power + floor(log(spread) / log(2))
  }
  active <- which(is.finite(spread_power))
  power <- if (!length(active)) 0 else max(-1074, min(1022, max(spread_power)))
  scale <- 2^power
  for (j in active) Z[, j] <- Z[, j] * 2^(column_power[j] - power)
  sv <- svd(Z, nu = 0L, nv = min(dim(Z)))
  keep <- which(sv$d > max(dim(X)) * .Machine$double.eps * max(sv$d))
  factor <- matrix(0, ncol(X), ncol(X))
  if (length(keep)) factor[, seq_along(keep)] <-
    sweep(sv$v[, keep, drop = FALSE], 2L, sv$d[keep] / sqrt(nrow(X) - 1L), "/")
  factor[!is.finite(spread_power), ] <- 0
  list(version = 2L, scale = scale, origin = origin, factor = factor,
    rank = as.integer(length(keep)))
}

.mahal_project <- function(X, ctx) {
  if (ctx$rank == 0L) return(matrix(0, nrow(X), ncol(ctx$factor)))
  active <- which(rowSums(abs(ctx$factor)) > 0)
  Z <- X[, active, drop = FALSE]
  exceptional <- rep(FALSE, nrow(X))
  for (j in seq_along(active)) {
    origin <- ctx$origin[active[j]]
    difference <- Z[, j] - origin
    recover <- !is.finite(difference) & ctx$scale > 1
    nonzero <- difference != 0
    difference <- difference / ctx$scale
    difference[recover] <- Z[recover, j] / ctx$scale - origin / ctx$scale
    exceptional <- exceptional | !is.finite(difference) | abs(difference) >= 2^1023 |
      (nonzero & abs(difference) < .Machine$double.xmin)
    Z[, j] <- difference
  }
  Z[exceptional, ] <- 0
  B <- ctx$factor[active, , drop = FALSE]
  Z <- .compensated_matprod(Z, B)
  # Normalize exceptional rows before projection, then restore the binary
  # ratio afterward. A factor below one can make a final projection finite
  # even when division by the training scale alone would overflow.
  for (i in which(exceptional)) {
    difference <- X[i, active] - ctx$origin[active]
    finite_difference <- all(is.finite(difference))
    amplitude <- if (finite_difference) max(abs(difference)) else
      max(abs(c(X[i, active], ctx$origin[active])))
    power <- if (amplitude == 0) 0 else
      max(-1074, min(1022, floor(log(amplitude) / log(2))))
    normalized <- if (finite_difference) difference / 2^power else
      X[i, active] / 2^power - ctx$origin[active] / 2^power
    value <- as.numeric(.compensated_matprod(matrix(normalized, nrow = 1L), B))
    remaining <- power - round(log(ctx$scale) / log(2))
    while (remaining != 0) {
      step <- max(-1074, min(1022, remaining))
      value <- value * 2^step
      remaining <- remaining - step
    }
    Z[i, ] <- value
  }
  if (any(!is.finite(Z)))
    stop("Mahalanobis projection is outside the numeric range; rescale covariates and targets.", call. = FALSE)
  Z
}

.mahal_eval_batch <- function(X, X0, ctx) {
  if (is.null(ctx)) ctx <- .mahal_fit(X)
  p <- ncol(X)
  if (!is.list(ctx) || !is.null(ctx$Si) ||
      !is.numeric(ctx$version) || length(ctx$version) != 1L ||
      !is.finite(ctx$version) || ctx$version != 2 ||
      !is.numeric(ctx$scale) || length(ctx$scale) != 1L ||
      !is.finite(ctx$scale) || ctx$scale <= 0 ||
      !is.numeric(ctx$origin) || length(ctx$origin) != p || any(!is.finite(ctx$origin)) ||
      !is.matrix(ctx$factor) || !is.numeric(ctx$factor) ||
      !identical(dim(ctx$factor), c(p, p)) || any(!is.finite(ctx$factor)) ||
      !is.numeric(ctx$rank) || length(ctx$rank) != 1L || !is.finite(ctx$rank) ||
      ctx$rank < 0 || ctx$rank > p || ctx$rank != floor(ctx$rank)) {
    stop(
      "`ctx` must be a current fitted Mahalanobis whitening context; refit with distance_mahalanobis()$setup(X).",
      call. = FALSE
    )
  }
  distance_l2()$eval_batch(.mahal_project(X, ctx), .mahal_project(X0, ctx))
}

# ── Public spec constructor ───────────────────────────────────────────────────

#' Mahalanobis distance spec (with one-time factorization)
#'
#' Returns a distance spec that fits the sample covariance geometry from the
#' training matrix X in `setup(X)`, then evaluates
#' \deqn{d(X_i, x_0) = \sqrt{(X_i - x_0)^\top \hat\Sigma^{-1}(X_i - x_0)}}
#' in `eval()` / `eval_batch()`. Requires at least two observations and one
#' covariate. The inverse in this formula is the Moore-Penrose inverse when
#' the centered design is numerically rank deficient.
#'
#' A whitening factor is computed once by SVD of the binary-normalized,
#' centered training design, avoiding an explicit covariance inverse and
#' squaring the design's condition number. Singular values greater than
#' `max(nrow(X), ncol(X)) * .Machine$double.eps * max(singular_values)`
#' are retained. Null directions contribute zero distance; an all-constant
#' training design gives zero distances. The cutoff uses variation in the
#' original covariate units; changing individual column units can therefore
#' change the retained numerical rank. The factor is reused through `ctx`,
#' including for subsets and new targets, following [distance_l2_scaled()].
#' Euclidean norms after projection use the selected C++ or R distance engine.
#' The explicit rank policy can change neighbors relative to older fits that
#' truncated covariance singular values after an inverse failed. Different
#' linear-algebra libraries can still differ at a numerical rank boundary.
#'
#' Particularly useful for DGPs with correlated covariates (e.g.\
#' `"dgp5_correlated"`) where L2 or diagonal-scaled L2 distances fail to
#' adapt to the off-diagonal covariance structure.
#'
#' @return A distance spec (list with `name`, `setup`, `eval`, `eval_batch`)
#'   for use with [cate_dnn()], [weights_dnn()], or [register_distance()].
#'   The string alias `"mahalanobis"` is registered at package load.
#'
#' @seealso [distance_l2()], [distance_l2_scaled()], [cate_dnn()]
#'
#' @examples
#' set.seed(1)
#' X  <- matrix(rnorm(100), 50, 2)
#' x0 <- c(0, 0)
#' dspec <- distance_mahalanobis()
#' ctx   <- dspec$setup(X)
#' head(dspec$eval(X, x0, ctx))
#'
#' @export
distance_mahalanobis <- function() {
  list(
    name  = "mahalanobis",
    setup = function(X) {
      .assert_no_na(X)
      X <- as.matrix(X)
      .mahal_fit(X)
    },
    eval = function(X, x0, ctx = NULL) {
      .assert_no_na(X, x0)
      X  <- as.matrix(X)
      x0 <- .validate_target_vector(x0, ncol(X), x0_arg = "x0", x_arg = "X")
      as.numeric(.mahal_eval_batch(X, matrix(x0, nrow = 1L), ctx))
    },
    eval_batch = function(X, X0, ctx = NULL) {
      .assert_no_na(X)
      .assert_no_na(X0)
      X  <- as.matrix(X)
      X0 <- .validate_ncol_relation(X, X0, y_arg = "X0", x_arg = "X", relation = "match")
      .mahal_eval_batch(X, X0, ctx)
    }
  )
}
