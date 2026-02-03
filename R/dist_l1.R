# R/dist_l1.R
# -----------------------------------------------------------------------------
# L1 (Manhattan) distances: fast C++ path with robust R fallbacks.

# --- internal R fallback ------------------------------------------------------
#' @noRd
.dist_l1_R <- function(X, x0) {
  X <- as.matrix(X)
  x0 <- as.numeric(x0)
  if (length(x0) != ncol(X)) stop("x0 length must equal ncol(X).")
  rowSums(abs(X - matrix(x0, nrow(X), ncol(X), byrow = TRUE)))
}

# -----------------------------------------------------------------------------
#' Fast L1 (Manhattan) distance to a single query point
#'
#' Computes \eqn{d_i = \|X_{i,\cdot} - x_0\|_1} for each row \eqn{i=1,\dots,n}.
#' If the compiled symbol `rustat_dist_l1()` is registered (via Rcpp), it is
#' used; otherwise a vectorised R fallback is used.
#'
#' @param X Numeric matrix \eqn{n \times p}. Each row is an observation.
#' @param x0 Numeric vector of length \eqn{p}. The query point.
#'
#' @return Numeric vector of length \eqn{n} with L1 distances.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(12), 4, 3)
#' x0 <- c(0, 1, -1)
#' dist_l1(X, x0)
#'
#' @seealso [dist_l2()], [dist_l1_scaled()], [distance_l1()], [order_by_l1()]
#' @export
dist_l1 <- function(X, x0) {
  .assert_no_na(X, x0)
  X <- as.matrix(X)
  x0 <- as.numeric(x0)
  if (length(x0) != ncol(X)) stop("x0 length must equal ncol(X).")

  cpp <- get0("rustat_dist_l1", mode = "function")
  if (!is.null(cpp)) {
    return(cpp(X, x0))
  }
  .dist_l1_R(X, x0)
}

# -----------------------------------------------------------------------------
#' Order rows by L1 distance to a point (stable ties)
#'
#' Returns distances and a stable ordering. If `rustat_order_by_l1()` exists,
#' it is used; otherwise `order(dx, seq_along(dx))` provides a deterministic
#' tie-break.
#'
#' @inheritParams dist_l1
#' @return A list with:
#' * `dx`: numeric vector (length \eqn{n})
#' * `ord`: integer vector of 1-based row indices
#'
#' @examples
#' set.seed(2)
#' X <- matrix(rnorm(10), 5, 2)
#' x0 <- c(0, 0)
#' order_by_l1(X, x0)
#'
#' @seealso [order_by_l1_scaled()], [dist_l1()]
#' @export
order_by_l1 <- function(X, x0) {
  .assert_no_na(X, x0)
  X <- as.matrix(X)
  x0 <- as.numeric(x0)

  cpp <- get0("rustat_order_by_l1", mode = "function")
  if (!is.null(cpp)) {
    return(cpp(X, x0))
  }

  dx <- dist_l1(X, x0)
  ord <- order(dx, seq_along(dx))
  list(dx = dx, ord = as.integer(ord))
}

# -----------------------------------------------------------------------------
#' L1 distance spec (single and batched evaluation)
#'
#' Returns a spec with `eval()` and `eval_batch()` closures for L1 distances.
#' The compiled paths (`rustat_dist_l1()`, `rustat_dist_l1_batch()`) are used
#' when available; otherwise R fallbacks are used.
#'
#' @return A list with:
#' * `name`: `"l1"`
#' * `eval(X, x0, ctx=NULL)`: numeric vector of length \eqn{n}
#' * `eval_batch(X, X0, ctx=NULL)`: numeric matrix \eqn{n \times m}
#'
#' @section Batch evaluation:
#' `eval_batch(X, X0)` expects an \eqn{m \times p} matrix of queries `X0`.
#'
#' @examples
#' set.seed(3)
#' X <- matrix(rnorm(12), 4, 3)
#' X0 <- matrix(rnorm(6), 2, 3)
#' dspec <- distance_l1()
#' dspec$eval(X, X0[1, ])
#' dspec$eval_batch(X, X0)
#'
#' @seealso [distance_l1_scaled()], [dist_l1()], [dist_l1_scaled()]
#' @export
distance_l1 <- function() {
  list(
    name = "l1",
    eval = function(X, x0, ctx = NULL) {
      dist_l1(X, x0)
    },
    eval_batch = function(X, X0, ctx = NULL) {
      X <- as.matrix(X)
      X0 <- as.matrix(X0)
      .assert_no_na(X)
      .assert_no_na(X0)
      if (ncol(X0) != ncol(X)) stop("ncol(X0) must match ncol(X).")

      cpp <- get0("rustat_dist_l1_batch", mode = "function")
      if (!is.null(cpp)) {
        return(cpp(X, X0))
      }

      m <- nrow(X0)
      out <- matrix(NA_real_, nrow(X), m)
      for (j in seq_len(m)) out[, j] <- dist_l1(X, X0[j, ])
      out
    }
  )
}

# -----------------------------------------------------------------------------
#' Scaled L1 distance (per-coordinate scaling)
#'
#' Computes \eqn{\|(X_{i,\cdot} - x_0) \odot \mathrm{inv\_s}\|_1} for each row.
#' If `rustat_dist_l1_scaled()` exists it is used; else inputs are rescaled
#' in R and passed to [dist_l1()].
#'
#' @inheritParams dist_l1
#' @param inv_s Numeric vector of length \eqn{p} with elementwise scales
#'   (often the reciprocals of column scale estimates).
#'
#' @return Numeric vector of length \eqn{n}.
#' @examples
#' X <- matrix(rnorm(12), 4, 3)
#' dist_l1_scaled(X, c(0, 0, 0), inv_s = c(1, 2, 0.5))
#' @export
dist_l1_scaled <- function(X, x0, inv_s) {
  .assert_no_na(X, x0)
  X <- as.matrix(X)
  x0 <- as.numeric(x0)
  inv_s <- as.numeric(inv_s)
  if (length(x0) != ncol(X) || length(inv_s) != ncol(X)) {
    stop("x0 and inv_s must be length p == ncol(X).")
  }

  cpp <- get0("rustat_dist_l1_scaled", mode = "function")
  if (!is.null(cpp)) {
    return(cpp(X, x0, inv_s))
  }

  Xs <- .scale_cols(X, inv_s)
  x0s <- x0 * inv_s
  dist_l1(Xs, x0s)
}

# -----------------------------------------------------------------------------
#' Scaled L1 distance spec (one-time setup + eval/batch)
#'
#' Builder that computes per-column scales once (`setup(X)`) and reuses them in
#' `eval()` / `eval_batch()`. Uses `rustat_dist_l1_scaled_batch()` if present.
#'
#' @param rule Scaling rule: `"sd"`, `"mad"`, `"iqr"`, `"range"`, or `"custom"`.
#' @param scale_vec Optional numeric vector of length \eqn{p} for
#'   `rule = "custom"`.
#' @param eps Small positive number to avoid division by zero in scaling.
#'
#' @return A list with `name`, `setup()`, `eval()`, and `eval_batch()`.
#' @examples
#' X <- matrix(rnorm(12), 4, 3)
#' X0 <- matrix(rnorm(6), 2, 3)
#' dspec <- distance_l1_scaled(rule = "sd")
#' ctx <- dspec$setup(X)
#' dspec$eval_batch(X, X0, ctx)
#' @export
distance_l1_scaled <- function(rule = c("sd", "mad", "iqr", "range", "custom"),
                               scale_vec = NULL, eps = 1e-12) {
  rule <- match.arg(rule)
  base <- distance_l1()

  list(
    name = paste0("l1_scaled[", rule, "]"),
    setup = function(X) {
      .assert_no_na(X)
      inv_s <- .compute_inv_scale(X, rule, scale_vec, eps)
      list(inv_s = inv_s, Xs = .scale_cols(as.matrix(X), inv_s))
    },
    eval = function(X, x0, ctx) {
      .assert_no_na(X, x0)
      inv_s <- if (!is.null(ctx) && !is.null(ctx$inv_s)) {
        ctx$inv_s
      } else {
        .compute_inv_scale(X, rule, scale_vec, eps)
      }

      cpp <- get0("rustat_dist_l1_scaled", mode = "function")
      if (!is.null(cpp)) {
        return(cpp(as.matrix(X), as.numeric(x0), as.numeric(inv_s)))
      }

      Xs <- if (!is.null(ctx) && !is.null(ctx$Xs)) ctx$Xs else .scale_cols(as.matrix(X), inv_s)
      x0s <- as.numeric(x0) * as.numeric(inv_s)
      base$eval(Xs, x0s, NULL)
    },
    eval_batch = function(X, X0, ctx) {
      .assert_no_na(X)
      .assert_no_na(X0)
      if (is.null(ctx) || is.null(ctx$inv_s)) {
        inv_s <- .compute_inv_scale(X, rule, scale_vec, eps)
        ctx <- list(inv_s = inv_s)
      }

      cpp <- get0("rustat_dist_l1_scaled_batch", mode = "function")
      if (!is.null(cpp)) {
        return(cpp(as.matrix(X), as.matrix(X0), as.numeric(ctx$inv_s)))
      }

      Xs <- if (!is.null(ctx$Xs)) ctx$Xs else .scale_cols(as.matrix(X), ctx$inv_s)
      X0s <- .scale_cols(as.matrix(X0), ctx$inv_s)
      base$eval_batch(Xs, X0s, NULL)
    }
  )
}

# -----------------------------------------------------------------------------
#' Order rows by scaled L1 (stable ties)
#'
#' Computes per-column scales (or uses `scale_vec` for `"custom"`), then prefers
#' compiled `rustat_order_by_l1_scaled()` if available. Otherwise computes
#' distances via [dist_l1_scaled()] and orders with a deterministic tie rule.
#'
#' @param X Numeric matrix \eqn{n \times p}.
#' @param x0 Numeric vector of length \eqn{p}.
#' @param rule Scaling rule: `"sd"`, `"mad"`, `"iqr"`, `"range"`, or `"custom"`.
#' @param scale_vec Optional numeric vector of length \eqn{p} for
#'   `rule = "custom"`.
#' @param eps Small positive number to avoid division by zero in scaling.
#'
#' @return A list with:
#' * `dx`: numeric vector (length \eqn{n})
#' * `ord`: integer vector of 1-based row indices
#'
#' @examples
#' set.seed(6)
#' X <- matrix(rnorm(10), 5, 2)
#' x0 <- c(0, 0)
#' order_by_l1_scaled(X, x0, rule = "sd")
#'
#' @seealso [order_by_l1()], [dist_l1_scaled()]
#' @export
order_by_l1_scaled <- function(X, x0,
                               rule = c("sd", "mad", "iqr", "range", "custom"),
                               scale_vec = NULL, eps = 1e-12) {
  .assert_no_na(X, x0)
  X <- as.matrix(X)
  x0 <- as.numeric(x0)
  rule <- match.arg(rule)

  inv_s <- if (rule == "custom") {
    .coerce_inv_s(scale_vec, ncol(X), eps)
  } else {
    .compute_inv_scale(X, rule, scale_vec, eps)
  }

  cpp <- get0("rustat_order_by_l1_scaled", mode = "function")
  if (!is.null(cpp)) {
    return(cpp(X, x0, as.numeric(inv_s)))
  }

  dx <- dist_l1_scaled(X, x0, inv_s)
  ord <- order(dx, seq_along(dx))
  list(dx = dx, ord = as.integer(ord))
}
