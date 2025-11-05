# R/dist_l2.R
# -----------------------------------------------------------------------------
# L2 (Euclidean) distances: fast C++ path with robust R fallbacks.

#' @noRd
.dist_l2_R <- function(X, x0, squared = FALSE) {
  X  <- as.matrix(X); x0 <- as.numeric(x0)
  if (length(x0) != ncol(X)) stop("x0 length must equal ncol(X).")
  Z <- X - matrix(x0, nrow(X), ncol(X), byrow = TRUE)
  s2 <- rowSums(Z * Z)
  if (isTRUE(squared)) s2 else sqrt(s2)
}

# -----------------------------------------------------------------------------
#' Fast L2 (Euclidean) distance to a single query point
#'
#' Computes \eqn{d_i = \|X_{i,\cdot} - x_0\|_2} for each row.
#'
#' @param X Numeric matrix \eqn{n \times p}.
#' @param x0 Numeric vector of length \eqn{p}.
#' @param squared Logical; if `TRUE`, return squared distances.
#' @return Numeric vector of length \eqn{n}.
#' @export
dist_l2 <- function(X, x0, squared = FALSE) {
  .assert_no_na(X, x0)
  X  <- as.matrix(X); x0 <- as.numeric(x0)
  if (length(x0) != ncol(X)) stop("x0 length must equal ncol(X).")

  cpp <- get0("rustat_dist_l2", mode = "function")
  if (!is.null(cpp)) return(cpp(X, x0, isTRUE(squared)))
  .dist_l2_R(X, x0, squared)
}

# -----------------------------------------------------------------------------
#' Order rows by L2 distance to a point (stable ties)
#'
#' @inheritParams dist_l2
#' @return A list with `dx` and `ord`.
#' @export
order_by_l2 <- function(X, x0, squared = FALSE) {
  .assert_no_na(X, x0)
  X  <- as.matrix(X); x0 <- as.numeric(x0)

  cpp <- get0("rustat_order_by_l2", mode = "function")
  if (!is.null(cpp)) return(cpp(X, x0, isTRUE(squared)))

  dx  <- dist_l2(X, x0, squared)
  ord <- order(dx, seq_along(dx))
  list(dx = dx, ord = as.integer(ord))
}

# -----------------------------------------------------------------------------
#' L2 distance spec (single and batched evaluation)
#'
#' Returns a spec with `eval()` and `eval_batch()` closures for L2 distances.
#' Uses compiled paths when present.
#'
#' @param squared Logical; if `TRUE`, return squared distances.
#' @return A list with:
#' * `name`: `"l2"` or `"l2_sq"`
#' * `eval(X, x0, ctx=NULL)`
#' * `eval_batch(X, X0, ctx=NULL)`
#'
#' @section Batch evaluation:
#' `eval_batch(X, X0)` expects an \eqn{m \times p} matrix of queries `X0`.
#'
#' @export
distance_l2 <- function(squared = FALSE) {
  list(
    name = if (isTRUE(squared)) "l2_sq" else "l2",
    eval = function(X, x0, ctx = NULL) dist_l2(X, x0, squared),
    eval_batch = function(X, X0, ctx = NULL) {
      X  <- as.matrix(X); X0 <- as.matrix(X0)
      .assert_no_na(X); .assert_no_na(X0)
      if (ncol(X0) != ncol(X)) stop("ncol(X0) must match ncol(X).")

      cpp <- get0("rustat_dist_l2_batch", mode = "function")
      if (!is.null(cpp)) return(cpp(X, X0, isTRUE(squared)))

      m <- nrow(X0)
      out <- matrix(NA_real_, nrow(X), m)
      for (j in seq_len(m)) out[, j] <- dist_l2(X, X0[j, ], squared)
      out
    }
  )
}

# -----------------------------------------------------------------------------
#' Scaled L2 distance (per-coordinate scaling)
#'
#' Computes \eqn{\|(X_{i,\cdot} - x_0) \odot \mathrm{inv\_s}\|_2} (or squared).
#'
#' @inheritParams dist_l2
#' @param inv_s Numeric vector of length \eqn{p} with elementwise scales.
#' @return Numeric vector of length \eqn{n}.
#' @export
dist_l2_scaled <- function(X, x0, inv_s, squared = FALSE) {
  .assert_no_na(X, x0)
  X  <- as.matrix(X); x0 <- as.numeric(x0); inv_s <- as.numeric(inv_s)
  if (length(x0) != ncol(X) || length(inv_s) != ncol(X)) {
    stop("x0 and inv_s must be length p == ncol(X).")
  }

  cpp <- get0("rustat_dist_l2_scaled", mode = "function")
  if (!is.null(cpp)) return(cpp(X, x0, inv_s, isTRUE(squared)))

  Xs  <- .scale_cols(X, inv_s)
  x0s <- x0 * inv_s
  dist_l2(Xs, x0s, squared = squared)
}

# -----------------------------------------------------------------------------
#' Scaled L2 distance spec (one-time setup + eval/batch)
#'
#' Builder that computes per-column scales once (`setup(X)`) and reuses them in
#' `eval()` / `eval_batch()`. Uses `rustat_dist_l2_scaled_batch()` if present.
#'
#' @param rule Scaling rule: `"sd"`, `"mad"`, `"iqr"`, `"range"`, or `"custom"`.
#' @param scale_vec Optional numeric vector of length \eqn{p} for
#'   `rule = "custom"`.
#' @param eps Small positive number to avoid division by zero in scaling.
#' @param squared Logical; if `TRUE`, return squared distances.
#'
#' @return A list with `name`, `setup()`, `eval()`, and `eval_batch()`.
#' @export
distance_l2_scaled <- function(rule = c("sd", "mad", "iqr", "range", "custom"),
                               scale_vec = NULL, eps = 1e-12,
                               squared = FALSE) {
  rule <- match.arg(rule)
  base <- distance_l2(squared = squared)

  list(
    name = paste0(if (isTRUE(squared)) "l2_scaled_sq[" else "l2_scaled[", rule, "]"),
    setup = function(X) {
      .assert_no_na(X)
      inv_s <- .compute_inv_scale(X, rule, scale_vec, eps)
      list(inv_s = inv_s, Xs = .scale_cols(as.matrix(X), inv_s))
    },
    eval = function(X, x0, ctx) {
      .assert_no_na(X, x0)
      inv_s <- if (!is.null(ctx) && !is.null(ctx$inv_s)) ctx$inv_s else
        .compute_inv_scale(X, rule, scale_vec, eps)

      cpp <- get0("rustat_dist_l2_scaled", mode = "function")
      if (!is.null(cpp)) return(cpp(as.matrix(X), as.numeric(x0), as.numeric(inv_s), isTRUE(squared)))

      Xs  <- if (!is.null(ctx) && !is.null(ctx$Xs)) ctx$Xs else .scale_cols(as.matrix(X), inv_s)
      x0s <- as.numeric(x0) * as.numeric(inv_s)
      base$eval(Xs, x0s, NULL)
    },
    eval_batch = function(X, X0, ctx) {
      .assert_no_na(X); .assert_no_na(X0)
      if (is.null(ctx) || is.null(ctx$inv_s)) {
        inv_s <- .compute_inv_scale(X, rule, scale_vec, eps)
        ctx <- list(inv_s = inv_s)
      }

      cpp <- get0("rustat_dist_l2_scaled_batch", mode = "function")
      if (!is.null(cpp)) return(cpp(as.matrix(X), as.matrix(X0), as.numeric(ctx$inv_s), isTRUE(squared)))

      Xs  <- if (!is.null(ctx$Xs)) ctx$Xs else .scale_cols(as.matrix(X), ctx$inv_s)
      X0s <- .scale_cols(as.matrix(X0), ctx$inv_s)
      base$eval_batch(Xs, X0s, NULL)
    }
  )
}

# -----------------------------------------------------------------------------
#' Order rows by scaled L2 (stable ties)
#'
#' Computes per-column scales (or uses `scale_vec` for `"custom"`), then prefers
#' compiled `rustat_order_by_l2_scaled()` if available. Otherwise computes
#' distances via [dist_l2_scaled()] and orders deterministically.
#'
#' @param X Numeric matrix \eqn{n \times p}.
#' @param x0 Numeric vector of length \eqn{p}.
#' @inheritParams distance_l2_scaled
#'
#' @return A list with `dx` and `ord`.
#' @export
order_by_l2_scaled <- function(X, x0,
                               rule = c("sd", "mad", "iqr", "range", "custom"),
                               scale_vec = NULL, eps = 1e-12,
                               squared = FALSE) {
  .assert_no_na(X, x0)
  X  <- as.matrix(X); x0 <- as.numeric(x0)
  rule <- match.arg(rule)

  inv_s <- if (rule == "custom") {
    .coerce_inv_s(scale_vec, ncol(X), eps)
  } else {
    .compute_inv_scale(X, rule, scale_vec, eps)
  }

  cpp <- get0("rustat_order_by_l2_scaled", mode = "function")
  if (!is.null(cpp)) return(cpp(X, x0, as.numeric(inv_s), isTRUE(squared)))

  dx  <- dist_l2_scaled(X, x0, inv_s, squared)
  ord <- order(dx, seq_along(dx))
  list(dx = dx, ord = as.integer(ord))
}

# --- internal helper ----------------------------------------------------------
#' @noRd
.coerce_inv_s <- function(scale_vec, p, eps) {
  inv_s <- as.numeric(scale_vec)
  if (length(inv_s) != p) stop("scale_vec must be length p when rule = 'custom'.")
  inv_s <- 1 / pmax(abs(inv_s), eps)
  inv_s
}
