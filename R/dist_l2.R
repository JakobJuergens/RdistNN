# L2 (Euclidean) distances: fast C++ path with robust R fallbacks.

# Stata reserves 2^1023 and larger doubles for missing values. Use the same
# intermediate boundary in every engine so exceptional arithmetic keeps ranks.
.l2_range_limit <- 2^1023

# Exceptional rows use a scaled sum of squares with scalar workspace. Recompute
# differences only here so ordinary values retain subtract-then-scale rounding.
.l2_range_row <- function(X, i, x0, inv_s, squared) {
  scale <- 0
  ssq <- 1
  for (j in seq_len(ncol(X))) {
    d <- X[i, j] - x0[j]
    if (!is.null(inv_s)) d <- .scaled_distance_difference(X[i, j], x0[j], inv_s[j])
    a <- abs(d)
    if (!is.finite(a)) return(a)
    if (a != 0) {
      if (scale < a) {
        ratio <- scale / a
        term <- ssq * (ratio * ratio)
        ssq <- 1 + term
        scale <- a
      } else {
        ratio <- a / scale
        ssq <- ssq + ratio * ratio
      }
    }
  }
  # Multiplying the normalized sum first recovers subnormal squared results
  # even when every individual square rounds to zero.
  if (isTRUE(squared)) (scale * ssq) * scale else scale * sqrt(ssq)
}

#' @noRd
.dist_l2_R <- function(X, x0, squared = FALSE, inv_s = NULL) {
  X <- as.matrix(X)
  x0 <- .validate_target_vector(x0, ncol(X), x0_arg = "x0", x_arg = "X", finite = TRUE)
  # Match native arithmetic, including double-precision accumulation order;
  # pre-scaling inputs or rowSums() can change distances and stable ranks.
  s2 <- numeric(nrow(X))
  exceptional <- FALSE
  for (j in seq_len(ncol(X))) {
    d <- X[, j] - x0[j]
    if (!is.null(inv_s)) {
      if (any(abs(d) >= .l2_range_limit)) {
        exceptional <- exceptional | (abs(d) >= .l2_range_limit)
      }
      d <- d * inv_s[j]
    }
    square <- d * d
    # Only allocate a row mask if a square can have underflowed. Nonfinite
    # squares propagate to s2 and are detected once after accumulation.
    if (any(square < .Machine$double.xmin, na.rm = TRUE)) {
      exceptional <- exceptional | (d != 0 & square < .Machine$double.xmin)
    }
    s2 <- s2 + square
  }
  exceptional <- exceptional | !is.finite(s2) | (s2 >= .l2_range_limit)
  out <- if (isTRUE(squared)) s2 else sqrt(s2)
  for (i in which(exceptional)) out[i] <- .l2_range_row(X, i, x0, inv_s, squared)
  out
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
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(12), 4, 3)
#' dist_l2(X, c(0, 0, 0))
#' @keywords internal
#' @noRd
dist_l2 <- function(X, x0, squared = FALSE) {
  .assert_no_na(X, x0)
  X <- as.matrix(X)
  x0 <- .validate_target_vector(x0, ncol(X), x0_arg = "x0", x_arg = "X", finite = TRUE)

  cpp <- .rdistnn_compiled_fun("rustat_dist_l2")
  if (!is.null(cpp)) {
    return(cpp(X, x0, isTRUE(squared)))
  }
  .dist_l2_R(X, x0, squared)
}

# -----------------------------------------------------------------------------
#' Order rows by L2 distance to a point (stable ties)
#'
#' @inheritParams dist_l2
#' @return A list with `dx` and `ord`.
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(10), 5, 2)
#' order_by_l2(X, c(0, 0))
#' @keywords internal
#' @noRd
order_by_l2 <- function(X, x0, squared = FALSE) {
  .assert_no_na(X, x0)
  X <- as.matrix(X)
  x0 <- .validate_target_vector(x0, ncol(X), x0_arg = "x0", x_arg = "X", finite = TRUE)

  cpp <- .rdistnn_compiled_fun("rustat_order_by_l2")
  if (!is.null(cpp)) {
    return(cpp(X, x0, isTRUE(squared)))
  }

  dx <- dist_l2(X, x0, squared)
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
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(12), 4, 3)
#' X0 <- matrix(rnorm(6), 2, 3)
#' dspec <- distance_l2()
#' dspec$eval_batch(X, X0)
#' @export
distance_l2 <- function(squared = FALSE) {
  list(
    name  = if (isTRUE(squared)) "l2_sq" else "l2",
    setup = function(X) list(),
    eval  = function(X, x0, ctx = NULL) dist_l2(X, x0, squared),
    eval_batch = function(X, X0, ctx = NULL) {
      X <- as.matrix(X)
      X0 <- as.matrix(X0)
      .assert_no_na(X)
      .assert_no_na(X0)
      X0 <- .validate_ncol_relation(X, X0, y_arg = "X0", x_arg = "X", relation = "match")

      cpp <- .rdistnn_compiled_fun("rustat_dist_l2_batch")
      if (!is.null(cpp)) {
        return(cpp(X, X0, isTRUE(squared)))
      }

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
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(12), 4, 3)
#' dist_l2_scaled(X, c(0, 0, 0), inv_s = c(1, 2, 0.5))
#' @keywords internal
#' @noRd
dist_l2_scaled <- function(X, x0, inv_s, squared = FALSE) {
  .assert_no_na(X, x0)
  X <- as.matrix(X)
  validated <- .validate_x0_inv_s(x0, inv_s, ncol(X), x_arg = "X", strict = TRUE, finite = TRUE)
  x0 <- validated$x0
  inv_s <- validated$inv_s

  cpp <- .rdistnn_compiled_fun("rustat_dist_l2_scaled")
  if (!is.null(cpp)) {
    return(cpp(X, x0, inv_s, isTRUE(squared)))
  }

  .dist_l2_R(X, x0, squared, inv_s)
}

# -----------------------------------------------------------------------------
#' Scaled L2 distance spec (one-time setup + eval/batch)
#'
#' Builder that computes per-column scales once (`setup(X)`) and reuses them in
#' `eval()` / `eval_batch()`. Uses a compiled fast path when available.
#'
#' @param rule Scaling rule: `"sd"`, `"mad"`, `"iqr"`, `"range"`, or `"custom"`.
#' @param scale_vec Optional numeric vector of length \eqn{p} for
#'   `rule = "custom"`.
#' @param eps Small positive number to avoid division by zero in scaling.
#' @param squared Logical; if `TRUE`, return squared distances.
#'
#' @return A list with `name`, `setup()`, `eval()`, and `eval_batch()`.
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(12), 4, 3)
#' X0 <- matrix(rnorm(6), 2, 3)
#' dspec <- distance_l2_scaled(rule = "sd")
#' ctx <- dspec$setup(X)
#' dspec$eval_batch(X, X0, ctx)
#' @export
distance_l2_scaled <- function(rule = c("sd", "mad", "iqr", "range", "custom"),
                               scale_vec = NULL, eps = 1e-12,
                               squared = FALSE) {
  rule <- match.arg(rule)

  list(
    name = paste0(if (isTRUE(squared)) "l2_scaled_sq[" else "l2_scaled[", rule, "]"),
    setup = function(X) {
      .assert_no_na(X)
      inv_s <- .compute_inv_scale(X, rule, scale_vec, eps)
      list(inv_s = inv_s, Xs = .scale_cols(as.matrix(X), inv_s))
    },
    eval = function(X, x0, ctx) {
      .assert_no_na(X, x0)
      X <- as.matrix(X)
      inv_s <- if (!is.null(ctx) && !is.null(ctx$inv_s)) {
        ctx$inv_s
      } else {
        .compute_inv_scale(X, rule, scale_vec, eps)
      }
      validated <- .validate_x0_inv_s(x0, inv_s, ncol(X), x_arg = "X", strict = TRUE, finite = TRUE)
      x0 <- validated$x0
      inv_s <- validated$inv_s

      cpp <- .rdistnn_compiled_fun("rustat_dist_l2_scaled")
      if (!is.null(cpp)) {
        return(cpp(X, x0, inv_s, isTRUE(squared)))
      }

      .dist_l2_R(X, x0, squared, inv_s)
    },
    eval_batch = function(X, X0, ctx) {
      .assert_no_na(X)
      .assert_no_na(X0)
      X <- as.matrix(X)
      X0 <- .validate_ncol_relation(X, X0, y_arg = "X0", x_arg = "X", relation = "match")
      if (is.null(ctx) || is.null(ctx$inv_s)) {
        inv_s <- .compute_inv_scale(X, rule, scale_vec, eps)
        ctx <- list(inv_s = inv_s)
      }
      if (!is.numeric(ctx$inv_s)) {
        stop("inv_s must be numeric.", call. = FALSE)
      }
      inv_s <- as.numeric(ctx$inv_s)
      if (.has_na_only(inv_s)) {
        stop("inv_s cannot contain NA values.", call. = FALSE)
      }
      if (.has_nonfinite_non_na(inv_s)) {
        stop("inv_s cannot contain non-finite values (`Inf`, `-Inf`, or `NaN`).", call. = FALSE)
      }
      if (length(inv_s) != ncol(X)) {
        stop("inv_s length must equal ncol(X).", call. = FALSE)
      }

      cpp <- .rdistnn_compiled_fun("rustat_dist_l2_scaled_batch")
      if (!is.null(cpp)) {
        return(cpp(X, X0, inv_s, isTRUE(squared)))
      }

      out <- matrix(NA_real_, nrow(X), nrow(X0))
      for (j in seq_len(nrow(X0))) out[, j] <- .dist_l2_R(X, X0[j, ], squared, inv_s)
      out
    }
  )
}

# -----------------------------------------------------------------------------
#' Order rows by scaled L2 (stable ties)
#'
#' Computes per-column scales (or uses `scale_vec` for `"custom"`), then
#' prefers a compiled fast path if available. Otherwise computes distances
#' via `dist_l2_scaled()` and orders deterministically.
#'
#' @param X Numeric matrix \eqn{n \times p}.
#' @param x0 Numeric vector of length \eqn{p}.
#' @inheritParams distance_l2_scaled
#'
#' @return A list with `dx` and `ord`.
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(10), 5, 2)
#' order_by_l2_scaled(X, c(0, 0), rule = "sd")
#' @keywords internal
#' @noRd
order_by_l2_scaled <- function(X, x0,
                               rule = c("sd", "mad", "iqr", "range", "custom"),
                               scale_vec = NULL, eps = 1e-12,
                               squared = FALSE) {
  .assert_no_na(X, x0)
  X <- as.matrix(X)
  x0 <- .validate_target_vector(x0, ncol(X), x0_arg = "x0", x_arg = "X", finite = TRUE)
  rule <- match.arg(rule)

  inv_s <- if (rule == "custom") {
    .coerce_inv_s(scale_vec, ncol(X), eps)
  } else {
    .compute_inv_scale(X, rule, scale_vec, eps)
  }

  cpp <- .rdistnn_compiled_fun("rustat_order_by_l2_scaled")
  if (!is.null(cpp)) {
    return(cpp(X, x0, as.numeric(inv_s), isTRUE(squared)))
  }

  dx <- dist_l2_scaled(X, x0, inv_s, squared)
  ord <- order(dx, seq_along(dx))
  list(dx = dx, ord = as.integer(ord))
}
