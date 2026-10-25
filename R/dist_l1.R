# L1 (Manhattan) distances: fast C++ path with robust R fallbacks.

# --- internal R fallback ------------------------------------------------------
#' @noRd
.dist_l1_R <- function(X, x0, inv_s = NULL) {
  X <- as.matrix(X)
  x0 <- .validate_target_vector(x0, ncol(X), x0_arg = "x0", x_arg = "X", finite = TRUE)
  # Match native arithmetic: subtract, scale, then accumulate in coordinate
  # order with double precision. rowSums() can break native ties differently.
  out <- numeric(nrow(X))
  exceptional <- FALSE
  for (j in seq_len(ncol(X))) {
    d <- X[, j] - x0[j]
    if (!is.null(inv_s)) {
      raw <- d
      d <- if (any(abs(d) >= 2^1023)) {
        .scaled_distance_difference(X[, j], x0[j], inv_s[j])
      } else {
        d * inv_s[j]
      }
      exceptional <- exceptional | abs(raw) >= 2^1023 |
        (raw != 0 & inv_s[j] != 0 & abs(d) < .Machine$double.xmin)
    }
    out <- out + abs(d)
  }
  for (i in which(exceptional)) {
    out[i] <- .scaled_distance_norm_row(X, i, x0, inv_s, l2 = FALSE)
  }
  out
}

# -----------------------------------------------------------------------------
#' Fast L1 (Manhattan) distance to a single query point
#'
#' Computes \eqn{d_i = \|X_{i,\cdot} - x_0\|_1} for each row \eqn{i=1,\dots,n}.
#' A compiled fast path is used when available; otherwise a vectorised R
#' fallback is used.
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
#' @seealso `dist_l2()`, [distance_l1_scaled()], [distance_l1()], [weights_matrix()]
#' @keywords internal
#' @noRd
dist_l1 <- function(X, x0) {
  .assert_no_na(X, x0)
  X <- as.matrix(X)
  x0 <- .validate_target_vector(x0, ncol(X), x0_arg = "x0", x_arg = "X", finite = TRUE)

  cpp <- .rdistnn_compiled_fun("rustat_dist_l1")
  if (!is.null(cpp)) {
    return(cpp(X, x0))
  }
  .dist_l1_R(X, x0)
}

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
#' L1 distance spec (single and batched evaluation)
#'
#' Returns a spec with `eval()` and `eval_batch()` closures for L1 distances.
#' Compiled fast paths are used when available; otherwise R fallbacks are
#' used.
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
#' @seealso [distance_l1_scaled()], [distance_eval()], [weights_matrix()]
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
      X0 <- .validate_ncol_relation(X, X0, y_arg = "X0", x_arg = "X", relation = "match")

      cpp <- .rdistnn_compiled_fun("rustat_dist_l1_batch")
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

# -----------------------------------------------------------------------------
#' Scaled L1 distance spec (one-time setup + eval/batch)
#'
#' Builder that computes per-column scales once (`setup(X)`) and reuses them in
#' `eval()` / `eval_batch()`. Uses a compiled fast path when available.
#'
#' @param rule Scaling rule: `"sd"`, `"mad"`, `"iqr"`, `"range"`, or `"custom"`.
#' @param scale_vec Optional numeric vector of length \eqn{p} for
#'   `rule = "custom"`.
#' @param eps Small positive number to avoid division by zero in scaling.
#'
#' @details Configuration is evaluated when the specification is created.
#' Rebinding a caller's configuration variables before first use does not alter
#' it. Setup retains per-column inverse scales and does not retain a scaled copy
#' of the training matrix.
#' @return A list with `name`, `setup()`, `eval()`, and `eval_batch()`.
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(12), 4, 3)
#' X0 <- matrix(rnorm(6), 2, 3)
#' dspec <- distance_l1_scaled(rule = "sd")
#' ctx <- dspec$setup(X)
#' dspec$eval_batch(X, X0, ctx)
#' @export
distance_l1_scaled <- function(rule = c("sd", "mad", "iqr", "range", "custom"),
                               scale_vec = NULL, eps = 1e-12) {
  rule <- .dnn_match_arg(rule)
  force(scale_vec)
  force(eps)

  list(
    name = paste0("l1_scaled[", rule, "]"),
    setup = function(X) {
      .assert_no_na(X)
      inv_s <- .compute_inv_scale(X, rule, scale_vec, eps)
      list(inv_s = inv_s)
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

      cpp <- .rdistnn_compiled_fun("rustat_dist_l1_scaled")
      if (!is.null(cpp)) {
        return(cpp(X, x0, inv_s))
      }

      .dist_l1_R(X, x0, inv_s)
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
        .dnn_abort("argument_invalid", "`inv_s` must be numeric.", argument = "inv_s")
      }
      inv_s <- as.numeric(ctx$inv_s)
      if (.has_na_only(inv_s)) {
        .dnn_abort("argument_invalid", "`inv_s` cannot contain NA values.", argument = "inv_s")
      }
      if (.has_nonfinite_non_na(inv_s)) {
        .dnn_abort("argument_invalid",
          "`inv_s` cannot contain non-finite values (`Inf`, `-Inf`, or `NaN`).",
          argument = "inv_s"
        )
      }
      if (length(inv_s) != ncol(X)) {
        .dnn_abort("argument_invalid", "`inv_s` length must equal ncol(X).", argument = "inv_s")
      }

      cpp <- .rdistnn_compiled_fun("rustat_dist_l1_scaled_batch")
      if (!is.null(cpp)) {
        return(cpp(X, X0, inv_s))
      }

      out <- matrix(NA_real_, nrow(X), nrow(X0))
      for (j in seq_len(nrow(X0))) out[, j] <- .dist_l1_R(X, X0[j, ], inv_s)
      out
    }
  )
}

# -----------------------------------------------------------------------------
