# dist_mahal.R
# Mahalanobis distance spec with one-time covariance estimation.
#
# Follows the same setup-then-eval convention as distance_l2_scaled():
#   setup(X) computes Si = cov(X)^{-1} once;
#   eval / eval_batch use Si from ctx.
#
# The string alias "mahalanobis" is registered at package load via
# .register_builtin_distances() in api-distances.R.

# ── Internal helpers ──────────────────────────────────────────────────────────

# Safe pseudoinverse via SVD; avoids dependency on MASS.
.mahal_solve_safe <- function(S, eps = 1e-10) {
  result <- tryCatch(solve(S), error = function(e) NULL)
  if (!is.null(result)) return(result)
  sv    <- svd(S)
  d_inv <- ifelse(sv$d > eps * max(sv$d, 1e-300), 1 / sv$d, 0)
  sv$v %*% diag(d_inv, nrow = length(d_inv)) %*% t(sv$u)
}

# Mahalanobis distances from a single target x0 given pre-computed Si.
.mahal_dist_one <- function(X, x0, Si) {
  Z  <- X - matrix(x0, nrow(X), ncol(X), byrow = TRUE)
  d2 <- rowSums((Z %*% Si) * Z)
  sqrt(pmax(d2, 0))
}

# ── Public spec constructor ───────────────────────────────────────────────────

#' Mahalanobis distance spec (with one-time covariance setup)
#'
#' Returns a distance spec that estimates the sample covariance from the
#' training matrix X in `setup(X)`, then evaluates
#' \deqn{d(X_i, x_0) = \sqrt{(X_i - x_0)^\top \hat\Sigma^{-1}(X_i - x_0)}}
#' in `eval()` / `eval_batch()`.  Requires at least two covariates.
#'
#' The covariance inverse is computed once in `setup()` and reused for all
#' evaluations via the `ctx` argument, following the same convention as
#' [distance_l2_scaled()].  When the sample covariance is (near-)singular a
#' pseudoinverse via SVD is used automatically; no external packages are
#' required.
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
      if (ncol(X) < 2L)
        stop("distance_mahalanobis() requires at least 2 covariates.", call. = FALSE)
      list(Si = .mahal_solve_safe(stats::cov(X)))
    },
    eval = function(X, x0, ctx = NULL) {
      .assert_no_na(X, x0)
      X  <- as.matrix(X)
      x0 <- .validate_target_vector(x0, ncol(X), x0_arg = "x0", x_arg = "X")
      Si <- if (!is.null(ctx) && !is.null(ctx$Si)) ctx$Si
            else .mahal_solve_safe(stats::cov(X))
      .mahal_dist_one(X, x0, Si)
    },
    eval_batch = function(X, X0, ctx = NULL) {
      .assert_no_na(X)
      .assert_no_na(X0)
      X  <- as.matrix(X)
      X0 <- .validate_ncol_relation(X, X0, y_arg = "X0", x_arg = "X", relation = "match")
      Si  <- if (!is.null(ctx) && !is.null(ctx$Si)) ctx$Si
             else .mahal_solve_safe(stats::cov(X))
      m   <- nrow(X0)
      out <- matrix(NA_real_, nrow(X), m)
      for (j in seq_len(m)) out[, j] <- .mahal_dist_one(X, X0[j, ], Si)
      out
    }
  )
}
