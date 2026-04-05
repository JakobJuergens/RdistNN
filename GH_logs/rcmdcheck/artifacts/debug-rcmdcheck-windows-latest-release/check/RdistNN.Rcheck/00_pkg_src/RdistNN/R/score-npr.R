#' Nonparametric regression **singleton** score from observed outcome \eqn{Y}
#' (order-1 pure score; no weights)
#'
#' Constructs a **pure-score** kernel for nonparametric regression that uses only
#' the observed response \eqn{Y}. For a singleton subset \eqn{S=\{i\}} the kernel
#' returns exactly \eqn{h(S) = Y_i}. This score is **order-1 only** (s = 1).
#' It performs no nearest-neighbor ordering, kernel weighting, or L-stat
#' aggregation. Any weighting/averaging of subset evaluations is handled outside
#' the score layer (e.g., by rank-/distance-weight modules).
#'
#' The arguments \code{x0}, \code{x_cols}, \code{s}, \code{k} are accepted for API
#' consistency and basic dimension checks; only \eqn{Y} is used to construct the
#' score. This factory **enforces** \code{s = 1}.
#'
#' @param x0 Numeric length-\eqn{p} target point. Used only to check that
#'   \code{length(x0)} matches \code{length(x_cols)}; it is not used in the kernel.
#' @param x_cols Integer or character indices of the \eqn{X} columns in \code{data},
#'   used for dimension checks with \code{x0}.
#' @param s Intended subset size (kernel order). **Must be 1**; otherwise an error
#'   is thrown.
#' @param k Integer included for interface compatibility with other scores; ignored.
#' @param y_col Optional name of the column in \code{data} that contains the response
#'   \eqn{Y}. If \code{NULL}, supply \code{y}.
#' @param y Optional numeric vector (length \eqn{n}) of \eqn{Y} values to use when
#'   \code{y_col} is \code{NULL}.
#'
#' @return An \code{ustat_score} object with:
#' \itemize{
#'   \item \code{$name}: a descriptive name,
#'   \item \code{$order}: always \code{1L},
#'   \item \code{$setup}: stores \code{g <- Y} in the context,
#'   \item \code{$eval}/\code{$eval_batch}: return \eqn{Y_i} for each singleton subset.
#' }
#'
#' @examples
#' set.seed(1)
#' n <- 50
#' p <- 2
#' dat <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
#' dat$Y <- with(dat, X1 - 2 * X2 + rnorm(n))
#'
#' sc <- score_np_y(
#'   x0 = c(0, 0),
#'   x_cols = c("X1", "X2"),
#'   s = 1L, k = 1L,
#'   y_col = "Y"
#' )
#' ctx <- score_setup(sc, dat)
#' subs <- matrix(seq_len(n), ncol = 1L) # singleton subsets
#' score_eval_batch(sc, subs, dat, ctx) # returns dat$Y
#'
#' @export
score_np_y <- function(
  x0, x_cols, s = 1, k = 1L,
  y_col = NULL, y = NULL
) {
  # XOR: exactly one of y_col or y
  if (is.null(y_col) == is.null(y)) {
    .abort_bad_state("Provide exactly one of 'y_col' or 'y' (mutually exclusive).")
  }
  # Enforce singleton score
  if (!is.numeric(s) || length(s) != 1L || s != 1L) {
    stop("`score_np_y` is a singleton (order-1) score: require `s = 1`.", call. = FALSE)
  }
  k <- .validate_scalar_int(k, arg = "k", min = 1L)

  new_score(
    name = "np_y_singleton_purescore",
    order = function(n, data, ctx) 1L,
    setup = function(data, ...) {
      X <- as.matrix(data[, x_cols, drop = FALSE])
      yy <- if (!is.null(y_col)) as.numeric(data[[y_col]]) else as.numeric(y)

      .validate_target_vector(x0, ncol(X), x0_arg = "x0", x_arg = "X")
      yy <- .validate_vector_nrow(yy, nrow(X), arg = "y", ref = "data")

      # Do NOT check finiteness here (fail at evaluation time per convention)
      list(g = yy, n = nrow(X), p = ncol(X), s = 1L)
    },

    # Runtime guards: fail at evaluation if non-finite anywhere (as in cate score)
    eval = function(idxs, data, ctx, rng) {
      if (any(!is.finite(ctx$g))) stop("Non-finite g encountered.", call. = FALSE)
      if (length(idxs) != 1L) stop("`score_np_y` expects a singleton subset in `eval()`.", call. = FALSE)
      as.numeric(ctx$g[as.integer(idxs[1L])])
    },
    eval_batch = function(subs, data, ctx, rng) {
      if (any(!is.finite(ctx$g))) stop("Non-finite g encountered.", call. = FALSE)
      if (!is.matrix(subs) || ncol(subs) != 1L) {
        stop("`score_np_y` expects `subs` to be an n x 1 matrix (singleton subsets).", call. = FALSE)
      }
      if (nrow(subs) == 0L) {
        return(numeric(0L))
      }

      if (exists("rustat_subset_row_sums_g", mode = "function")) {
        return(rustat_subset_row_sums_g(subs, ctx$g))
      }
      as.numeric(ctx$g[subs[, 1L]])
    },
    randomized = FALSE,
    symmetric = TRUE,
    ptr = NULL # pure-score layer has no compiled pointer
  )
}
