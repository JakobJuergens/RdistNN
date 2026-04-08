#' CATE score from a precomputed pseudo-outcome \eqn{\rho} (pure score; no weights)
#'
#' Constructs a **pure-score** kernel for CATE experiments that uses only the
#' precomputed pseudo-outcome \eqn{\rho}. Per observation \eqn{i}, the stored
#' contribution is \deqn{g_i \equiv \rho_i\,,} and for a subset
#' \eqn{S \subset \{1,\ldots,n\}} of size \eqn{s} the kernel evaluates to the
#' **sum of contributions**
#' \deqn{h(S) \;=\; \sum_{i \in S} g_i \;=\; \sum_{i \in S} \rho_i\,.}
#'
#' This score does **not** perform any nearest-neighbor ordering, kernel
#' weighting, or L-stat aggregation. Any weighting/averaging of subset
#' evaluations is handled outside the score layer (e.g., by the rank-/distance-
#' weight modules). The arguments \code{x0}, \code{x_cols}, \code{s}, \code{k}
#' are accepted for API consistency and basic dimension checks; only \eqn{\rho}
#' is used to construct the score.
#'
#' @param x0 Numeric length-\eqn{p} target point. Used only to check that
#'   \code{length(x0)} matches \code{length(x_cols)}; it is not used in the
#'   kernel itself.
#' @param x_cols Integer or character indices of the \eqn{X} columns in
#'   \code{data}, used for dimension checks with \code{x0}.
#' @param s Integer subset size (kernel order). Recorded in the score metadata;
#'   the pure kernel internally computes \eqn{\sum_{i \in S} \rho_i} for any
#'   subset \eqn{S} of size \eqn{s}.
#' @param k Integer included for interface compatibility with other scores; it
#'   is ignored here.
#' @param rho_col Optional name of the column in \code{data} that contains the
#'   pseudo-outcome \eqn{\rho}. If \code{NULL}, supply \code{rho}.
#' @param rho Optional numeric vector (length \eqn{n}) of \eqn{\rho} values to
#'   use when \code{rho_col} is \code{NULL}.
#'
#' @return An \code{ustat_score} object with:
#' \itemize{
#'   \item \code{$name}: a descriptive name,
#'   \item \code{$order}: the kernel order \eqn{s},
#'   \item \code{$setup}: stores \code{g <- rho} in the context,
#'   \item \code{$eval}/\code{$eval_batch}: return \eqn{\sum_{i \in S} g_i}
#'         for each subset \eqn{S} (no weights applied).
#' }
#'
#' @examples
#' demo <- cate_demo_data(n = 200, p = 2, design = "hetero", seed = 1)
#' sc <- score_cate_rho(
#'   x0 = rep(0, 2), x_cols = paste0("X", 1:2),
#'   s = 5, k = 1, rho_col = "rho_oracle"
#' )
#' ctx <- score_setup(sc, demo$data)
#' subs <- cate_demo_subsets(n = nrow(demo$data), K = 3, s = 5, seed = 1)
#' score_eval_batch(sc, subs, demo$data, ctx) # three subset evaluations
#'
#' @export
score_cate_rho <- function(
  x0, x_cols, s, k = 1L,
  rho_col = NULL, rho = NULL
) {
  # XOR: exactly one of rho_col or rho
  if (is.null(rho_col) == is.null(rho)) {
    stop("Provide exactly one of 'rho_col' or 'rho' (mutually exclusive).", call. = FALSE)
  }
  s <- .validate_scalar_int(s, arg = "s", min = 1L)
  k <- .validate_scalar_int(k, arg = "k", min = 1L)

  new_score(
    name = "cate_rho_purescore",
    order = function(n, data, ctx) as.integer(s),

    # in score_cate_rho(...)

    setup = function(data, ...) {
      X <- .dnn_col_mat(x_cols, data = data, arg = "x_cols", out_arg = "X")
      n <- nrow(X)
      r <- if (!is.null(rho_col)) {
        .dnn_col_vec(rho_col, data = data, n = n, arg = "rho_col", ref = "data")
      } else {
        .validate_vector_nrow(rho, n, arg = "rho", ref = "data")
      }

      .validate_target_vector(x0, ncol(X), x0_arg = "x0", x_arg = "X")

      # DO NOT check finiteness here (tests expect failure at evaluation)
      list(g = r, n = nrow(X), p = ncol(X), s = as.integer(s))
    },

    # keep these guards (they make it fail at evaluation time):
    eval = function(idxs, data, ctx, rng) {
      gi <- ctx$g[idxs]
      if (any(!is.finite(gi))) stop("Non-finite g encountered.", call. = FALSE)
      as.numeric(sum(gi))
    },
    eval_batch = function(subs, data, ctx, rng) {
      gi <- ctx$g
      if (any(!is.finite(gi))) stop("Non-finite g encountered.", call. = FALSE)
      if (!is.matrix(subs)) subs <- as.matrix(subs)
      if (nrow(subs) == 0L) {
        return(numeric(0L))
      }

      if (exists("rustat_subset_row_sums_g", mode = "function")) {
        return(rustat_subset_row_sums_g(subs, gi))
      }
      rowSums(matrix(gi[subs], nrow = nrow(subs), ncol = ncol(subs)))
    },
    randomized = FALSE,
    symmetric = TRUE,
    ptr = NULL # no compiled pointer in the pure-score layer
  )
}
