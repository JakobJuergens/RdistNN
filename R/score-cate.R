#' CATE score from precomputed pseudo-outcome rho (pure score; no weights)
#'
#' This is a pure-score implementation. It computes the per-observation
#' contributions \eqn{g_i = \rho_i} and stores them in `ctx$g`.
#' No ordering, no DNN/kernels, no L-stat weights are touched here.
#' Any weighting/aggregation happens outside the score layer.
#'
#' @param x0 numeric length-p target point (checked for dimension only; not used here).
#' @param x_cols names/indices of X columns in `data` (used for dimension checks).
#' @param s integer subset size (reported as the kernel order; not used internally).
#' @param k integer (present only for API compatibility; ignored by this pure score).
#' @param rho_col optional name of a column in `data` holding rho.
#' @param rho optional numeric vector (length n) of rho; used if `rho_col` is NULL.
#'
#' @return An `ustat_score` whose setup places `g` in the context; `eval` and
#'         `eval_batch` only index/sum `g` (no weights).
#' @export
score_cate_rho <- function(
    x0, x_cols, s, k = 1L,
    rho_col = NULL, rho = NULL
) {
  # XOR: exactly one of rho_col or rho
  if (is.null(rho_col) == is.null(rho)) {
    stop("Provide exactly one of 'rho_col' or 'rho' (mutually exclusive).")
  }
  if (!is.numeric(s) || length(s) != 1L || s < 1L) {
    stop("`s` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(k) || length(k) != 1L || k < 1L) {
    stop("`k` must be a positive integer.", call. = FALSE)
  }

  new_score(
    name  = "cate_rho_purescore",
    order = function(n, data, ctx) as.integer(s),

    # in score_cate_rho(...)

    setup = function(data, ...) {
      X <- as.matrix(data[, x_cols, drop = FALSE])
      r <- if (!is.null(rho_col)) as.numeric(data[[rho_col]]) else as.numeric(rho)

      if (length(x0) != ncol(X)) stop("x0 length must equal ncol(X).", call. = FALSE)
      if (length(r)  != nrow(X)) stop("rho length must equal nrow(data).", call. = FALSE)

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
      rowSums(matrix(gi[subs], nrow = nrow(subs)))
    },


    randomized = FALSE,
    symmetric  = TRUE,
    ptr        = NULL           # no compiled pointer in the pure-score layer
  )
}
