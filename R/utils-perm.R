#' Convert an index order to ranks (inverse permutation)
#'
#' @param ord Integer vector (values `1..n`) giving increasing order.
#' @return Integer ranks vector `r` where `r[i]` is the rank of item `i`.
#' @examples
#' ord_to_ranks(c(3L, 1L, 2L))
#' @keywords internal
#' @noRd
ord_to_ranks <- function(ord) {
  ord <- as.integer(ord)
  n <- length(ord)
  if (!setequal(ord, seq_len(n))) stop("`ord` must be a permutation of 1..n.", call. = FALSE)
  r <- integer(n)
  r[ord] <- seq_len(n)
  r
}

# Canonical numeric row keys use original covariates in their supplied column
# order. Observation position resolves only rows identical on every covariate.
.covariate_tie_key <- function(X, row_ids = seq_len(nrow(X))) {
  columns <- unname(as.list(as.data.frame(X)))
  ord_to_ranks(do.call(order, c(columns, list(row_ids, method = "radix"))))
}

.order_with_ties <- function(dx, X = NULL, tie_key = NULL) {
  if (is.null(tie_key)) {
    tie_key <- if (is.null(X)) seq_along(dx) else .covariate_tie_key(X)
  }
  order(dx, tie_key, method = "radix")
}
