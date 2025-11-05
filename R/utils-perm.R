# R/permutation-utils.R
# -----------------------------------------------------------------------------

#' Convert an index order to ranks (inverse permutation)
#'
#' @param ord Integer vector (values `1..n`) giving increasing order.
#' @return Integer ranks vector `r` where `r[i]` is the rank of item `i`.
#' @export
ord_to_ranks <- function(ord) {
  ord <- as.integer(ord)
  n <- length(ord)
  if (!setequal(ord, seq_len(n))) stop("`ord` must be a permutation of 1..n.", call. = FALSE)
  r <- integer(n)
  r[ord] <- seq_len(n)
  r
}
