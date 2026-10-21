# Independent scalar references for diagnostics contracts.
.diag_top_share <- function(w, k) {
  k_use <- max(1L, min(length(w), as.integer(k)))
  sum(sort(as.numeric(w), decreasing = TRUE)[seq_len(k_use)])
}

.diag_weighted_median <- function(x, w) {
  .diag_weighted_quantile(x, w, prob = 0.5)
}

.diag_weighted_quantile <- function(x, w, prob = 0.5) {
  ord <- order(x)
  x_ord <- x[ord]
  w_ord <- w[ord]
  cw <- cumsum(w_ord) / sum(w_ord)
  x_ord[which(cw >= prob)[1L]]
}

.diag_weighted_sd <- function(x, w) {
  mu <- sum(x * w) / sum(w)
  sqrt(sum(w * (x - mu)^2) / sum(w))
}
