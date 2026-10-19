# Canonical exact rank law shared with C++ and Mata. Normalize each
# order-statistic PMF at its mode, then average and normalize once in rank order.
# Recurrence state retains a binary exponent after emitted values underflow.
.dnn_exact_component_R <- function(n, s, l) {
  p <- numeric(n - s + 1L)
  last <- n - s + l
  mode <- max(l, min(last, floor((as.double(n) * (l - 1)) / (s - 1)) + 1))
  p[mode - l + 1L] <- 1
  for (direction in c(1, -1)) {
    mantissa <- .5
    exponent <- 1
    rank <- mode
    end <- if (direction == 1) last else l
    while (rank != end) {
      ratio <- if (direction == 1) {
        (rank / (rank - l + 1)) * ((n - rank - s + l) / (n - rank))
      } else {
        ((rank - l) / (rank - 1)) * ((n - rank + 1) / (n - rank - s + l + 1))
      }
      mantissa <- mantissa * ratio
      while (mantissa < .5) {
        mantissa <- mantissa * 2
        exponent <- exponent - 1
      }
      while (mantissa >= 1) {
        mantissa <- mantissa / 2
        exponent <- exponent + 1
      }
      rank <- rank + direction
      p[rank - l + 1L] <- if (exponent < -1074) 0 else mantissa * 2^exponent
    }
  }
  p / .compensated_sum_R(p)
}

.dnn_exact_table_R <- function(n, k, s) {
  if (k == s) return(rep(1 / n, n))
  high <- low <- numeric(n)
  for (l in seq_len(k)) {
    value <- .dnn_exact_component_R(n, s, l)
    at <- seq.int(l, n - s + l)
    state <- .compensated_add_R(high[at], low[at], value)
    high[at] <- state$high
    low[at] <- state$low
  }
  out <- (high + low) / k
  out / .compensated_sum_R(out)
}

.dnn_exact_table <- function(n, k, s) {
  cpp <- .rdistnn_compiled_fun("rustat_exact_rank_table")
  if (is.null(cpp)) .dnn_exact_table_R(n, k, s) else cpp(n, k, s)
}
