# Shared finite double-precision accumulation. Inputs are canonical numeric
# arrays validated by callers. Keep this arithmetic aligned with the native
# header and the Stata companion; covariance policy belongs to its own module.
.compensated_state_R <- function(x) {
  high <- low <- 0
  for (value in x) {
    next_high <- high + value
    if (!is.finite(next_high)) {
      high <- next_high
      low <- 0
      next
    }
    low <- low + if (abs(high) >= abs(value)) {
      (high - next_high) + value
    } else {
      (value - next_high) + high
    }
    high <- next_high
  }
  c(high, low)
}

.compensated_sum_R <- function(x) {
  state <- .compensated_state_R(x)
  state[1L] + state[2L]
}

.compensated_sum <- function(x) {
  cpp <- .rdistnn_compiled_fun("rustat_compensated_sum")
  if (!is.null(cpp)) return(cpp(as.numeric(x)))
  .compensated_sum_R(x)
}

.compensated_crossprod_R <- function(x, y) {
  out <- matrix(0, ncol(x), ncol(y))
  for (j in seq_len(ncol(y))) {
    for (k in seq_len(ncol(x))) out[k, j] <- .compensated_sum_R(x[, k] * y[, j])
  }
  out
}

.compensated_crossprod <- function(x, y) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  if (nrow(x) != nrow(y)) .abort_bad_state("Weighted signals require aligned rows.")
  cpp <- .rdistnn_compiled_fun("rustat_compensated_crossprod")
  out <- if (is.null(cpp)) .compensated_crossprod_R(x, y) else cpp(x, y)
  if (!is.null(colnames(x)) || !is.null(colnames(y))) {
    dimnames(out) <- list(colnames(x), colnames(y))
  }
  out
}

# Sensitive projections need the same dot-product order for every row shape.
# Reuse the native reduction; vectorize independent rows in the R fallback.
.compensated_matprod <- function(x, y) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  if (ncol(x) != nrow(y)) .abort_bad_state("Matrix products require aligned inner dimensions.")
  cpp <- .rdistnn_compiled_fun("rustat_compensated_crossprod")
  if (!is.null(cpp)) {
    out <- cpp(t(x), y)
  } else {
    out <- matrix(0, nrow(x), ncol(y))
    for (j in seq_len(ncol(y))) {
      high <- low <- numeric(nrow(x))
      for (k in seq_len(ncol(x))) {
        state <- .compensated_add_R(high, low, x[, k] * y[k, j])
        high <- state$high
        low <- state$low
      }
      out[, j] <- high + low
    }
  }
  if (!is.null(rownames(x)) || !is.null(colnames(y))) {
    dimnames(out) <- list(rownames(x), colnames(y))
  }
  out
}

# Vector update of independent compensated states; also used by rank tables.
.compensated_add_R <- function(high, low, value) {
  next_high <- high + value
  finite <- is.finite(next_high)
  first <- finite & abs(high) >= abs(value)
  second <- finite & !first
  error <- numeric(length(high))
  value <- rep_len(value, length(high))
  error[first] <- (high[first] - next_high[first]) + value[first]
  error[second] <- (value[second] - next_high[second]) + high[second]
  low <- low + error
  low[!finite] <- 0
  list(high = next_high, low = low)
}

# Merge a bounded block of states, in the same order as the native helper.
.compensated_merge_R <- function(left_high, right_high, left_low, right_low, offset) {
  high <- low <- numeric(length(left_high))
  for (value in list(left_high, right_high, left_low, right_low, offset)) {
    state <- .compensated_add_R(high, low, value)
    high <- state$high
    low <- state$low
  }
  high + low
}

# Two n-element forward arrays and fixed-size reverse/merge buffers. This
# avoids full-length temporary vectors for every final compensated addition.
.compensated_deletions_R <- function(z, table, offset) {
  n <- length(z)
  high <- low <- numeric(n)
  left_high <- left_low <- 0
  if (n > 1L) for (r in 2L:n) {
    value <- table[r - 1L] * z[r - 1L]
    next_high <- left_high + value
    if (is.finite(next_high)) {
      left_low <- left_low + if (abs(left_high) >= abs(value)) {
        (left_high - next_high) + value
      } else {
        (value - next_high) + left_high
      }
    } else {
      left_low <- 0
    }
    left_high <- next_high
    high[r] <- left_high
    low[r] <- left_low
  }
  right_high <- right_low <- 0
  end <- n
  while (end > 0L) {
    start <- max(1L, end - 2047L)
    indices <- seq.int(start, end)
    suffix_high <- suffix_low <- numeric(length(indices))
    for (r in seq.int(end, start)) {
      at <- r - start + 1L
      suffix_high[at] <- right_high
      suffix_low[at] <- right_low
      if (r > 1L) {
        value <- table[r - 1L] * z[r]
        next_high <- right_high + value
        if (is.finite(next_high)) {
          right_low <- right_low + if (abs(right_high) >= abs(value)) {
            (right_high - next_high) + value
          } else {
            (value - next_high) + right_high
          }
        } else {
          right_low <- 0
        }
        right_high <- next_high
      }
    }
    high[indices] <- .compensated_merge_R(high[indices], suffix_high,
      low[indices], suffix_low, offset)
    end <- start - 1L
  }
  high
}
