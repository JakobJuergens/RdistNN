# Shared covariance arithmetic. Callers own statistical normalization and result
# classes; this module owns centering, numeric range and symmetric assembly.
.covariance_symmetrize <- function(V) {
  for (i in seq_len(nrow(V) - 1L)) {
    j <- seq.int(i + 1L, nrow(V))
    a <- V[i, j]
    b <- V[j, i]
    same <- (a >= 0) == (b >= 0)
    mid <- a / 2 + b / 2
    mid[same & !is.na(same)] <- a[same & !is.na(same)] +
      (b[same & !is.na(same)] - a[same & !is.na(same)]) / 2
    mid[!is.na(a == b) & a == b] <- a[!is.na(a == b) & a == b]
    V[i, j] <- V[j, i] <- mid
  }
  V
}

# BLAS remains the ordinary route. Revisit cancellation-sensitive off-diagonals
# with ordered products, including the residual of mean subtraction.
.covariance_crossprod <- function(x, low) {
  V <- crossprod(x)
  norm <- sqrt(diag(V))
  for (j in seq_len(ncol(x) - 1L)) for (k in seq.int(j + 1L, ncol(x))) {
    if (abs(V[j, k]) <= sqrt(.Machine$double.eps) * norm[j] * norm[k]) {
      value <- x[, j] * x[, k]
      if (!is.null(low)) value <- c(value, x[, j] * low[, k],
        low[, j] * x[, k], low[, j] * low[, k])
      V[j, k] <- V[k, j] <- .compensated_sum(value)
    }
  }
  V
}

# Binary range normalization can erase a small cross-product even when the
# original product and final covariance are finite. Revisit only sensitive
# entries of exceptional full-covariance requests, using original-scale data.
# Retain the normalized answer when direct products need range protection.
.covariance_recover <- function(V, scaled, original, exponent, factor) {
  if (is.null(original)) return(V)
  norm <- sqrt(diag(scaled))
  columns <- new.env(parent = emptyenv())
  column <- function(j) {
    key <- as.character(j)
    cached <- columns[[key]]
    if (!is.null(cached)) return(if (cached$available) cached else NULL)
    z <- original$x[, j]
    low <- if (is.null(original$low)) numeric(length(z)) else original$low[, j]
    if (original$center) {
      z <- z - .numeric_signal_level(matrix(z, ncol = 1L))
      mean <- .compensated_sum(z) / length(z)
      if (!is.finite(mean)) {
        columns[[key]] <- list(available = FALSE)
        return(NULL)
      }
      difference <- z - mean
      back <- difference - z
      low <- (z - (difference - back)) - (mean + back)
      z <- difference
    }
    columns[[key]] <- list(high = z, low = low, available = TRUE)
    columns[[key]]
  }
  for (j in seq_len(ncol(V) - 1L)) for (k in seq.int(j + 1L, ncol(V))) {
    if ((exponent[j] == 0 && exponent[k] == 0) ||
        !is.finite(V[j, j]) || !is.finite(V[k, k]) ||
        abs(scaled[j, k]) > sqrt(.Machine$double.eps) * norm[j] * norm[k]) next
    left <- column(j)
    right <- column(k)
    if (is.null(left) || is.null(right)) next
    a <- c(left$high, left$high, left$low, left$low)
    b <- c(right$high, right$low, right$high, right$low)
    value <- a * b
    # Nonzero subnormal products can lose a finite cross-entry through their
    # separate rounding even when none becomes zero. Keep normalized recovery.
    if (any(!is.finite(value)) ||
        any(a != 0 & b != 0 & abs(value) < .Machine$double.xmin)) next
    value <- .compensated_sum(value)
    if (is.finite(value) && is.finite(factor * value)) {
      V[j, k] <- V[k, j] <- factor * value
    }
  }
  V
}

.covariance_matrix <- function(x, factor = 1, center = TRUE,
                               diagonal = FALSE, allow_failure = FALSE) {
  x <- as.matrix(x)
  n <- nrow(x)
  q <- ncol(x)
  exponent <- numeric(q)
  invalid <- colSums(!is.finite(x)) > 0L
  x[, invalid] <- 0
  # Keep ordinary jobs on vectorized cross-products. Only exceptional columns
  # require exact binary normalization before centering or squaring. The gate
  # is shared with Mata, whose finite range ends below IEEE double's maximum.
  limit <- sqrt(2^1020 / n / max(1, factor))
  amplitude <- vapply(seq_len(q), function(j) max(abs(x[, j])), numeric(1L))
  needs_scale <- amplitude >= limit |
    (amplitude > 0 & amplitude < sqrt(.Machine$double.xmin))
  original <- if (!diagonal && q > 1L && any(needs_scale)) {
    list(x = x, low = NULL, center = center)
  } else {
    NULL
  }
  for (j in seq_len(q)) {
    if (needs_scale[j]) {
      exponent[j] <- max(-1074, min(1022, floor(log(amplitude[j]) / log(2))))
      x[, j] <- x[, j] / 2^exponent[j]
    }
    if (center) x[, j] <- x[, j] - .numeric_signal_level(x[, j, drop = FALSE])
  }
  low <- NULL
  if (center) {
    means <- colMeans(x)
    if (!diagonal && q > 1L) low <- matrix(0, n, q)
    for (j in seq_len(q)) {
      difference <- x[, j] - means[j]
      if (!is.null(low)) {
        back <- difference - x[, j]
        low[, j] <- (x[, j] - (difference - back)) - (means[j] + back)
      }
      x[, j] <- difference
    }
  }
  varying <- colSums(x != 0) > 0L
  # Small differences around an ordinary offset may need normalization even
  # when the original column did not. Do this before any quadratic product.
  for (j in seq_len(q)) {
    spread <- max(abs(x[, j]))
    if (spread > 0 && spread < sqrt(.Machine$double.xmin)) {
      if (!diagonal && q > 1L && is.null(original)) {
        original <- list(x = x, low = low, center = FALSE)
      }
      extra <- max(-1074, min(1022, floor(log(spread) / log(2))))
      x[, j] <- x[, j] / 2^extra
      if (!is.null(low)) low[, j] <- low[, j] / 2^extra
      exponent[j] <- exponent[j] + extra
    }
  }
  # Use the same diagonal reduction for marginal-only and joint requests.
  variance <- factor * colSums(x * x)
  restore <- function(value, power) {
    half <- trunc(power / 2)
    (value * 2^half) * 2^(power - half)
  }
  variance <- restore(variance, 2 * exponent)
  invalid <- invalid | !is.finite(variance) | variance < 0 | (varying & variance == 0)
  if (any(invalid) && !allow_failure) {
    .dnn_abort("numerical_range",
      "Covariance is outside the numeric range; rescale the signals or request point estimates only."
    )
  }
  variance[invalid] <- NA_real_
  if (diagonal) return(unname(variance))
  scaled <- .covariance_crossprod(x, low)
  V <- restore(factor * scaled, outer(exponent, exponent, "+"))
  diag(V) <- variance
  V[invalid, ] <- NA_real_
  V[, invalid] <- NA_real_
  .covariance_recover(V, scaled, original, exponent, factor)
}
