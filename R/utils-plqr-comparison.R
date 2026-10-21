# Exact sign of sum(D * theta) + g - y for finite binary64 inputs. Separate
# positive/negative base-2^26 limbs retain products too large or small for a
# double and cancel without discarding a residual. This is the same bounded
# integer arithmetic used by the native and Stata comparison owners.
.plqr_exact_comparison <- function(D, theta, g, y) {
  a <- c(D, g, y)
  b <- c(theta, 1, -1)
  if (any(!is.finite(a)) || any(!is.finite(b))) {
    .abort_bad_state("PLQR threshold comparisons require finite inputs.")
  }
  left <- .numeric_binary_parts(a)
  right <- .numeric_binary_parts(b)
  base <- 2^26
  origin <- 2 * (-1074) - 2 * 52
  # Finite products are below 2^2048; at most INT_MAX+2 terms sum below
  # 2^2080. Three spare limbs keep the fixed storage bound explicit.
  limb_count <- ceiling((2048 + 32 - origin) / 26) + 3L
  limbs <- matrix(0, limb_count, 2L)
  for (k in which(a != 0 & b != 0)) {
    x <- abs(left$value[k]) * 2^52
    z <- abs(right$value[k]) * 2^52
    xs <- zs <- numeric(3L)
    for (i in seq_len(3L)) {
      next_x <- floor(x / base)
      next_z <- floor(z / base)
      xs[i] <- x - next_x * base
      zs[i] <- z - next_z * base
      x <- next_x
      z <- next_z
    }
    side <- if ((a[k] < 0) == (b[k] < 0)) 1L else 2L
    start <- left$power[k] + right$power[k] - 104 - origin
    for (i in seq_len(3L)) for (j in seq_len(3L)) {
      product <- xs[i] * zs[j]
      high <- floor(product / base)
      pieces <- c(product - high * base, high)
      bit <- start + 26 * (i + j - 2L)
      for (part in seq_len(2L)) {
        if (pieces[part] == 0) {
          bit <- bit + 26
          next
        }
        index <- floor(bit / 26) + 1L
        piece <- pieces[part] * 2^(bit %% 26)
        while (piece != 0) {
          value <- limbs[index, side] + piece
          carry <- floor(value / base)
          limbs[index, side] <- value - carry * base
          piece <- carry
          index <- index + 1L
        }
        bit <- bit + 26
      }
    }
  }
  different <- which(limbs[, 1L] != limbs[, 2L])
  if (!length(different)) return(TRUE)
  top <- different[length(different)]
  limbs[top, 1L] > limbs[top, 2L]
}

# Keep the ordinary matrix-product comparison only when a conservative
# forward-error envelope separates its margin from zero. Range-risk and
# uncertain rows use the exact sign, including inclusive exact equality.
.plqr_below_threshold <- function(D, theta, g, y) {
  threshold <- as.vector(D %*% theta) + g
  bound <- abs(g) + abs(y)
  for (j in seq_len(ncol(D))) {
    product <- D[, j] * theta[j]
    bound <- bound + abs(product)
  }
  L <- ncol(D) + 2
  envelope <- 8 * (L + 1) * .Machine$double.eps * bound + 32 * L * 2^-1074
  margin <- threshold - y
  exceptional <- !is.finite(bound) | bound >= 2^1022 | !is.finite(threshold) |
    !is.finite(margin) | abs(margin) <= envelope
  below <- y <= threshold
  for (i in which(exceptional)) {
    below[i] <- .plqr_exact_comparison(D[i, ], theta, g[i], y[i])
  }
  below
}
