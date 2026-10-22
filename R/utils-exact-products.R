# Exact nonnegativity of a sum of products of aligned finite binary64 vectors.
# Callers validate operands. Separate positive/negative base-2^26 limbs retain
# out-of-range products and cancel without discarding a represented residual.
.exact_sum_products_nonnegative <- function(a, b) {
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
