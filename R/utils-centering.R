# Shared observed-reference policy for finite signal matrices. Callers put
# their preferred local observation first; zero preserves inexact differences.
.numeric_signal_level <- function(signal, weights = NULL, weight_mass = NULL) {
  if (!nrow(signal)) return(rep(0, ncol(signal)))
  lo <- apply(signal, 2L, min)
  hi <- apply(signal, 2L, max)
  level <- signal[1L, ]
  # Prefer a common observed origin even for mixed signs. The conservative
  # range gate is shared with Mata and avoids overflowing extreme differences.
  safe <- is.finite(hi - level) & is.finite(lo - level) &
    abs(hi - level) < 2^1023 & abs(lo - level) < 2^1023
  fallback <- ifelse(lo > 0, lo, ifelse(hi < 0, hi, 0))
  level <- ifelse(safe, level, fallback)
  # A translated weighted residual can underflow before its subnormal common
  # origin is restored. Zero preserves direct separately rounded products.
  amplitude <- pmax(abs(lo), abs(hi))
  level[amplitude < .Machine$double.xmin] <- 0
  # When every weighted signal product is subnormal, translating a normal
  # origin can create variation even for zero survivors or constant signals.
  # A maximum-weight bound leaves ordinary jobs with tiny tails unchanged.
  # A subnormal nonzero mass also uses original products, including when
  # large signals amplify rounded block coefficients. Exact zero mass stays
  # centered. Recomputed callers may pass an explicit mass summary.
  if (!is.null(weights)) {
    weight_amplitude <- if (length(weights)) max(abs(weights)) else 0
    if (is.null(weight_mass)) weight_mass <- .compensated_sum(weights)
    underflow <- weight_amplitude * amplitude < .Machine$double.xmin
    if (weight_mass != 0) {
      underflow <- underflow | abs(weight_mass) < .Machine$double.xmin |
        abs(weight_mass) * amplitude < .Machine$double.xmin
    }
    level[underflow] <- 0
  }
  for (j in seq_len(ncol(signal))) {
    difference <- signal[, j] - level[j]
    back <- difference - signal[, j]
    error <- (signal[, j] - (difference - back)) - (level[j] + back)
    # Do not erase small, independently represented signals while translating
    # a mixed-sign column (for example 2^54, 1, -2^54).
    if (any(error != 0)) level[j] <- 0
  }
  level
}
