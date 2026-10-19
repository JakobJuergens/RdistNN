# Shared observed-reference policy for finite signal matrices. Callers put
# their preferred local observation first; zero preserves inexact differences.
.numeric_signal_level <- function(signal) {
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
