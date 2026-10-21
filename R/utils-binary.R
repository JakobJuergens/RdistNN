# Exact binary decomposition for bounded range recoveries. Nonzero values
# normalize to magnitude [1,2), except the shared upper cap permits values <4.
.numeric_binary_parts <- function(value) {
  power <- numeric(length(value))
  nonzero <- value != 0
  power[nonzero] <- pmax(-1074, pmin(1022,
    floor(log(abs(value[nonzero])) / log(2))))
  mantissa <- value / 2^power
  # Correct log rounding at exact binary boundaries to match native ilogb.
  lower <- nonzero & abs(mantissa) < 1 & power > -1074
  while (any(lower)) {
    power[lower] <- power[lower] - 1
    mantissa[lower] <- mantissa[lower] * 2
    lower <- nonzero & abs(mantissa) < 1 & power > -1074
  }
  upper <- nonzero & abs(mantissa) >= 2 & power < 1022
  while (any(upper)) {
    power[upper] <- power[upper] + 1
    mantissa[upper] <- mantissa[upper] / 2
    upper <- nonzero & abs(mantissa) >= 2 & power < 1022
  }
  list(value = mantissa, power = power)
}
