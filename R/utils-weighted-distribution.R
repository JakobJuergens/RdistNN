# Canonical discrete CDF decisions: sort at the caller, then use the same
# Neumaier total and prefix arithmetic in C++, R and Mata. No boundary epsilon.
.weighted_cdf_index_R <- function(w, prob) {
  mass <- .compensated_sum_R(w)
  if (!is.finite(mass) || mass <= 0) return(NA_integer_)
  high <- low <- 0
  for (i in seq_along(w)) {
    next_high <- high + w[i]
    low <- low + if (abs(high) >= abs(w[i])) {
      (high - next_high) + w[i]
    } else {
      (w[i] - next_high) + high
    }
    high <- next_high
    if ((high + low) / mass >= prob) return(i)
  }
  length(w)
}

.weighted_cdf_index <- function(w, prob) {
  cpp <- .rdistnn_compiled_fun("rustat_weighted_cdf_index")
  if (!is.null(cpp)) return(cpp(w, prob))
  .weighted_cdf_index_R(w, prob)
}

# Population weighted SD, computed directly as a norm. Its square need not be
# representable. Binary scaling and an observed origin preserve small variation
# around a large level; inactive observations must not determine that scale.
.weighted_sd <- function(x, w) {
  active <- w > 0
  if (any(w < 0) || !any(active)) return(NA_real_)
  x <- x[active]
  w <- w[active]
  mass <- .compensated_sum(w)
  if (!is.finite(mass) || mass <= 0) return(NA_real_)
  amplitude <- max(abs(x))
  if (amplitude == 0) return(0)
  scale <- 2^max(-1074, min(1022, floor(log(amplitude) / log(2))))
  z <- x / scale
  z <- z - z[1L]
  weights <- w / mass
  center <- .compensated_sum(z * weights)
  d <- sqrt(weights) * (z - center)
  spread <- max(abs(d))
  if (spread == 0) return(0)
  z <- d / spread
  norm <- spread * sqrt(.compensated_sum(z * z))
  norm * scale
}
