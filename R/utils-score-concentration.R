# Signed weighted residuals and their five-largest absolute-mass share.
# right is already rounded by the caller. Preserve ordinary finite arithmetic;
# exceptional finite-operand differences retain magnitude in binary exponents.
.dnn_score_contributions <- function(w, left, right) {
  if (length(w) != length(left) || length(w) != length(right)) {
    .abort_bad_state("Score concentration requires aligned weights and residual operands.")
  }
  contribution <- numeric(length(w))
  result <- function(available, share = NA_real_) {
    list(available = available, share = share, contribution = contribution)
  }
  if (any(!is.finite(w)) || any(w < 0)) {
    contribution[] <- NA_real_
    return(result(FALSE))
  }
  active <- which(w > 0)
  if (!length(active)) return(result(TRUE))
  weights <- w[active]
  l <- left[active]
  r <- right[active]
  valid <- is.finite(l) & is.finite(r)
  signed <- rep(NA_real_, length(active))
  signed[valid] <- weights[valid] * (l[valid] - r[valid])
  contribution[active] <- signed
  magnitude <- abs(signed)
  total <- sum(magnitude)
  top_share <- function(value, mass) {
    sum(sort(value, decreasing = TRUE)[seq_len(min(5L, length(value)))]) / mass
  }
  if (all(is.finite(magnitude)) && is.finite(total)) {
    if (total == 0) return(result(TRUE))
    return(result(TRUE, top_share(magnitude, total)))
  }

  mantissa <- exponent <- numeric(length(weights))
  ordinary <- which(is.finite(magnitude) & magnitude > 0)
  if (length(ordinary)) {
    parts <- .numeric_binary_parts(magnitude[ordinary])
    mantissa[ordinary] <- parts$value
    exponent[ordinary] <- parts$power
  }
  for (i in which(valid & !is.finite(signed))) {
    parts <- .numeric_binary_parts(c(l[i], r[i]))
    power <- max(parts$power)
    if (l[i] == 0) power <- parts$power[2L]
    if (r[i] == 0) power <- parts$power[1L]
    difference <- l[i] / 2^power - r[i] / 2^power
    weight <- .numeric_binary_parts(weights[i])
    value <- difference * weight$value
    if (!is.finite(value) || value == 0) next
    parts <- .numeric_binary_parts(value)
    mantissa[i] <- abs(parts$value)
    exponent[i] <- power + weight$power + parts$power
    # Keep the first multiplication normal when the final result is tiny,
    # so a subnormal result is rounded only by the final multiplication.
    first_power <- max(-1022, min(1022, exponent[i]))
    restored <- (parts$value * 2^first_power) * 2^(exponent[i] - first_power)
    contribution[active[i]] <- if (is.finite(restored)) restored else NA_real_
  }
  if (any(!valid) || any(valid & !is.finite(signed) & mantissa == 0)) {
    return(result(FALSE))
  }
  nonzero <- which(mantissa > 0)
  if (!length(nonzero)) return(result(TRUE))
  delta <- exponent - max(exponent[nonzero])
  keep <- nonzero[delta[nonzero] >= -1074]
  magnitude <- numeric(length(weights))
  magnitude[keep] <- mantissa[keep] * 2^delta[keep]
  result(TRUE, top_share(magnitude, sum(magnitude)))
}
