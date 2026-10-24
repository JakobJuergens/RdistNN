# A scalar correction is enough for the bounded exceptional-scale sums. Keep
# the operation order identical to Mata, including rounding each square first.
.scale_sum <- function(x) {
  total <- correction <- 0
  for (value in x) {
    next_total <- total + value
    correction <- correction + if (abs(total) >= abs(value)) {
      (total - next_total) + value
    } else {
      (value - next_total) + total
    }
    total <- next_total
  }
  total + correction
}

# Use common conservative range gates before either host's quadratic sum can
# fail. Binary normalization is exact and preserves small local differences.
.sample_sd <- function(x) {
  n <- length(x)
  amplitude <- max(abs(x))
  if (amplitude < sqrt(2^1022 / n)) {
    z <- x - .numeric_signal_level(matrix(x, ncol = 1L))
    spread <- max(abs(z - mean(z)))
    if (spread == 0) return(0)
    if (spread >= sqrt(.Machine$double.xmin) * sqrt(n)) return(stats::sd(z))
  }
  exponent <- max(-1074, min(1022, floor(log(amplitude) / log(2))))
  scale <- 2^exponent
  z <- x / scale
  z <- z - z[1L]
  z <- z - .scale_sum(z) / n
  sqrt(.scale_sum(z * z) / (n - 1L)) * scale
}

.col_sd <- function(X) apply(X, 2L, .sample_sd)
.col_mad <- function(X) apply(X, 2L, stats::mad, constant = 1.4826)
.col_iqr <- function(X) apply(X, 2L, stats::IQR)
.col_rng <- function(X) apply(X, 2L, function(z) diff(range(z)))

# Common exceptional subtraction for scaled L1/L2. Binary halving is exact
# for huge finite operands and retains difference rounding before scaling.
.scaled_distance_difference <- function(x, target, inverse) {
  d <- x - target
  recover <- abs(d) >= 2^1023 & is.finite(x) & is.finite(target) & abs(inverse) < 1
  if (any(recover)) {
    half <- x[recover] * .5 - target * .5
    quarter <- abs(half) >= 2^1023
    half[quarter] <- (x[recover][quarter] * .25 - target * .25) * (4 * inverse)
    half[!quarter] <- half[!quarter] * (inverse + inverse)
    d[recover] <- half
    d[!recover] <- d[!recover] * inverse
    return(d)
  }
  d * inverse
}

# Split restoration avoids double rounding at the final subnormal boundary.
.distance_power2 <- function(value, power) {
  if (value == 0) return(value)
  half <- trunc(power / 2)
  (value * 2^half) * 2^(power - half)
}

# Exceptional rows retain a binary exponent until after norm accumulation.
# Scratch remains scalar apart from the existing two-element decomposition.
.scaled_distance_norm_row <- function(X, i, x0, inv_s, l2, squared = FALSE) {
  scale_m <- total <- 0
  scale_e <- 0
  active <- 0L
  rounded_max <- 0
  for (j in seq_len(ncol(X))) {
    raw <- X[i, j] - x0[j]
    shift <- 0
    if (abs(raw) >= 2^1023) {
      raw <- X[i, j] * .25 - x0[j] * .25
      shift <- 2
    }
    if (raw == 0 || inv_s[j] == 0) next
    if (!is.finite(raw) || !is.finite(inv_s[j])) return(abs(raw * inv_s[j]))
    active <- active + 1L
    rounded_max <- max(rounded_max, abs(.scaled_distance_difference(X[i, j], x0[j], inv_s[j])))
    parts <- .numeric_binary_parts(c(abs(raw), abs(inv_s[j])))
    mantissa <- parts$value
    exponent <- parts$power
    while (any(mantissa >= 1)) {
      upper <- mantissa >= 1
      mantissa[upper] <- mantissa[upper] / 2
      exponent[upper] <- exponent[upper] + 1
    }
    m <- mantissa[1L] * mantissa[2L]
    e <- exponent[1L] + exponent[2L] + shift
    if (m < .5) {
      m <- m * 2
      e <- e - 1
    }
    if (scale_m == 0) {
      scale_m <- m
      scale_e <- e
      total <- 1
    } else if (e > scale_e || (e == scale_e && m > scale_m)) {
      ratio <- .distance_power2(scale_m / m, scale_e - e)
      square <- ratio * ratio
      term <- total * if (l2) square else ratio
      total <- 1 + term
      scale_m <- m
      scale_e <- e
    } else {
      ratio <- .distance_power2(m / scale_m, e - scale_e)
      square <- ratio * ratio
      total <- total + if (l2) square else ratio
    }
  }
  if (scale_m == 0) return(0)
  if (active == 1L) return(if (squared) rounded_max * rounded_max else rounded_max)
  if (squared) return(.distance_power2((scale_m * total) * scale_m, 2 * scale_e))
  max(.distance_power2(scale_m * if (l2) sqrt(total) else total, scale_e), rounded_max)
}

.assert_no_na <- function(X, x0 = NULL) {
  X_check <- if (is.data.frame(X)) as.matrix(X) else X
  if (!is.numeric(X_check)) {
    .dnn_abort("argument_invalid", "`X` must be numeric.", argument = "X")
  }
  if (.has_na_only(X_check)) {
    .dnn_abort("argument_invalid", "NA values are not supported: `X` contains NA.", argument = "X")
  }
  if (.has_nonfinite_non_na(X_check)) {
    .dnn_abort("argument_invalid",
      "Non-finite values are not supported: `X` contains Inf, -Inf, or NaN.",
      argument = "X"
    )
  }
  if (!is.null(x0)) {
    if (!is.numeric(x0)) {
      .dnn_abort("argument_invalid", "`x0` must be numeric.", argument = "x0")
    }
    if (.has_na_only(x0)) {
      .dnn_abort("argument_invalid", "NA values are not supported: `x0` contains NA.", argument = "x0")
    }
    if (.has_nonfinite_non_na(x0)) {
      .dnn_abort("argument_invalid",
        "Non-finite values are not supported: `x0` contains Inf, -Inf, or NaN.",
        argument = "x0"
      )
    }
  }
  invisible(NULL)
}

.coerce_inv_s <- function(scale_vec, p, eps) {
  if (!is.numeric(scale_vec)) {
    .dnn_abort("argument_invalid", "`scale_vec` must be numeric when rule = 'custom'.", argument = "scale_vec")
  }
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    .dnn_abort("argument_invalid", "`eps` must be a positive finite numeric scalar.", argument = "eps")
  }
  inv_s <- as.numeric(scale_vec)
  if (length(inv_s) != p) .dnn_abort(
    "argument_invalid",
    "`scale_vec` must be length p when rule = 'custom'.",
    argument = "scale_vec"
  )
  if (.has_na_only(inv_s)) {
    .dnn_abort("argument_invalid", "`scale_vec` cannot contain NA values.", argument = "scale_vec")
  }
  if (.has_nonfinite_non_na(inv_s)) {
    .dnn_abort("argument_invalid",
      "`scale_vec` cannot contain non-finite values (`Inf`, `-Inf`, or `NaN`).",
      argument = "scale_vec"
    )
  }
  1 / pmax(abs(inv_s), eps)
}

.compute_inv_scale <- function(X, rule = c("sd", "mad", "iqr", "range", "custom"),
                               scale_vec = NULL, eps = 1e-12) {
  rule <- .dnn_match_arg(rule)
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    .dnn_abort("argument_invalid", "`eps` must be a positive finite numeric scalar.", argument = "eps")
  }
  if (identical(rule, "custom")) {
    return(.coerce_inv_s(scale_vec, ncol(X), eps))
  }
  if (identical(rule, "sd") && nrow(X) < 2L) {
    .dnn_abort("argument_invalid", "Scaling rule 'sd' requires at least two training observations.", argument = "X")
  }

  s <- switch(rule,
    sd = .col_sd(X),
    mad = .col_mad(X),
    iqr = .col_iqr(X),
    range = .col_rng(X)
  )
  if (length(s) != ncol(X)) .dnn_abort("internal_invariant",
    "Scale vector length must match ncol(X).",
    hint = "Report a reproducible example of this internal inconsistency."
  )
  if (any(!is.finite(s))) {
    .dnn_abort("numerical_range", "Fitted scales are nonfinite; rescale covariates.")
  }
  s <- pmax(s, eps)
  1 / s
}
