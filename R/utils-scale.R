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
    spread <- max(abs(x - mean(x)))
    if (spread == 0) return(0)
    if (spread >= sqrt(.Machine$double.xmin) * sqrt(n)) return(stats::sd(x))
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
    d[recover] <- (x[recover] * .5 - target * .5) * (inverse + inverse)
    d[!recover] <- d[!recover] * inverse
    return(d)
  }
  d * inverse
}

.assert_no_na <- function(X, x0 = NULL) {
  X_check <- if (is.data.frame(X)) as.matrix(X) else X
  if (!is.numeric(X_check)) {
    stop("`X` must be numeric.", call. = FALSE)
  }
  if (.has_na_only(X_check)) {
    stop("NA values are not supported: `X` contains NA.", call. = FALSE)
  }
  if (.has_nonfinite_non_na(X_check)) {
    stop("Non-finite values are not supported: `X` contains Inf, -Inf, or NaN.", call. = FALSE)
  }
  if (!is.null(x0)) {
    if (!is.numeric(x0)) {
      stop("`x0` must be numeric.", call. = FALSE)
    }
    if (.has_na_only(x0)) {
      stop("NA values are not supported: `x0` contains NA.", call. = FALSE)
    }
    if (.has_nonfinite_non_na(x0)) {
      stop("Non-finite values are not supported: `x0` contains Inf, -Inf, or NaN.", call. = FALSE)
    }
  }
  invisible(NULL)
}

.coerce_inv_s <- function(scale_vec, p, eps) {
  if (!is.numeric(scale_vec)) {
    stop("scale_vec must be numeric when rule = 'custom'.", call. = FALSE)
  }
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("eps must be a positive finite numeric scalar.", call. = FALSE)
  }
  inv_s <- as.numeric(scale_vec)
  if (length(inv_s) != p) stop("scale_vec must be length p when rule = 'custom'.", call. = FALSE)
  if (.has_na_only(inv_s)) {
    stop("scale_vec cannot contain NA values.", call. = FALSE)
  }
  if (.has_nonfinite_non_na(inv_s)) {
    stop("scale_vec cannot contain non-finite values (`Inf`, `-Inf`, or `NaN`).", call. = FALSE)
  }
  1 / pmax(abs(inv_s), eps)
}

.compute_inv_scale <- function(X, rule = c("sd", "mad", "iqr", "range", "custom"),
                               scale_vec = NULL, eps = 1e-12) {
  rule <- match.arg(rule)
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("eps must be a positive finite numeric scalar.", call. = FALSE)
  }
  if (identical(rule, "custom")) {
    return(.coerce_inv_s(scale_vec, ncol(X), eps))
  }
  if (identical(rule, "sd") && nrow(X) < 2L) {
    stop("Scaling rule 'sd' requires at least two training observations.", call. = FALSE)
  }

  s <- switch(rule,
    sd = .col_sd(X),
    mad = .col_mad(X),
    iqr = .col_iqr(X),
    range = .col_rng(X)
  )
  if (length(s) != ncol(X)) stop("Scale vector length must match ncol(X).", call. = FALSE)
  if (any(!is.finite(s))) {
    stop("Fitted scales are nonfinite; rescale covariates.", call. = FALSE)
  }
  s <- pmax(s, eps)
  1 / s
}

# Lightweight, readable column scaling (yes, it copies; that’s OK for simplicity)
.scale_cols <- function(M, inv_s) {
  sweep(M, 2L, as.numeric(inv_s), `*`)
}
