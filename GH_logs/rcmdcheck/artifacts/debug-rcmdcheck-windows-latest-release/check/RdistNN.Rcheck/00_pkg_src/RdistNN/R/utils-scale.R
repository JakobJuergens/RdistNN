.col_sd <- function(X) apply(X, 2L, stats::sd)
.col_mad <- function(X) apply(X, 2L, stats::mad, constant = 1.4826)
.col_iqr <- function(X) apply(X, 2L, stats::IQR)
.col_rng <- function(X) apply(X, 2L, function(z) diff(range(z)))

.assert_no_na <- function(X, x0 = NULL) {
  if (anyNA(X)) {
    stop("NA values are not supported: `X` contains NA.", call. = FALSE)
  }
  if (!is.null(x0) && anyNA(x0)) {
    stop("NA values are not supported: `x0` contains NA.", call. = FALSE)
  }
}

.coerce_inv_s <- function(scale_vec, p, eps) {
  inv_s <- as.numeric(scale_vec)
  if (length(inv_s) != p) stop("scale_vec must be length p when rule = 'custom'.", call. = FALSE)
  1 / pmax(abs(inv_s), eps)
}

.compute_inv_scale <- function(X, rule = c("sd", "mad", "iqr", "range", "custom"),
                               scale_vec = NULL, eps = 1e-12) {
  rule <- match.arg(rule)
  if (identical(rule, "custom")) {
    return(.coerce_inv_s(scale_vec, ncol(X), eps))
  }

  s <- switch(rule,
    sd = .col_sd(X),
    mad = .col_mad(X),
    iqr = .col_iqr(X),
    range = .col_rng(X)
  )
  if (length(s) != ncol(X)) stop("Scale vector length must match ncol(X).")
  s <- pmax(s, eps)
  1 / s
}

# Lightweight, readable column scaling (yes, it copies; that’s OK for simplicity)
.scale_cols <- function(M, inv_s) {
  sweep(M, 2L, as.numeric(inv_s), `*`)
}
