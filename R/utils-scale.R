# R/utils-scale.R
# -----------------------------------------------------------------------------

.col_sd <- function(X) apply(X, 2L, stats::sd)
.col_mad <- function(X) apply(X, 2L, stats::mad, constant = 1.4826)
.col_iqr <- function(X) apply(X, 2L, stats::IQR)
.col_rng <- function(X) apply(X, 2L, function(z) diff(range(z)))

.assert_no_na <- function(X, x0 = NULL) {
  if (anyNA(X) || (!is.null(x0) && anyNA(x0))) {
    stop("NA values are not supported (yet).", call. = FALSE)
  }
}

.compute_inv_scale <- function(X, rule = c("sd", "mad", "iqr", "range", "custom"),
                               scale_vec = NULL, eps = 1e-12) {
  rule <- match.arg(rule)
  s <- switch(rule,
    sd = .col_sd(X),
    mad = .col_mad(X),
    iqr = .col_iqr(X),
    range = .col_rng(X),
    custom = {
      if (is.null(scale_vec)) stop("scale_vec must be supplied for rule='custom'.")
      as.numeric(scale_vec)
    }
  )
  if (length(s) != ncol(X)) stop("Scale vector length must match ncol(X).")
  s <- pmax(s, eps)
  1 / s
}

# Lightweight, readable column scaling (yes, it copies; that’s OK for simplicity)
.scale_cols <- function(M, inv_s) {
  sweep(M, 2L, as.numeric(inv_s), `*`)
}
