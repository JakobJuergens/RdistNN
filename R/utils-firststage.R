.cate_validate_and_clip_propensity <- function(ehat, n, clip = 1e-8, arg = "ehat") {
  if (!is.numeric(clip) || length(clip) != 1L || is.na(clip) ||
      !is.finite(clip) || clip <= 0 || clip >= 0.5) {
    .abort_bad_arg("clip", "must be a finite numeric scalar in (0, 0.5)")
  }
  if (!is.numeric(ehat)) {
    .abort_bad_arg(arg, "must be numeric")
  }
  if (length(ehat) != n) {
    stop(sprintf("`%s` must have length equal to `y` (%d); got %d.", arg, n, length(ehat)), call. = FALSE)
  }
  if (.has_na_only(ehat)) {
    .abort_bad_arg(arg, "must not contain missing values")
  }
  if (.has_nonfinite_non_na(ehat)) {
    .abort_bad_arg(arg, "must contain only finite values")
  }
  if (any(ehat < 0 | ehat > 1)) {
    .abort_bad_arg(arg, "must lie in [0, 1] before clipping")
  }

  pmin(pmax(as.numeric(ehat), as.numeric(clip)), 1 - as.numeric(clip))
}

# R-learner style pseudo-outcome:
#   rho_i = (D_i - ehat_i) * (Y_i - mhat_i) / (ehat_i * (1 - ehat_i))
# ehat, mhat come from YOUR code; no fitting happens here.
cate_rho_rlearner <- function(y, d, ehat, mhat, clip = 1e-8) {
  n <- length(y)
  if (length(d) != n)
    stop(sprintf("`d` must have length equal to `y` (%d); got %d.", n, length(d)), call. = FALSE)
  if (length(mhat) != n)
    stop(sprintf("`mhat` must have length equal to `y` (%d); got %d.", n, length(mhat)), call. = FALSE)
  e_used <- .cate_validate_and_clip_propensity(ehat, n = n, clip = clip)
  (d - e_used) * (y - mhat) / (e_used * (1 - e_used))
}

# AIPW / doubly-robust pseudo-outcome (optional alternative):
#   tau_i = m1hat(X_i) - m0hat(X_i)
#           + D_i * (Y_i - m1hat(X_i))/ehat_i
#           - (1 - D_i) * (Y_i - m0hat(X_i))/(1 - ehat_i)
cate_rho_aipw <- function(y, d, m1hat, m0hat, ehat, clip = 1e-8) {
  n <- length(y)
  if (length(d) != n)
    stop(sprintf("`d` must have length equal to `y` (%d); got %d.", n, length(d)), call. = FALSE)
  if (length(m1hat) != n)
    stop(sprintf("`m1hat` must have length equal to `y` (%d); got %d.", n, length(m1hat)), call. = FALSE)
  if (length(m0hat) != n)
    stop(sprintf("`m0hat` must have length equal to `y` (%d); got %d.", n, length(m0hat)), call. = FALSE)
  e_used <- .cate_validate_and_clip_propensity(ehat, n = n, clip = clip)
  m1hat - m0hat + d * (y - m1hat) / e_used -
    (1 - d) * (y - m0hat) / (1 - e_used)
}
