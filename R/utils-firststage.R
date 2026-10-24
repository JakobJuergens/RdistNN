.validate_propensity_clip <- function(clip) {
  if (!is.numeric(clip) || length(clip) != 1L || is.na(clip) ||
      !is.finite(clip) || clip <= 0 || clip >= 0.5 || 1 - clip >= 1) {
    .abort_bad_arg("clip", "must be a finite numeric scalar in (0, 0.5) with 1 - clip < 1")
  }
  as.numeric(clip)
}

.cate_validate_and_clip_propensity <- function(ehat, n, clip = 1e-8, arg = "ehat") {
  clip <- .validate_propensity_clip(clip)
  if (!is.numeric(ehat)) {
    .abort_bad_arg(arg, "must be numeric")
  }
  if (length(ehat) != n) {
    .dnn_abort("argument_invalid",
      sprintf("`%s` must have length equal to `y` (%d); got %d.", arg, n, length(ehat)),
      argument = arg
    )
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

# AIPW / doubly-robust pseudo-outcome:
#   tau_i = m1hat(X_i) - m0hat(X_i)
#           + D_i * (Y_i - m1hat(X_i))/ehat_i
#           - (1 - D_i) * (Y_i - m0hat(X_i))/(1 - ehat_i)
cate_rho_aipw <- function(y, d, m1hat, m0hat, ehat, clip = 1e-8) {
  n <- length(y)
  if (length(d) != n)
    .dnn_abort("argument_invalid",
      sprintf("`d` must have length equal to `y` (%d); got %d.", n, length(d)),
      argument = "d"
    )
  if (length(m1hat) != n)
    .dnn_abort("argument_invalid",
      sprintf("`m1hat` must have length equal to `y` (%d); got %d.", n, length(m1hat)),
      argument = "m1hat"
    )
  if (length(m0hat) != n)
    .dnn_abort("argument_invalid",
      sprintf("`m0hat` must have length equal to `y` (%d); got %d.", n, length(m0hat)),
      argument = "m0hat"
    )
  e_used <- .cate_validate_and_clip_propensity(ehat, n = n, clip = clip)
  m1hat - m0hat + d * (y - m1hat) / e_used -
    (1 - d) * (y - m0hat) / (1 - e_used)
}
