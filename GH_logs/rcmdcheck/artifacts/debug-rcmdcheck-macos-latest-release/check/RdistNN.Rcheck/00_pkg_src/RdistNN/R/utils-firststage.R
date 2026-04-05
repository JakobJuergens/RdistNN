# R-learner style pseudo-outcome:
#   rho_i = (D_i - ehat_i) * (Y_i - mhat_i) / (ehat_i * (1 - ehat_i))
# ehat, mhat come from YOUR code; no fitting happens here.
cate_rho_rlearner <- function(y, d, ehat, mhat, clip = NULL) {
  stopifnot(length(y) == length(d), length(ehat) == length(y), length(mhat) == length(y))
  if (!is.null(clip)) {
    ehat <- pmin(pmax(ehat, clip), 1 - clip)
  }
  (d - ehat) * (y - mhat) / (ehat * (1 - ehat))
}

# AIPW / doubly-robust pseudo-outcome (optional alternative):
#   tau_i = m1hat(X_i) - m0hat(X_i)
#           + D_i * (Y_i - m1hat(X_i))/ehat_i
#           - (1 - D_i) * (Y_i - m0hat(X_i))/(1 - ehat_i)
cate_rho_aipw <- function(y, d, m1hat, m0hat, ehat, clip = NULL) {
  stopifnot(length(y) == length(d), length(m1hat) == length(y), length(m0hat) == length(y), length(ehat) == length(y))
  if (!is.null(clip)) {
    ehat <- pmin(pmax(ehat, clip), 1 - clip)
  }
  m1hat - m0hat + d * (y - m1hat) / ehat - (1 - d) * (y - m0hat) / (1 - ehat)
}
