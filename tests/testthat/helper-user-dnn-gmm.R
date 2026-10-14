# Shared GMM model fixtures for estimator, interface, validation, and method tests.
.test_dr_iv_model <- function(q) {
  gmm_model_dr_iv(
    z_cols = paste0("Z", seq_len(q)),
    pi_cols = paste0("piZ", seq_len(q), "_oracle"),
    muY_col = "muY_oracle",
    muD_col = "muD_oracle"
  )
}

.test_mean2_gmm_model <- function() {
  new_gmm_model(
    name = "mean2",
    setup = function(sample, n, ...) {
      y <- as.matrix(sample[, c("Y1", "Y2"), drop = FALSE])
      list(y = y)
    },
    moment_obs = function(sample, theta, ctx) {
      sweep(ctx$y, 2L, as.numeric(theta), FUN = "-")
    },
    jacobian_obs = function(sample, theta, ctx) {
      n <- nrow(ctx$y)
      out <- array(0, dim = c(n, 2L, 2L))
      out[, 1L, 1L] <- -1
      out[, 2L, 2L] <- -1
      out
    },
    start = function(sample, ctx, w, x0) {
      as.numeric(crossprod(ctx$y, w))
    },
    coef_names = c("theta1", "theta2"),
    moment_names = c("mean1", "mean2")
  )
}
