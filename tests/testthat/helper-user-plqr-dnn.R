# Shared helpers for plqr_dnn() user-level tests.

.make_plqr_data <- function(n = 80L, p = 2L, p_d = 2L, m = 4L, tau = 0.5, seed = 1L) {
  .test_rng_scope({
    demo <- plqr_demo_data(n = n, p = p, p_d = p_d, tau = tau, seed = 1)
    df <- demo$data
    X <- as.matrix(df[paste0("X", seq_len(p))])
    D <- as.matrix(df[paste0("D", seq_len(p_d))])
    mD <- as.matrix(df[paste0("mD", seq_len(p_d))])
    g <- as.numeric(df$g_oracle)
    Y <- as.numeric(df$Y)
    X0 <- matrix(rnorm(m * p), m, p)
    list(
      df = df, X = X, D = D, mD = mD, g = g, Y = Y, X0 = X0,
      theta_eval = demo$truth$theta_true, tau = tau,
      x_cols = paste0("X", seq_len(p)),
      D_cols = paste0("D", seq_len(p_d)),
      mD_cols = paste0("mD", seq_len(p_d))
    )
  }, seed = seed)
}

.make_plqr_inference_data <- function(n = 90L, p = 2L, p_d = 2L, m = 2L, tau = 0.5, seed = 818L) {
  .test_rng_scope({
    demo <- plqr_demo_data(n = n, p = p, p_d = p_d, tau = tau, seed = 818)
    df <- demo$data
    X <- as.matrix(df[paste0("X", seq_len(p))])
    D <- as.matrix(df[paste0("D", seq_len(p_d))])
    mD <- as.matrix(df[paste0("mD", seq_len(p_d))])
    g <- as.numeric(df$g_oracle)
    Y <- as.numeric(df$Y)
    X0 <- matrix(rnorm(m * p), m, p)
    density <- 0.75 + 0.1 * seq_len(n) / n
    df$density_oracle <- density
    list(
      df = df, X = X, D = D, mD = mD, g = g, Y = Y, X0 = X0,
      density = density, theta_eval = demo$truth$theta_true, tau = tau,
      x_cols = paste0("X", seq_len(p)),
      D_cols = paste0("D", seq_len(p_d)),
      mD_cols = paste0("mD", seq_len(p_d))
    )
  }, seed = seed)
}

.manual_plqr_sandwich <- function(d, w) {
  d_resid <- d$D - d$mD
  S <- d$tau * (1 - d$tau) * crossprod(d_resid, d_resid * w)
  J <- crossprod(d_resid, d$D * as.numeric(d$density * w))
  J_inv <- solve(J)
  omega <- J_inv %*% S %*% t(J_inv)
  vcov <- omega / (1 / sum(w^2))
  list(S = S, J = J, omega = omega, vcov = vcov)
}
