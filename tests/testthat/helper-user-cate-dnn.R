# Shared helpers for cate_dnn user-level tests.

.make_data <- function(n = 40L, p = 3L, m = 4L) {
  X <- matrix(rnorm(n * p), n, p)
  X0 <- matrix(rnorm(m * p), m, p)
  beta <- seq_len(p) / p
  Y <- drop(X %*% beta) + rnorm(n, sd = 0.2)
  list(X = X, X0 = X0, Y = Y, beta = beta)
}

.make_df_data <- function(n = 60L, p = 2L, m = 5L, seed = 7L) {
  set.seed(seed)
  demo <- cate_demo_data(n = n, p = p, design = "hetero", sigma = 0.3, seed = seed)
  df <- demo$data
  df$m1_oracle <- df$m0_oracle + df$tau_oracle
  x_cols <- paste0("X", seq_len(p))
  X0 <- matrix(rnorm(m * p), m, p)

  Y <- df$Y
  T_ <- df$D
  m0 <- df$m0_oracle
  m1 <- df$m1_oracle
  e <- df$e_oracle
  m_oracle <- df$m_oracle
  den <- pmax(e * (1 - e), 1e-8)
  psi_aipw <- (T_ * (Y - m1) - (1 - T_) * (Y - m0)) / den + (m1 - m0)
  psi_robinson <- (Y - m_oracle) * (T_ - e) / den

  list(
    df = df,
    x_cols = x_cols,
    X0 = X0,
    psi_aipw = psi_aipw,
    psi_robinson = psi_robinson
  )
}
