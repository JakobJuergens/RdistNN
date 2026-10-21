# Shared deterministic NP/TDNN wrapper fixtures.

.make_np_data <- function(n = 40L, p = 2L, m = 4L, seed = 123L) {
  .test_rng_scope({
    X <- matrix(rnorm(n * p), n, p)
    X0 <- matrix(rnorm(m * p), m, p)
    y <- X[, 1] - 0.5 * X[, min(2L, p)] + rnorm(n, sd = 0.25)
    list(X = X, X0 = X0, y = y)
  }, seed = seed)
}

.make_tdnn_data <- function(n = 50L, p = 2L, m = 4L, seed = 123L) {
  .test_rng_scope({
    X  <- matrix(rnorm(n * p), n, p)
    X0 <- matrix(rnorm(m * p), m, p)
    y  <- X[, 1L] - 0.5 * X[, min(2L, p)] + rnorm(n, sd = 0.25)
    list(X = X, X0 = X0, y = y)
  }, seed = seed)
}

.make_np_iface_data <- function(n = 50L, p = 2L, m = 4L, seed = 42L) {
  .test_rng_scope({
    X  <- matrix(rnorm(n * p), n, p)
    X0 <- matrix(rnorm(m * p), m, p)
    y  <- X[, 1L] - 0.5 * X[, min(2L, p)] + rnorm(n, sd = 0.25)
    x_cols <- paste0("X", seq_len(p))
    df <- as.data.frame(X)
    names(df) <- x_cols
    df$Y <- y
    list(X = X, X0 = X0, y = y, df = df, x_cols = x_cols)
  }, seed = seed)
}

.make_tdnn_iface_data <- function(n = 55L, p = 2L, m = 4L, seed = 42L) {
  .test_rng_scope({
    X  <- matrix(rnorm(n * p), n, p)
    X0 <- matrix(rnorm(m * p), m, p)
    y  <- X[, 1L] - 0.5 * X[, min(2L, p)] + rnorm(n, sd = 0.25)
    x_cols <- paste0("X", seq_len(p))
    df <- as.data.frame(X)
    names(df) <- x_cols
    df$Y <- y
    list(X = X, X0 = X0, y = y, df = df, x_cols = x_cols)
  }, seed = seed)
}
