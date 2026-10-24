# sim-iv-data.R
# Data-generating processes for the local IV-GMM simulation study.
#
# Public API:
#   sim_iv_data() -- draw one sample from a named IV DGP
#
# The DGPs are intentionally estimator-free. They provide ordinary (Y, D, Z, X)
# data plus oracle nuisance columns for the doubly-robust IV moment
#
#   (Z_k - pi_k(X)) * (Y - muY(X) - (D - muD(X)) * theta).

# -- internal helpers ---------------------------------------------------------

.sim_iv_m0 <- function(X) {
  1 + X[, 1L] - 0.5 * X[, 2L]
}

.sim_iv_default_x0 <- function(design, kx) {
  switch(design,
    dgp1_iv_justid_linear = rep(0, kx),
    dgp2_iv_homog_overid = c(1, rep(0, max(0L, kx - 1L))),
    dgp3_iv_heterog_overid = c(1, rep(0, max(0L, kx - 1L))),
    dgp4_iv_weak_local = c(1, rep(0, max(0L, kx - 1L))),
    dgp5_iv_correlated_instruments = c(0.75, rep(0, max(0L, kx - 1L))),
    dgp6_iv_nonlinear_firststage = c(1, rep(0, max(0L, kx - 1L))),
    rep(0, kx)
  )
}

.sim_iv_tau_linear <- function(xx) {
  xx <- as.matrix(xx)
  beta <- c(0.5, -0.25, rep(0, max(0L, ncol(xx) - 2L)))
  as.numeric(xx %*% beta[seq_len(ncol(xx))])
}

.sim_iv_tau_nonlinear <- function(xx) {
  xx <- as.matrix(xx)
  0.5 + 0.75 * sin(1.25 * xx[, 1L])
}

.sim_iv_tau_curved <- function(xx) {
  xx <- as.matrix(xx)
  sin(1.5 * xx[, 1L]) + 0.4 * xx[, 1L]^2
}

.sim_iv_binary_pi <- function(X, q, clip) {
  X <- as.matrix(X)
  out <- matrix(NA_real_, nrow(X), q)
  for (j in seq_len(q)) {
    x2 <- if (ncol(X) >= 2L) X[, 2L] else 0
    lin <- -0.15 + 0.35 * X[, 1L] - 0.25 * x2 + 0.20 * (j - 1L)
    out[, j] <- pmin(pmax(stats::plogis(lin), clip), 1 - clip)
  }
  out
}

.sim_iv_continuous_pi <- function(X, q) {
  X <- as.matrix(X)
  out <- matrix(NA_real_, nrow(X), q)
  for (j in seq_len(q)) {
    x2 <- if (ncol(X) >= 2L) X[, 2L] else 0
    out[, j] <- 0.25 * X[, 1L] - 0.15 * x2 + 0.10 * (j - 1L)
  }
  out
}

.sim_iv_mixed_pi <- function(X, q, clip) {
  pi_mat <- .sim_iv_continuous_pi(X, q)
  pi_mat[, 1L] <- .sim_iv_binary_pi(X, 1L, clip)[, 1L]
  pi_mat
}

.sim_iv_d_base <- function(X) {
  X <- as.matrix(X)
  x2 <- if (ncol(X) >= 2L) X[, 2L] else 0
  -0.15 + 0.30 * X[, 1L] - 0.20 * x2
}

.sim_iv_gamma_const <- function(X, q, strength = 0.75) {
  X <- as.matrix(X)
  base <- strength * exp(-0.15 * (seq_len(q) - 1L))
  matrix(base, nrow(X), q, byrow = TRUE)
}

.sim_iv_gamma_weak <- function(X, q, x0, local_strength = 0.08,
                               far_strength = 0.85, radius = 0.6) {
  X <- as.matrix(X)
  d2 <- rowSums((X - matrix(x0, nrow(X), ncol(X), byrow = TRUE))^2)
  scale <- local_strength + (far_strength - local_strength) *
    (1 - exp(-d2 / (2 * radius^2)))
  base <- exp(-0.15 * (seq_len(q) - 1L))
  scale * matrix(base, nrow(X), q, byrow = TRUE)
}

.sim_iv_gamma_nonlinear <- function(X, q) {
  X <- as.matrix(X)
  x2 <- if (ncol(X) >= 2L) X[, 2L] else 0
  base <- 0.55 + 0.25 * stats::plogis(X[, 1L] - 0.5 * x2)
  base * matrix(exp(-0.20 * (seq_len(q) - 1L)), nrow(X), q, byrow = TRUE)
}

.sim_iv_corr_matrix <- function(q, rho = 0.85) {
  mat <- matrix(rho, q, q)
  diag(mat) <- 1
  mat
}

.sim_iv_draw_z <- function(pi_mat, instrument_type, corr = NULL) {
  n <- nrow(pi_mat)
  q <- ncol(pi_mat)
  Z <- matrix(NA_real_, n, q)

  bin_cols <- instrument_type == "binary"
  cont_cols <- instrument_type == "continuous"

  if (any(bin_cols)) {
    for (j in which(bin_cols)) {
      Z[, j] <- stats::rbinom(n, size = 1L, prob = pi_mat[, j])
    }
  }

  if (any(cont_cols)) {
    cont_idx <- which(cont_cols)
    n_cont <- length(cont_idx)
    if (!is.null(corr) && n_cont > 1L) {
      R <- chol(corr)
      eps <- matrix(stats::rnorm(n * n_cont), n, n_cont) %*% R
    } else {
      eps <- matrix(stats::rnorm(n * n_cont), n, n_cont)
    }
    Z[, cont_idx] <- pi_mat[, cont_idx, drop = FALSE] + eps
  }

  colnames(Z) <- paste0("Z", seq_len(q))
  Z
}

.sim_iv_muD <- function(pi_mat, gamma_mat, d_base, instrument_type, corr = NULL) {
  q <- ncol(pi_mat)
  n <- nrow(pi_mat)
  muD <- numeric(n)

  bin_idx <- which(instrument_type == "binary")
  cont_idx <- which(instrument_type == "continuous")

  if (length(bin_idx) == 0L) {
    Sigma <- if (is.null(corr)) diag(q) else corr
    var_term <- rowSums((gamma_mat %*% Sigma) * gamma_mat)
    return(stats::pnorm(d_base / sqrt(1 + var_term)))
  }

  grids <- expand.grid(rep(list(c(0, 1)), length(bin_idx)))
  grids <- as.matrix(grids)

  cont_var <- rep(0, n)
  if (length(cont_idx) > 0L) {
    Sigma_cont <- if (is.null(corr)) diag(length(cont_idx)) else corr[cont_idx, cont_idx, drop = FALSE]
    G_cont <- gamma_mat[, cont_idx, drop = FALSE]
    cont_var <- rowSums((G_cont %*% Sigma_cont) * G_cont)
  }

  for (row in seq_len(nrow(grids))) {
    z_bin <- grids[row, ]
    pi_bin <- pi_mat[, bin_idx, drop = FALSE]
    prob <- rep(1, n)
    resid_shift <- rep(0, n)

    for (h in seq_along(bin_idx)) {
      pj <- pi_bin[, h]
      zj <- z_bin[h]
      prob <- prob * if (zj == 1) pj else (1 - pj)
      resid_shift <- resid_shift + gamma_mat[, bin_idx[h]] * (zj - pj)
    }

    muD <- muD + prob * stats::pnorm((d_base + resid_shift) / sqrt(1 + cont_var))
  }

  muD
}

.sim_iv_assemble <- function(X, Z, D, Y, pi_mat, muD, theta_true,
                             theta_fn, pi_fn, muD_fn, muY_fn, x0,
                             instrument_type) {
  p <- ncol(X)
  q <- ncol(Z)
  x_names <- paste0("X", seq_len(p))
  z_names <- paste0("Z", seq_len(q))
  pi_names <- paste0("piZ", seq_len(q), "_oracle")

  df <- data.frame(X, Y = as.numeric(Y), D = as.integer(D), Z,
                   pi_mat,
                   muY_oracle = as.numeric(muY_fn(X)),
                   muD_oracle = as.numeric(muD),
                   theta_oracle = as.numeric(theta_true),
                   check.names = FALSE)
  names(df)[seq_len(p)] <- x_names
  names(df)[p + 2L + seq_len(q)] <- z_names
  names(df)[p + 2L + q + seq_len(q)] <- pi_names

  theta_x0 <- as.numeric(theta_fn(matrix(x0, nrow = 1L)))

  list(
    data = df,
    x0 = as.numeric(x0),
    q = as.integer(q),
    theta_x0 = theta_x0,
    truth = list(
      theta_fn = theta_fn,
      pi_fn = pi_fn,
      muD_fn = muD_fn,
      muY_fn = muY_fn,
      instrument_type = instrument_type
    )
  )
}

.sim_iv_draw_common <- function(n, kx, q, sigma, clip, seed, x0,
                                theta_fn, pi_fn, gamma_fn,
                                instrument_type, corr = NULL,
                                direct_fn = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- matrix(stats::rnorm(n * kx), n, kx)
  pi_mat <- pi_fn(X)
  gamma_mat <- gamma_fn(X)
  d_base <- .sim_iv_d_base(X)

  Z <- .sim_iv_draw_z(pi_mat, instrument_type, corr = corr)
  Z_resid <- Z - pi_mat
  probD <- stats::pnorm(d_base + rowSums(Z_resid * gamma_mat))
  D <- stats::rbinom(n, size = 1L, prob = probD)

  muD_fn <- function(xx) {
    xx <- as.matrix(xx)
    pi_x <- pi_fn(xx)
    gamma_x <- gamma_fn(xx)
    .sim_iv_muD(pi_x, gamma_x, .sim_iv_d_base(xx), instrument_type, corr = corr)
  }

  muD <- muD_fn(X)
  theta_true <- theta_fn(X)
  direct <- if (is.null(direct_fn)) rep(0, n) else direct_fn(X, Z, pi_mat)
  eps <- 0.7 * (D - probD) + stats::rnorm(n, sd = sigma)
  Y <- .sim_iv_m0(X) + theta_true * D + direct + eps

  muY_fn <- function(xx) {
    xx <- as.matrix(xx)
    .sim_iv_m0(xx) + theta_fn(xx) * muD_fn(xx)
  }

  .sim_iv_assemble(
    X = X, Z = Z, D = D, Y = Y, pi_mat = pi_mat, muD = muD,
    theta_true = theta_true, theta_fn = theta_fn, pi_fn = pi_fn,
    muD_fn = muD_fn, muY_fn = muY_fn, x0 = x0,
    instrument_type = instrument_type
  )
}

# -- DGP catalog --------------------------------------------------------------

.sim_iv_dgp1 <- function(n, kx, sigma, clip, seed, x0) {
  q <- 1L
  pi_fn <- function(xx) .sim_iv_binary_pi(xx, q = q, clip = clip)
  gamma_fn <- function(xx) .sim_iv_gamma_const(xx, q = q, strength = 0.85)
  .sim_iv_draw_common(
    n, kx, q, sigma, clip, seed, x0,
    theta_fn = .sim_iv_tau_linear,
    pi_fn = pi_fn,
    gamma_fn = gamma_fn,
    instrument_type = "binary"
  )
}

.sim_iv_dgp2 <- function(n, kx, q, sigma, clip, seed, x0) {
  pi_fn <- function(xx) .sim_iv_binary_pi(xx, q = q, clip = clip)
  gamma_fn <- function(xx) .sim_iv_gamma_const(xx, q = q, strength = 0.80)
  .sim_iv_draw_common(
    n, kx, q, sigma, clip, seed, x0,
    theta_fn = .sim_iv_tau_nonlinear,
    pi_fn = pi_fn,
    gamma_fn = gamma_fn,
    instrument_type = rep("binary", q)
  )
}

.sim_iv_dgp3 <- function(n, kx, q, sigma, clip, seed, x0) {
  pi_fn <- function(xx) .sim_iv_binary_pi(xx, q = q, clip = clip)
  gamma_fn <- function(xx) .sim_iv_gamma_const(xx, q = q, strength = 0.75)
  direct_fn <- function(X, Z, pi_mat) {
    local <- exp(-rowSums((X - matrix(x0, nrow(X), ncol(X), byrow = TRUE))^2) / (2 * 0.75^2))
    delta <- seq(-0.45, 0.45, length.out = ncol(Z))
    as.numeric(local * rowSums((Z - pi_mat) * matrix(delta, nrow(X), ncol(Z), byrow = TRUE)))
  }
  .sim_iv_draw_common(
    n, kx, q, sigma, clip, seed, x0,
    theta_fn = .sim_iv_tau_nonlinear,
    pi_fn = pi_fn,
    gamma_fn = gamma_fn,
    instrument_type = rep("binary", q),
    direct_fn = direct_fn
  )
}

.sim_iv_dgp4 <- function(n, kx, q, sigma, clip, seed, x0) {
  pi_fn <- function(xx) .sim_iv_binary_pi(xx, q = q, clip = clip)
  gamma_fn <- function(xx) .sim_iv_gamma_weak(xx, q = q, x0 = x0)
  .sim_iv_draw_common(
    n, kx, q, sigma, clip, seed, x0,
    theta_fn = .sim_iv_tau_nonlinear,
    pi_fn = pi_fn,
    gamma_fn = gamma_fn,
    instrument_type = rep("binary", q)
  )
}

.sim_iv_dgp5 <- function(n, kx, q, sigma, clip, seed, x0, rho = 0.85) {
  pi_fn <- function(xx) .sim_iv_continuous_pi(xx, q = q)
  gamma_fn <- function(xx) .sim_iv_gamma_const(xx, q = q, strength = 0.70)
  corr <- .sim_iv_corr_matrix(q, rho = rho)
  .sim_iv_draw_common(
    n, kx, q, sigma, clip, seed, x0,
    theta_fn = .sim_iv_tau_curved,
    pi_fn = pi_fn,
    gamma_fn = gamma_fn,
    instrument_type = rep("continuous", q),
    corr = corr
  )
}

.sim_iv_dgp6 <- function(n, kx, q, sigma, clip, seed, x0) {
  pi_fn <- function(xx) .sim_iv_mixed_pi(xx, q = q, clip = clip)
  gamma_fn <- function(xx) .sim_iv_gamma_nonlinear(xx, q = q)
  instrument_type <- c("binary", rep("continuous", max(0L, q - 1L)))
  .sim_iv_draw_common(
    n, kx, q, sigma, clip, seed, x0,
    theta_fn = .sim_iv_tau_nonlinear,
    pi_fn = pi_fn,
    gamma_fn = gamma_fn,
    instrument_type = instrument_type
  )
}

#' Draw one sample from a named IV-GMM simulation DGP
#'
#' Generates a synthetic local IV-GMM dataset with oracle nuisance quantities for
#' the doubly-robust IV moment
#' \deqn{
#'   (Z_k - \pi_k(X))\{Y - \mu_Y(X) - (D - \mu_D(X))\theta\}.
#' }
#' Use the returned data with [gmm_dnn()], [gmm_model_dr_iv()] or the
#' simulation methods listed by [sim_iv_methods()]. No model is fit here.
#'
#' @param design Character DGP name. One of `"dgp1_iv_justid_linear"`,
#'   `"dgp2_iv_homog_overid"`, `"dgp3_iv_heterog_overid"`,
#'   `"dgp4_iv_weak_local"`, `"dgp5_iv_correlated_instruments"`, or
#'   `"dgp6_iv_nonlinear_firststage"`.
#' @param n Integer sample size.
#' @param kx Integer covariate dimension, at least 2.
#' @param q Integer number of instruments. The simulation DGPs support
#'   `1 <= q <= 3`. The just-identified design always returns `q = 1`.
#' @param sigma Noise standard deviation in `Y`.
#' @param clip Clipping bound in `(0, 0.5)` for binary instrument probabilities.
#' @param seed Optional integer seed, used locally without changing the caller's
#'   RNG state. `NULL` uses the current stream.
#' @param x0 Optional numeric length-`kx` evaluation point. When omitted, a
#'   design-specific interior point is used.
#' @param ... Extra arguments forwarded to the underlying DGP. Currently `rho`
#'   is used by `"dgp5_iv_correlated_instruments"`.
#'
#' @return A list with:
#' \itemize{
#'   \item `data`: a `data.frame` with `X1..Xkx`, `Y`, `D`, `Z1..Zq`, oracle
#'     nuisance columns `piZ*_oracle`, `muY_oracle`, `muD_oracle`, and the
#'     target column `theta_oracle`.
#'   \item `x0`: numeric length-`kx` evaluation point.
#'   \item `q`: realized number of instruments.
#'   \item `theta_x0`: true local IV target at `x0`.
#'   \item `truth`: named list of oracle functions.
#' }
#'
#' @examples
#' dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 200, kx = 2, q = 2, seed = 1)
#' head(dgp$data)
#' dgp$theta_x0
#'
#' @export
sim_iv_data <- function(design, n, kx = 2L, q = 2L, sigma = 0.5,
                        clip = 1e-3, seed = NULL, x0 = NULL, ...) {
  n <- .validate_scalar_int(n, arg = "n", min = 1L)
  kx <- .validate_scalar_int(kx, arg = "kx", min = 2L)
  q <- .validate_scalar_int(q, arg = "q", min = 1L)
  if (q > 3L) {
    .abort_bad_arg("q", "must be <= 3 for the Phase-A IV simulation DGPs")
  }
  if (!is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma) || sigma < 0) {
    .abort_bad_arg("sigma", "must be a non-negative finite numeric scalar")
  }
  sigma <- as.numeric(sigma)
  if (!is.numeric(clip) || length(clip) != 1L || !is.finite(clip) || clip <= 0 || clip >= 0.5) {
    .abort_bad_arg("clip", "must be a finite numeric scalar in (0, 0.5)")
  }
  clip <- as.numeric(clip)

  valid <- c(
    "dgp1_iv_justid_linear",
    "dgp2_iv_homog_overid",
    "dgp3_iv_heterog_overid",
    "dgp4_iv_weak_local",
    "dgp5_iv_correlated_instruments",
    "dgp6_iv_nonlinear_firststage"
  )
  if (!is.character(design) || length(design) != 1L || !design %in% valid) {
    .dnn_abort("argument_invalid", paste0(
      "Unknown IV design: '", design, "'. Choose from: ",
      paste(valid, collapse = ", "), "."), argument = "design")
  }
  if (is.null(x0)) {
    x0 <- .sim_iv_default_x0(design, kx)
  } else {
    x0 <- .validate_target_vector(x0, p = kx, x0_arg = "x0", x_arg = "X", finite = TRUE)
  }

  if (!is.null(seed)) {
    restore_rng <- .rdistnn_rng_restore()
    on.exit(restore_rng(), add = TRUE)
  }
  switch(design,
    dgp1_iv_justid_linear = .sim_iv_dgp1(n, kx, sigma, clip, seed, x0, ...),
    dgp2_iv_homog_overid = .sim_iv_dgp2(n, kx, q, sigma, clip, seed, x0, ...),
    dgp3_iv_heterog_overid = .sim_iv_dgp3(n, kx, q, sigma, clip, seed, x0, ...),
    dgp4_iv_weak_local = .sim_iv_dgp4(n, kx, q, sigma, clip, seed, x0, ...),
    dgp5_iv_correlated_instruments = .sim_iv_dgp5(n, kx, q, sigma, clip, seed, x0, ...),
    dgp6_iv_nonlinear_firststage = .sim_iv_dgp6(n, kx, q, sigma, clip, seed, x0, ...)
  )
}
