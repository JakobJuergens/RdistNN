# PLQR demo data (oracle nuisances; no internal fitting)
# + test-friendly injection hooks
# + multi-scenario generator for unit tests
# -------------------------------------------------------------------

# --- Internal: Asymmetric Laplace draws (location 0, asymmetry = tau) ----
# Ensures the tau-quantile of the returned draws is zero, with a density
# that is continuous at the quantile (unlike the historical two-piece
# generator, whose density jumped at zero for tau != 0.5). The scale is
# parameterized so that tau = 0.5 draws match the historical generator
# exactly, giving the closed-form threshold density
#   f(0) = 4 * tau * (1 - tau) / sigma
# returned by .plqr_ald_density0().
.r_ald <- function(n, tau = 0.5, sigma = 1) {
  n <- .validate_scalar_int(n, arg = "n", min = 1L)
  tau <- .validate_scalar_prob(tau, arg = "tau")
  if (!is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
    .abort_bad_arg("sigma", "must be a positive finite numeric scalar")
  }
  sigma <- as.numeric(sigma)
  scale <- sigma / 4
  u <- runif(n)
  ifelse(
    u < tau,
    (scale / (1 - tau)) * log(u / tau),
    -(scale / tau) * log((1 - u) / (1 - tau))
  )
}

# Density of the .r_ald() error at its tau-quantile (zero).
.plqr_ald_density0 <- function(tau, sigma) {
  4 * tau * (1 - tau) / sigma
}

# Full density function of the .r_ald() error, vectorized in v.
.plqr_ald_density_fn <- function(tau, sigma) {
  f0 <- .plqr_ald_density0(tau, sigma)
  scale <- sigma / 4
  function(v) {
    v <- as.numeric(v)
    ifelse(
      v <= 0,
      f0 * exp((1 - tau) * v / scale),
      f0 * exp(-tau * v / scale)
    )
  }
}

#' Minimal PLQR demo data (oracle nuisances; no internal fitting) (experimental)
#'
#' Generates data for a partially linear quantile regression (PLQR) model
#' \deqn{
#'   Q_{\tau}(Y \mid D, U) \;=\; D^{\top}\,\theta_{0} \;+\; g_{0}(U) \, ,
#' }
#' with a user-chosen design for the nonparametric shape nuisance
#' \eqn{g_{0}}; the demo specializes to \eqn{U = X}. Errors follow an
#' Asymmetric Laplace distribution (ALD) at level \eqn{\tau}, so the
#' conditional \eqn{\tau}-quantile equals the right-hand side.
#'
#' The function also returns the oracle
#' \eqn{m_{D}(U) \equiv \mathbb{E}[D \mid U]} to build the orthogonal PLQR score
#' \deqn{
#'   \psi(W;\theta)
#'   \;=\;
#'   (\tau - \mathbf{1}\{\,Y \le D^{\top}\theta + g_{0}(U)\,\})\,
#'   (D - m_{D}(U))\, .
#' }
#'
#' @section Lifecycle:
#' **Experimental.** The PLQR API and inference interfaces may change.
#' @param n Integer; sample size.
#' @param p Positive integer number of features in \eqn{X}.
#' @param p_d Integer; dimension of \eqn{D} (number of partially linear regressors).
#' @param tau Quantile level in \eqn{(0,1)}.
#' @param sigma Positive ALD scale for the regression error (larger means noisier \eqn{Y}).
#' @param design Character; one of \code{"linear"}, \code{"sinusoidal"}, \code{"bump"} for \eqn{g_{0}}.
#' @param seed Optional integer seed, used locally without changing the caller's
#'   RNG state. `NULL` uses the current stream.
#' @param x0 Optional target point in \eqn{\mathbb{R}^{p}} for localization (defaults to the zero vector).
#' @param inject Optional, one of \code{NULL}, \code{"na"}, \code{"inf"}; if provided,
#'   inject a single non-finite into \code{inject_col} at row \code{idx_inject}.
#'   Useful for testing runtime guards.
#' @param inject_col Character; which column to corrupt when \code{inject} is not \code{NULL}.
#'   One of \code{"Y"}, \code{"g_oracle"}, \code{"D1"}, \code{"mD1"}, \code{"rD1"}.
#' @param idx_inject Integer row index to corrupt (default 3).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{data}: a \code{data.frame} with columns \code{X1..Xp}, \code{D1..Dp_d}, \code{Y};
#'         oracle nuisances \code{g_oracle} (the value of \eqn{g_{0}(X)}),
#'         \code{mD1..mDp_d} (components of \eqn{m_{D}(X)}),
#'         residuals \code{rD1..rDp_d} (components of \eqn{D - m_{D}(X)}),
#'         and \code{density_oracle}, the true conditional density of
#'         \eqn{Y} at the quantile threshold. Because the demo's ALD errors
#'         are independent of \eqn{(X, D)} and homoskedastic, this density
#'         is constant, \eqn{f(0) = 4\tau(1-\tau)/\sigma}.
#'   \item \code{x0}: the target point used.
#'   \item \code{truth}: a list with \code{theta_true} (length \code{p_d}),
#'         \code{g_fn(x)} and \code{mD_fn(x)} closures to evaluate ground truth at new \eqn{X},
#'         and \code{density_fn(v)}, the error density at signed distance
#'         \eqn{v} from the quantile threshold (\code{density_fn(0)} is the
#'         \code{density_oracle} value).
#'   \item \code{g_x0}: the true value \eqn{g_{0}(x_{0})}.
#' }
#'
#' @examples
#' demo <- plqr_demo_data(
#'   n = 300, p = 2, p_d = 1, tau = 0.5,
#'   design = "sinusoidal", seed = 1
#' )
#' head(demo$data)
#' demo$truth$theta_true
#' demo$g_x0
#'
#' @export
plqr_demo_data <- function(n = 300,
                           p = 2,
                           p_d = 1,
                           tau = 0.5,
                           sigma = 0.5,
                           design = c("linear", "sinusoidal", "bump"),
                           seed = NULL,
                           x0 = NULL,
                           inject = NULL,
                           inject_col = c("Y", "g_oracle", "D1", "mD1", "rD1")[1],
                           idx_inject = 3L) {
  n <- .validate_scalar_int(n, arg = "n", min = 1L)
  p <- .validate_scalar_int(p, arg = "p", min = 1L)
  p_d <- .validate_scalar_int(p_d, arg = "p_d", min = 1L)
  tau <- .validate_scalar_prob(tau, arg = "tau")
  if (!is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
    .abort_bad_arg("sigma", "must be a positive finite numeric scalar")
  }
  sigma <- as.numeric(sigma)
  design <- .dnn_match_arg(design)
  if (is.null(x0)) {
    x0 <- rep(0, p)
  } else {
    x0 <- .validate_target_vector(x0, p = p, x0_arg = "x0", x_arg = "X", finite = TRUE)
  }
  if (!is.null(seed)) {
    restore_rng <- .rdistnn_rng_restore()
    on.exit(restore_rng(), add = TRUE)
    set.seed(seed)
  }

  # Covariates
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", seq_len(p))

  # True theta (partially linear effect at quantile tau)
  theta_true <- rep(0.5, p_d)

  # Nonparametric shape nuisance g_0 (evaluated at U = X in the demo)
  g_fn <- switch(design,
    linear = function(xx) {
      xx <- as.matrix(xx)
      0.5 * xx[, 1]
    },
    sinusoidal = function(xx) {
      xx <- as.matrix(xx)
      0.5 * sin(xx[, 1]) + if (ncol(xx) >= 2) 0.25 * xx[, 2]^2 else 0
    },
    bump = {
      center <- c(0.5, -0.5, rep(0, max(0, p - 2L)))[seq_len(p)]
      function(xx) {
        xx <- as.matrix(xx)
        d2 <- rowSums((xx - rep(center, each = nrow(xx)))^2)
        0.2 + 0.8 * exp(-0.5 * d2)
      }
    }
  )
  gX <- g_fn(X)

  # Partially linear regressors D with oracle E(D | X) = mD(X):
  # simple linear mD with different loadings per column of D
  Gamma <- matrix(0, nrow = p, ncol = p_d)
  for (j in seq_len(p_d)) {
    loadings <- c(0.6, -0.4, rep(0, max(0, p - 2L))) + (j - 1) * 0.2
    Gamma[, j] <- loadings[seq_len(p)]
  }
  mD <- X %*% Gamma # n x p_d
  V <- matrix(rnorm(n * p_d, sd = 0.6), n, p_d)
  D <- mD + V
  colnames(D) <- paste0("D", seq_len(p_d))

  # Error with tau-quantile zero
  eps <- .r_ald(n, tau = tau, sigma = sigma)

  # Outcome: Y = D^T theta_true + g_0(X) + eps  with Q_tau(eps | X, D) = 0
  Y <- as.vector(D %*% theta_true) + gX + eps

  # Tidy outputs
  df <- data.frame(X, D, Y = as.numeric(Y), g_oracle = as.numeric(gX))
  for (j in seq_len(p_d)) {
    df[[paste0("mD", j)]] <- mD[, j]
    df[[paste0("rD", j)]] <- D[, j] - mD[, j]
  }
  # Constant across observations: the errors are independent of (X, D).
  df$density_oracle <- rep(.plqr_ald_density0(tau, sigma), n)

  # Optional injection of non-finites for guard testing
  if (!is.null(inject)) {
    if (!inject %in% c("na", "inf")) {
      .abort_bad_arg("inject", "must be one of NULL, 'na', or 'inf'", quoted = FALSE)
    }
    idx_inject <- .validate_scalar_int(idx_inject, arg = "idx_inject", min = 1L)
    if (idx_inject > n) {
      .abort_bad_arg("idx_inject", sprintf("must be <= n (got idx_inject=%d, n=%d)", idx_inject, n))
    }
    allowed_inject_cols <- c(
      "Y", "g_oracle", "density_oracle",
      paste0("D", seq_len(p_d)[1]),
      paste0("mD", seq_len(p_d)[1]),
      paste0("rD", seq_len(p_d)[1])
    )
    if (!inject_col %in% allowed_inject_cols) {
      .dnn_abort("argument_invalid", sprintf(
        "Column '%s' (specified by `inject_col`) is not a valid injection target.\nValid options: %s",
        inject_col, paste(allowed_inject_cols, collapse = ", ")
      ), argument = "inject_col")
    }
    if (identical(inject, "na")) df[[inject_col]][idx_inject] <- NA_real_
    if (identical(inject, "inf")) df[[inject_col]][idx_inject] <- Inf
  }

  # Ground-truth closures
  mD_fn <- function(xx) {
    xx <- as.matrix(xx)
    xx %*% Gamma
  }

  list(
    data = df,
    x0 = as.numeric(x0),
    truth = list(
      theta_true = as.numeric(theta_true),
      g_fn = function(xx) g_fn(as.matrix(xx)),
      mD_fn = mD_fn,
      density_fn = .plqr_ald_density_fn(tau, sigma)
    ),
    g_x0 = as.numeric(g_fn(matrix(x0, nrow = 1)))
  )
}

#' Quick helper: sample a batch of subsets (for U-stat / localization experiments) (experimental)
#' @section Lifecycle:
#' **Experimental.** The PLQR API and inference interfaces may change.
#'
#' @inheritParams cate_demo_subsets
#' @inherit cate_demo_subsets return
#' @examples
#' plqr_demo_subsets(n = 10, K = 4, s = 3, seed = 1)
#' @export
plqr_demo_subsets <- function(n, K = 1000L, s = 8L, seed = NULL) {
  .rdistnn_demo_subsets(n = n, K = K, s = s, seed = seed)
}
