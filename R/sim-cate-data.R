# sim-cate-data.R
# Data-generating processes for the CATE simulation study.
#
# Three DGPs from the paper "Robust Inference for Conditional Z-Estimation via
# Distributional Nearest Neighbors":
#
#   dgp1_linear          -- linear tau, logistic propensity; correctly specified
#                           for OLS-interaction benchmark
#   dgp2_nonlinear       -- sinusoidal tau; same propensity; misspecified for OLS
#   dgp3_sparse_overlap  -- sinusoidal tau; propensity near 0/1 everywhere
#                           except a Gaussian neighbourhood of x0
#
# All DGPs supply oracle nuisances (e, m0, m1, tau, rho) so that
# DGP isolation is possible; first-stage estimation error is deferred.
#
# Public API: sim_cate_data()   (exported)
# Internals:  .sim_dgp1(), .sim_dgp2(), .sim_dgp3()   (not exported)

# ── internal helpers ─────────────────────────────────────────────────────────

.sim_cate_m0 <- function(X) 1 + X[, 1] - X[, 2]

.sim_cate_pseudo <- function(D, Y, e, m) (D - e) * (Y - m) / (e * (1 - e))

.sim_cate_assemble <- function(X, D, Y, e_true, m0_true, tau_true,
                                tau_fn, e_fn, x0) {
  m1_true <- m0_true + tau_true
  m_true  <- m0_true + e_true * tau_true
  rho     <- .sim_cate_pseudo(D, Y, e_true, m_true)
  p       <- ncol(X)

  df <- data.frame(
    X,
    D          = D,
    Y          = as.numeric(Y),
    e_oracle   = as.numeric(e_true),
    m_oracle   = as.numeric(m_true),
    m0_oracle  = as.numeric(m0_true),
    m1_oracle  = as.numeric(m1_true),
    tau_oracle = as.numeric(tau_true),
    rho_oracle = as.numeric(rho)
  )
  colnames(df)[seq_len(p)] <- paste0("X", seq_len(p))

  list(
    data   = df,
    x0     = as.numeric(x0),
    tau_x0 = as.numeric(tau_fn(matrix(x0, nrow = 1L))),
    truth  = list(
      tau_fn = tau_fn,
      e_fn   = e_fn,
      m0_fn  = .sim_cate_m0,
      m_fn   = function(xx) {
        xx <- as.matrix(xx)
        .sim_cate_m0(xx) + e_fn(xx) * tau_fn(xx)
      }
    )
  )
}

# ── DGP 1: linear CATE ───────────────────────────────────────────────────────

.sim_dgp1 <- function(n, kx = 2L, sigma = 0.5, clip = 1e-3, seed = NULL,
                       x0 = rep(0, kx)) {
  if (!is.null(seed)) set.seed(seed)
  X <- matrix(stats::rnorm(n * kx), n, kx)

  beta_tau <- c(0.5, -0.3, rep(0, max(0L, kx - 2L)))
  tau_fn   <- function(xx) {
    xx <- as.matrix(xx)
    as.numeric(xx %*% beta_tau[seq_len(ncol(xx))])
  }

  gamma <- c(0.5, -0.5, rep(0, max(0L, kx - 2L)))
  e_fn  <- function(xx) {
    xx <- as.matrix(xx)
    pmin(pmax(stats::plogis(as.numeric(xx %*% gamma[seq_len(ncol(xx))])),
              clip), 1 - clip)
  }

  tau_true <- tau_fn(X)
  e_true   <- e_fn(X)
  m0_true  <- .sim_cate_m0(X)
  D        <- stats::rbinom(n, 1L, e_true)
  Y        <- m0_true + D * tau_true + stats::rnorm(n, sd = sigma)

  .sim_cate_assemble(X, D, Y, e_true, m0_true, tau_true, tau_fn, e_fn, x0)
}

# ── DGP 2: nonlinear (sinusoidal) CATE ───────────────────────────────────────

.sim_dgp2 <- function(n, kx = 2L, sigma = 0.5, clip = 1e-3, seed = NULL,
                       x0 = rep(0, kx)) {
  if (!is.null(seed)) set.seed(seed)
  X <- matrix(stats::rnorm(n * kx), n, kx)

  tau_fn <- function(xx) {
    xx <- as.matrix(xx)
    0.5 + 0.5 * sin(xx[, 1])
  }

  gamma <- c(0.5, -0.5, rep(0, max(0L, kx - 2L)))
  e_fn  <- function(xx) {
    xx <- as.matrix(xx)
    pmin(pmax(stats::plogis(as.numeric(xx %*% gamma[seq_len(ncol(xx))])),
              clip), 1 - clip)
  }

  tau_true <- tau_fn(X)
  e_true   <- e_fn(X)
  m0_true  <- .sim_cate_m0(X)
  D        <- stats::rbinom(n, 1L, e_true)
  Y        <- m0_true + D * tau_true + stats::rnorm(n, sd = sigma)

  .sim_cate_assemble(X, D, Y, e_true, m0_true, tau_true, tau_fn, e_fn, x0)
}

# ── DGP 3: sparse overlap ────────────────────────────────────────────────────

# e(x) interpolates between 0.5 (at x0) and logistic(c * gamma'x) (far from
# x0) via a Gaussian weight w(x) = exp(-||x - x0||^2 / (2 r^2)).  Tests
# methods that do not adapt to the local effective sample size.

.sim_dgp3 <- function(n, kx = 2L, sigma = 0.5, clip = 1e-3, seed = NULL,
                       x0 = rep(0, kx), r = 1.0, c_sharp = 3.0) {
  if (!is.null(seed)) set.seed(seed)
  X <- matrix(stats::rnorm(n * kx), n, kx)

  tau_fn <- function(xx) {
    xx <- as.matrix(xx)
    0.5 + 0.5 * sin(xx[, 1])
  }

  gamma <- c(0.5, -0.5, rep(0, max(0L, kx - 2L)))
  e_fn  <- function(xx) {
    xx  <- as.matrix(xx)
    d2  <- rowSums((xx - rep(x0, each = nrow(xx)))^2)
    w   <- exp(-d2 / (2 * r^2))
    lin <- as.numeric(xx %*% gamma[seq_len(ncol(xx))])
    e   <- w * 0.5 + (1 - w) * stats::plogis(c_sharp * lin)
    pmin(pmax(e, clip), 1 - clip)
  }

  tau_true <- tau_fn(X)
  e_true   <- e_fn(X)
  m0_true  <- .sim_cate_m0(X)
  D        <- stats::rbinom(n, 1L, e_true)
  Y        <- m0_true + D * tau_true + stats::rnorm(n, sd = sigma)

  .sim_cate_assemble(X, D, Y, e_true, m0_true, tau_true, tau_fn, e_fn, x0)
}

# ── Public dispatcher ─────────────────────────────────────────────────────────

#' Draw one sample from a named CATE simulation DGP
#'
#' Generates a synthetic dataset for the CATE simulation study, including
#' oracle nuisance quantities \eqn{e}, \eqn{m_0}, \eqn{m_1}, and the
#' Robinson pseudo-outcome \eqn{\rho}, so that first-stage estimation error
#' can be isolated from the target inference step.
#'
#' Three designs are available:
#' \describe{
#'   \item{`dgp1_linear`}{Linear CATE \eqn{\tau(x) = \beta'x}; logistic
#'     propensity. Correctly specified for the OLS-interaction benchmark.}
#'   \item{`dgp2_nonlinear`}{Sinusoidal CATE \eqn{\tau(x) = 0.5 + 0.5\sin(x_1)};
#'     same propensity. Misspecified for OLS; primary horse-race DGP.}
#'   \item{`dgp3_sparse_overlap`}{Same \eqn{\tau} as DGP 2; propensity
#'     interpolates between 0.5 at `x0` and a near-boundary logistic far from
#'     `x0`. Tests estimators that do not adapt to local effective sample size.}
#' }
#'
#' @param design Character; one of `"dgp1_linear"`, `"dgp2_nonlinear"`,
#'   `"dgp3_sparse_overlap"`.
#' @param n Integer sample size.
#' @param kx Integer covariate dimension (must be \eqn{\ge 2}).
#' @param sigma Noise standard deviation in \eqn{Y}.
#' @param clip Propensity clipping bound in \eqn{(0, 0.5)}.
#' @param seed Optional integer seed for reproducibility.
#' @param x0 Optional numeric length-`kx` evaluation point (defaults to the
#'   origin).
#' @param ... Extra arguments forwarded to the underlying DGP (e.g. `r` and
#'   `c_sharp` for `dgp3_sparse_overlap`).
#'
#' @return A list with:
#' \itemize{
#'   \item `data`: a `data.frame` with columns `X1..Xkx`, `D`, `Y`, and oracle
#'     nuisances `e_oracle`, `m_oracle`, `m0_oracle`, `m1_oracle`,
#'     `tau_oracle`, `rho_oracle`.
#'   \item `x0`: numeric length-`kx` evaluation point.
#'   \item `tau_x0`: true treatment effect at `x0`.
#'   \item `truth`: named list of closures `tau_fn`, `e_fn`, `m0_fn`, `m_fn`.
#' }
#'
#' @examples
#' dgp <- sim_cate_data("dgp2_nonlinear", n = 200, kx = 2, seed = 1)
#' head(dgp$data)
#' dgp$tau_x0
#'
#' @export
sim_cate_data <- function(design, n, kx = 2L, sigma = 0.5, clip = 1e-3,
                           seed = NULL, x0 = NULL, ...) {
  stopifnot(kx >= 2L, clip > 0, clip < 0.5)
  if (is.null(x0)) x0 <- rep(0, kx)
  switch(design,
    dgp1_linear         = .sim_dgp1(n, kx, sigma, clip, seed, x0, ...),
    dgp2_nonlinear      = .sim_dgp2(n, kx, sigma, clip, seed, x0, ...),
    dgp3_sparse_overlap = .sim_dgp3(n, kx, sigma, clip, seed, x0, ...),
    stop("Unknown design: '", design,
         "'. Choose from: dgp1_linear, dgp2_nonlinear, dgp3_sparse_overlap.")
  )
}
