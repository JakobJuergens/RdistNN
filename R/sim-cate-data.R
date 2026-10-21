# sim-cate-data.R
# Data-generating processes for the CATE simulation study.
#
# Six DGPs from the companion conditional Z-estimation paper:
#
#   dgp1_linear          -- linear tau, logistic propensity; correctly specified
#                           for OLS-interaction benchmark
#   dgp2_nonlinear       -- sharper sinusoidal tau; same propensity;
#                           misspecified for OLS
#   dgp3_sparse_overlap  -- same sharper tau; overlap concentrated in a tight
#                           Gaussian neighbourhood of x0
#   dgp4_clustered       -- same sharper tau; three-component mixture covariates
#                           with a dense inner spike and two remote clusters;
#                           tests density adaptivity
#   dgp5_correlated      -- curved tau; strongly correlated covariates (rho=0.85);
#                           x0 = (0.75, 0); tests shape adaptivity
#   dgp6_nuisance_efficiency
#                         -- correctly specified logistic/linear nuisances with
#                            nonlinear tau; isolates full-sample vs cross-fit
#                            nuisance-efficiency tradeoff
#
# All DGPs supply oracle nuisances (e, m0, m1, tau) so that
# DGP isolation is possible; first-stage estimation error is deferred.
#
# Public API: sim_cate_data()   (exported)
# Internals:  .sim_dgp1(), ..., .sim_dgp6()   (not exported)

# ── internal helpers ─────────────────────────────────────────────────────────

.sim_cate_m0 <- function(X) 1 + X[, 1] - X[, 2]

.sim_cate_aipw <- function(D, Y, e, m0, m1) {
  m1 - m0 + D * (Y - m1) / e - (1 - D) * (Y - m0) / (1 - e)
}

.sim_cate_tau_shared <- function(xx) {
  xx <- as.matrix(xx)
  0.5 + 0.75 * sin(1.5 * xx[, 1L])
}

.sim_cate_assemble <- function(X, D, Y, e_true, m0_true, tau_true,
                                tau_fn, e_fn, x0, sigma,
                                m0_fn = .sim_cate_m0) {
  m1_true <- m0_true + tau_true
  m_true  <- m0_true + e_true * tau_true
  aipw    <- .sim_cate_aipw(D, Y, e_true, m0_true, m1_true)
  p       <- ncol(X)
  e_x0    <- as.numeric(e_fn(matrix(x0, nrow = 1L)))
  score_sd_x0 <- as.numeric(sigma * sqrt(1 / e_x0 + 1 / (1 - e_x0)))

  df <- data.frame(
    X,
    D          = D,
    Y          = as.numeric(Y),
    e_oracle   = as.numeric(e_true),
    m_oracle   = as.numeric(m_true),
    m0_oracle  = as.numeric(m0_true),
    m1_oracle  = as.numeric(m1_true),
    tau_oracle = as.numeric(tau_true),
    aipw_oracle = as.numeric(aipw)
  )
  colnames(df)[seq_len(p)] <- paste0("X", seq_len(p))

  list(
    data   = df,
    x0     = as.numeric(x0),
    tau_x0 = as.numeric(tau_fn(matrix(x0, nrow = 1L))),
    truth  = list(
      tau_fn = tau_fn,
      e_fn   = e_fn,
      m0_fn  = m0_fn,
      sigma = as.numeric(sigma),
      score_sd_x0 = score_sd_x0,
      m_fn   = function(xx) {
        xx <- as.matrix(xx)
        m0_fn(xx) + e_fn(xx) * tau_fn(xx)
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

  .sim_cate_assemble(
    X, D, Y, e_true, m0_true, tau_true, tau_fn, e_fn, x0,
    sigma = sigma
  )
}

# ── DGP 2: nonlinear (sinusoidal) CATE ───────────────────────────────────────

.sim_dgp2 <- function(n, kx = 2L, sigma = 0.5, clip = 1e-3, seed = NULL,
                       x0 = c(1, rep(0, max(0L, kx - 1L)))) {
  if (!is.null(seed)) set.seed(seed)
  X <- matrix(stats::rnorm(n * kx), n, kx)

  tau_fn <- .sim_cate_tau_shared

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

  .sim_cate_assemble(
    X, D, Y, e_true, m0_true, tau_true, tau_fn, e_fn, x0,
    sigma = sigma
  )
}

# ── DGP 3: sparse overlap ────────────────────────────────────────────────────

# e(x) interpolates between 0.5 (at x0) and logistic(c * gamma'x) (far from
# x0) via a Gaussian weight w(x) = exp(-||x - x0||^2 / (2 r^2)).  Tests
# methods that do not adapt to the local effective sample size.

.sim_dgp3 <- function(n, kx = 2L, sigma = 0.5, clip = 1e-3, seed = NULL,
                       x0 = c(1, rep(0, max(0L, kx - 1L))), r = 0.55,
                       c_sharp = 5.5) {
  if (!is.null(seed)) set.seed(seed)
  X <- matrix(stats::rnorm(n * kx), n, kx)

  tau_fn <- .sim_cate_tau_shared

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

  .sim_cate_assemble(
    X, D, Y, e_true, m0_true, tau_true, tau_fn, e_fn, x0,
    sigma = sigma
  )
}

# ── DGP 4: clustered (mixture) covariate design ──────────────────────────────

# X is a three-component mixture: a dense inner spike N(x0, inner_sd^2 * I)
# centred at x0 (weight 0.8) plus two diffuse remote clusters away from x0
# (weights 0.1 each).  Tests density adaptivity: the DNN s-NN effective
# bandwidth automatically shrinks at x0 because the local density is high,
# while a fixed-bandwidth (Silverman) kernel uses a bandwidth driven by the
# global spread and therefore over-smooths at x0.  tau and e are the same
# functional forms as dgp2_nonlinear to allow direct comparison; only the
# covariate distribution changes.

.sim_dgp4 <- function(n, kx = 2L, sigma = 0.5, clip = 1e-3, seed = NULL,
                       x0 = c(1, rep(0, max(0L, kx - 1L))),
                       mix_probs = c(0.8, 0.1, 0.1), inner_sd = 0.20,
                       outer_sd = 1.25,
                       outer_mu_pos = c(4.5, 2.0, rep(0, max(0L, kx - 2L))),
                       outer_mu_neg = c(-2.5, -2.0, rep(0, max(0L, kx - 2L)))) {
  if (!is.null(seed)) set.seed(seed)

  if (!is.numeric(mix_probs) || length(mix_probs) != 3L ||
      any(!is.finite(mix_probs)) || any(mix_probs <= 0)) {
    .abort_bad_arg("mix_probs", "must be a positive finite numeric vector of length 3")
  }
  if (!is.numeric(inner_sd) || length(inner_sd) != 1L || !is.finite(inner_sd) || inner_sd <= 0) {
    .abort_bad_arg("inner_sd", "must be a positive finite numeric scalar")
  }
  if (!is.numeric(outer_sd) || length(outer_sd) != 1L || !is.finite(outer_sd) || outer_sd <= 0) {
    .abort_bad_arg("outer_sd", "must be a positive finite numeric scalar")
  }
  outer_mu_pos <- .validate_target_vector(outer_mu_pos, p = kx,
                                          x0_arg = "outer_mu_pos", x_arg = "X", finite = TRUE)
  outer_mu_neg <- .validate_target_vector(outer_mu_neg, p = kx,
                                          x0_arg = "outer_mu_neg", x_arg = "X", finite = TRUE)

  mix_probs <- as.numeric(mix_probs / sum(mix_probs))
  cluster <- sample.int(3L, size = n, replace = TRUE, prob = mix_probs)
  n_inner <- sum(cluster == 1L)
  n_pos   <- sum(cluster == 2L)
  n_neg   <- sum(cluster == 3L)

  X <- matrix(0, n, kx)
  if (n_inner > 0L)
    X[cluster == 1L, ] <- matrix(stats::rnorm(n_inner * kx, sd = inner_sd),
                                 n_inner, kx) +
                          matrix(x0, n_inner, kx, byrow = TRUE)
  if (n_pos > 0L)
    X[cluster == 2L, ] <- matrix(stats::rnorm(n_pos * kx, sd = outer_sd),
                                 n_pos, kx) +
                          matrix(outer_mu_pos, n_pos, kx, byrow = TRUE)
  if (n_neg > 0L)
    X[cluster == 3L, ] <- matrix(stats::rnorm(n_neg * kx, sd = outer_sd),
                                 n_neg, kx) +
                          matrix(outer_mu_neg, n_neg, kx, byrow = TRUE)

  tau_fn <- .sim_cate_tau_shared

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

  .sim_cate_assemble(
    X, D, Y, e_true, m0_true, tau_true, tau_fn, e_fn, x0,
    sigma = sigma
  )
}

# ── DGP 5: correlated covariates ─────────────────────────────────────────────

# X ~ N(0, Sigma) where Sigma = [[1, rho], [rho, 1]] (high positive correlation
# between X1 and X2; higher dimensions are iid N(0,1)).  Tests shape adaptivity:
# the L2-ball neighbourhood around x0 is geometrically misaligned with the
# elliptical isocontours of the covariate distribution, causing the Nadaraya-
# Watson (Silverman diagonal) kernel and DNN-L2 to include a biased mix of
# observations.  DNN with Mahalanobis distance (dnn_oracle_mahal) forms the
# correct ellipsoidal neighbourhood and achieves lower bias.
# x0 is set off-centre (x1 = 0.75) so neighbourhood asymmetry is non-trivial:
# at x0 = (0,0) the symmetry of the design would cancel the bias.

.sim_dgp5 <- function(n, kx = 2L, sigma = 0.5, clip = 1e-3, seed = NULL,
                       x0 = c(0.75, rep(0, max(0L, kx - 1L))), rho = 0.85) {
  if (!is.null(seed)) set.seed(seed)

  # Correlated draw: X1 = Z1, X2 = rho*Z1 + sqrt(1-rho^2)*Z2; Xj = Zj for j>2.
  Z <- matrix(stats::rnorm(n * kx), n, kx)
  X <- Z
  if (kx >= 2L)
    X[, 2L] <- rho * Z[, 1L] + sqrt(1 - rho^2) * Z[, 2L]

  tau_fn <- function(xx) {
    xx <- as.matrix(xx)
    sin(2 * xx[, 1L]) + 0.5 * xx[, 1L]^2
  }

  e_fn <- function(xx) {
    xx  <- as.matrix(xx)
    lin <- xx[, 1L]
    if (ncol(xx) >= 2L) lin <- lin - xx[, 2L]
    pmin(pmax(stats::plogis(lin), clip), 1 - clip)
  }

  tau_true <- tau_fn(X)
  e_true   <- e_fn(X)
  m0_true  <- .sim_cate_m0(X)
  D        <- stats::rbinom(n, 1L, e_true)
  Y        <- m0_true + D * tau_true + stats::rnorm(n, sd = sigma)

  .sim_cate_assemble(
    X, D, Y, e_true, m0_true, tau_true, tau_fn, e_fn, x0,
    sigma = sigma
  )
}

# ── DGP 6: nuisance-efficiency design ────────────────────────────────────────

# Propensity and the marginal outcome nuisance are aligned with
# simple benchmark models.  The DNN simulation methods now use AIPW potential
# outcome nuisances, so this design mainly remains a full-sample versus
# cross-fit precision stress test rather than a fully correctly specified
# AIPW-nuisance design.  tau remains nonlinear so the target inference problem
# is still genuinely local.

.sim_dgp6 <- function(n, kx = 2L, sigma = 0.5, clip = 1e-3, seed = NULL,
                       x0 = c(1, rep(0, max(0L, kx - 1L)))) {
  if (!is.null(seed)) set.seed(seed)
  X <- matrix(stats::rnorm(n * kx), n, kx)

  tau_fn <- function(xx) {
    xx <- as.matrix(xx)
    0.5 + 0.60 * sin(1.5 * xx[, 1L])
  }

  e_fn <- function(xx) {
    xx <- as.matrix(xx)
    lin <- 1.25 * (xx[, 1L] - 1) - 0.75 * xx[, 2L]
    pmin(pmax(stats::plogis(lin), clip), 1 - clip)
  }

  m_fn <- function(xx) {
    xx <- as.matrix(xx)
    1 + 1.5 * (xx[, 1L] - 1) - xx[, 2L]
  }

  m0_fn <- function(xx) {
    xx <- as.matrix(xx)
    m_fn(xx) - e_fn(xx) * tau_fn(xx)
  }

  tau_true <- tau_fn(X)
  e_true   <- e_fn(X)
  m0_true  <- m0_fn(X)
  D        <- stats::rbinom(n, 1L, e_true)
  Y        <- m0_true + D * tau_true + stats::rnorm(n, sd = sigma)

  .sim_cate_assemble(X, D, Y, e_true, m0_true, tau_true,
                     tau_fn, e_fn, x0, m0_fn = m0_fn, sigma = sigma)
}

# ── Public dispatcher ─────────────────────────────────────────────────────────

#' Draw one sample from a named CATE simulation DGP
#'
#' Generates a synthetic dataset for the CATE simulation study, including
#' oracle nuisance quantities \eqn{e}, \eqn{m_0}, \eqn{m_1}, the AIPW score used
#' by the general DNN CATE benchmarks.
#'
#' Six designs are available:
#' \describe{
#'   \item{`dgp1_linear`}{Linear CATE \eqn{\tau(x) = \beta'x}; logistic
#'     propensity. Correctly specified for the OLS-interaction benchmark.}
#'   \item{`dgp2_nonlinear`}{Sharpened nonlinear CATE
#'     \eqn{\tau(x) = 0.5 + 0.75\sin(1.5 x_1)};
#'     same propensity. Default
#'     target is the interior point \eqn{x_0 = (1,0,\ldots)} so the local CATE
#'     is visibly nonlinear around the evaluation point.}
#'   \item{`dgp3_sparse_overlap`}{Same \eqn{\tau} as DGP 2; propensity
#'     interpolates between 0.5 at `x0` and a much sharper near-boundary
#'     logistic far from `x0`, using a tighter Gaussian overlap pocket. Default
#'     target is the same interior point \eqn{x_0 = (1,0,\ldots)}.}
#'   \item{`dgp4_clustered`}{Same \eqn{\tau} and propensity as DGP 2; covariates
#'     from a three-component mixture: 80\% \eqn{N(x_0, 0.2^2 I)} (dense spike at
#'     `x0`) plus two diffuse 10\% remote clusters centered at
#'     \eqn{(4.5,2,0,\ldots)} and \eqn{(-2.5,-2,0,\ldots)} with covariance
#'     \eqn{1.25^2 I}. Observations are independent draws from this mixture;
#'     the name does not indicate clustered sampling. Default target is
#'     \eqn{x_0 = (1,0,\ldots)} and the inner cluster is recentered there.}
#'   \item{`dgp5_correlated`}{Covariates \eqn{X \sim N(0,\Sigma)} with
#'     \eqn{\Sigma_{12} = \rho = 0.85}; \eqn{\tau(x) = \sin(2x_1) + 0.5x_1^2};
#'     evaluation point at \eqn{x_0 = (0.75, 0)}. The `dnn_oracle_mahal`
#'     benchmark uses Mahalanobis distance to account for off-diagonal covariance.}
#'   \item{`dgp6_nuisance_efficiency`}{Propensity
#'     \eqn{e(x) = \operatorname{logit}^{-1}(1.25(x_1 - 1) - 0.75x_2)} and
#'     marginal regression \eqn{m(x) = 1 + 1.5(x_1 - 1) - x_2} are aligned with
#'     simple benchmark nuisance models, while
#'     \eqn{\tau(x) = 0.5 + 0.60\sin(1.5 x_1)} remains nonlinear. Tests the
#'     full-sample and cross-fitted nuisance workflows. Default target
#'     is \eqn{x_0 = (1,0,\ldots)}.}
#' }
#'
#' @param design Character; one of `"dgp1_linear"`, `"dgp2_nonlinear"`,
#'   `"dgp3_sparse_overlap"`, `"dgp4_clustered"`, `"dgp5_correlated"`, or
#'   `"dgp6_nuisance_efficiency"`.
#' @param n Integer sample size.
#' @param kx Integer covariate dimension (must be \eqn{\ge 2}).
#' @param sigma Noise standard deviation in \eqn{Y}.
#' @param clip Propensity clipping bound in \eqn{(0, 0.5)}.
#' @param seed Optional integer seed, used locally without changing the caller's
#'   RNG state. `NULL` uses the current stream.
#' @param x0 Optional numeric length-`kx` evaluation point. Defaults are
#'   DGP-specific interior points: \eqn{(0,\ldots,0)} for `dgp1_linear`,
#'   \eqn{(1,0,\ldots)} for `dgp2_nonlinear`, `dgp3_sparse_overlap`, and
#'   `dgp4_clustered`, and `dgp6_nuisance_efficiency`, and
#'   \eqn{(0.75,0,\ldots)} for `dgp5_correlated`.
#' @param ... Extra arguments forwarded to the underlying DGP (e.g. `r` and
#'   `c_sharp` for `dgp3_sparse_overlap`).
#'
#' @return A list with:
#' \itemize{
#'   \item `data`: a `data.frame` with columns `X1..Xkx`, `D`, `Y`, and oracle
#'     nuisances `e_oracle`, `m_oracle`, `m0_oracle`, `m1_oracle`,
#'     `tau_oracle`, and `aipw_oracle`.
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
  n <- .validate_scalar_int(n, arg = "n", min = 1L)
  kx <- .validate_scalar_int(kx, arg = "kx", min = 2L)
  if (!is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma) || sigma < 0) {
    .abort_bad_arg("sigma", "must be a non-negative finite numeric scalar")
  }
  sigma <- as.numeric(sigma)
  if (!is.numeric(clip) || length(clip) != 1L || !is.finite(clip) || clip <= 0 || clip >= 0.5) {
    .abort_bad_arg("clip", "must be a finite numeric scalar in (0, 0.5)")
  }
  clip <- as.numeric(clip)
  if (!is.null(x0)) {
    x0 <- .validate_target_vector(x0, p = kx, x0_arg = "x0", x_arg = "X", finite = TRUE)
  }

  if (!is.null(seed)) {
    restore_rng <- .rdistnn_rng_restore()
    on.exit(restore_rng(), add = TRUE)
  }

  # Each DGP has its own natural default for x0; resolve here so callers that
  # pass x0 = NULL still get the DGP-appropriate evaluation point rather than a
  # generic rep(0, kx) that would be wrong for, e.g., dgp5_correlated (x0=(0.75,0)).
  switch(design,
    dgp1_linear         = .sim_dgp1(n, kx, sigma, clip, seed,
                            x0 = if (is.null(x0)) rep(0, kx) else x0, ...),
    dgp2_nonlinear      = .sim_dgp2(n, kx, sigma, clip, seed,
                            x0 = if (is.null(x0)) c(1, rep(0, max(0L, kx - 1L))) else x0, ...),
    dgp3_sparse_overlap = .sim_dgp3(n, kx, sigma, clip, seed,
                            x0 = if (is.null(x0)) c(1, rep(0, max(0L, kx - 1L))) else x0, ...),
    dgp4_clustered      = .sim_dgp4(n, kx, sigma, clip, seed,
                            x0 = if (is.null(x0)) c(1, rep(0, max(0L, kx - 1L))) else x0, ...),
    dgp5_correlated     = .sim_dgp5(n, kx, sigma, clip, seed,
                            x0 = if (is.null(x0)) c(0.75, rep(0, max(0L, kx - 1L))) else x0,
                            ...),
    dgp6_nuisance_efficiency = .sim_dgp6(n, kx, sigma, clip, seed,
                            x0 = if (is.null(x0)) c(1, rep(0, max(0L, kx - 1L))) else x0,
                            ...),
    stop("Unknown design: '", design,
         "'. Choose from: dgp1_linear, dgp2_nonlinear, dgp3_sparse_overlap, ",
         "dgp4_clustered, dgp5_correlated, dgp6_nuisance_efficiency.", call. = FALSE)
  )
}
