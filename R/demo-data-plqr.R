# R/demo-data-plqr.R

# --- Internal: Asymmetric Laplace draws (location 0, scale = sigma, asymmetry = tau) ----
# Ensures Q_tau(U) = 0, i.e., the tau-quantile of U is zero.
.r_ald <- function(n, tau = 0.5, sigma = 1) {
  stopifnot(0 < tau, tau < 1, sigma > 0)
  u <- runif(n)
  out <- ifelse(
    u < tau,
    sigma * ((1 - tau) * log(u / tau)),
    -sigma * (tau * log((1 - u) / (1 - tau)))
  )
  out
}

#' Minimal PLQR demo data (oracle nuisances; no internal fitting)
#'
#' Generates data for a partially linear quantile regression (PLQR) model
#' \deqn{
#'   Q_{\tau}(Y \mid X, D) \;=\; D^{\top}\,\beta_{0}(\tau) \;+\; q_{\tau}(X) \, ,
#' }
#' with a user-chosen design for \eqn{q_{\tau}}. Errors follow an Asymmetric
#' Laplace distribution (ALD) at level \eqn{\tau}, so the conditional
#' \eqn{\tau}-quantile equals the right-hand side.
#'
#' The function also returns the oracle
#' \eqn{m_{D}(X) \equiv \mathbb{E}[D \mid X]} to build the orthogonal PLQR score
#' \deqn{
#'   \psi(W;\beta)
#'   \;=\;
#'   (\tau - \mathbf{1}\{\,Y \le D^{\top}\beta + q_{\tau}(X)\,\})\,
#'   (D - m_{D}(X))\, .
#' }
#'
#' @param n Integer; sample size.
#' @param p Integer; number of features in \eqn{X} (recommend \eqn{p \ge 2}).
#' @param p_d Integer; dimension of \eqn{D} (number of partially linear regressors).
#' @param tau Quantile level in \eqn{(0,1)}.
#' @param sigma Positive ALD scale for the regression error (larger means noisier \eqn{Y}).
#' @param design Character; one of \code{"linear"}, \code{"sinusoidal"}, \code{"bump"} for \eqn{q_{\tau}(X)}.
#' @param seed Optional integer seed.
#' @param x0 Optional target point in \eqn{\mathbb{R}^{p}} for localization (defaults to the zero vector).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{data}: a \code{data.frame} with columns \code{X1..Xp}, \code{D1..Dp_d}, \code{Y};
#'         oracle nuisances \code{q_oracle} (the value of \eqn{q_{\tau}(X)}),
#'         \code{mD1..mDp_d} (components of \eqn{m_{D}(X)}),
#'         and residuals \code{rD1..rDp_d} (components of \eqn{D - m_{D}(X)}).
#'   \item \code{x0}: the target point used.
#'   \item \code{truth}: a list with \code{beta_true} (length \code{p_d}),
#'         \code{q_fn(x)} and \code{mD_fn(x)} closures to evaluate ground truth at new \eqn{X}.
#'   \item \code{q_x0}: the true value \eqn{q_{\tau}(x_{0})}.
#' }
#'
#' @examples
#' demo <- plqr_demo_data(
#'   n = 300, p = 2, p_d = 1, tau = 0.5,
#'   design = "sinusoidal", seed = 1
#' )
#' head(demo$data)
#' demo$truth$beta_true
#' demo$q_x0
#'
#' @export
plqr_demo_data <- function(n = 300,
                           p = 2,
                           p_d = 1,
                           tau = 0.5,
                           sigma = 0.5,
                           design = c("linear", "sinusoidal", "bump"),
                           seed = NULL,
                           x0 = NULL) {
  stopifnot(p >= 1L, p_d >= 1L, 0 < tau, tau < 1, sigma > 0)
  design <- match.arg(design)
  if (!is.null(seed)) set.seed(seed)
  if (is.null(x0)) x0 <- rep(0, p)

  # Covariates
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", seq_len(p))

  # True beta (partially linear effect at quantile tau)
  beta_true <- rep(0.5, p_d)

  # Nonparametric quantile component q_tau(X)
  q_fn <- switch(design,
    linear = function(xx) 0.5 * xx[, 1],
    sinusoidal = function(xx) 0.5 * sin(xx[, 1]) + 0.25 * (if (ncol(xx) >= 2) xx[, 2]^2 else 0),
    bump = {
      center <- c(0.5, -0.5, rep(0, max(0, p - 2L)))
      function(xx) {
        xx <- as.matrix(xx)
        d2 <- rowSums((xx - rep(center, each = nrow(xx)))^2)
        0.2 + 0.8 * exp(-0.5 * d2)
      }
    }
  )
  qX <- q_fn(X)

  # Partially linear regressors D with oracle E(D given X) = mD(X)
  # We'll use a simple linear mD(X) with different loadings per column of D.
  Gamma <- matrix(0, nrow = p, ncol = p_d)
  for (j in seq_len(p_d)) {
    g <- c(0.6, -0.4, rep(0, max(0, p - 2L))) + (j - 1) * 0.2
    Gamma[, j] <- g[seq_len(p)]
  }
  mD <- X %*% Gamma # n x p_d matrix, oracle E(D given X)
  V <- matrix(rnorm(n * p_d, sd = 0.6), n, p_d)
  D <- mD + V
  colnames(D) <- paste0("D", seq_len(p_d))

  # Error with tau-quantile zero
  U <- .r_ald(n, tau = tau, sigma = sigma)

  # Outcome satisfying PLQR structure:
  # Y = D^T beta_true + q_tau(X) + U   with Q_tau(U  given  X, D) = 0
  Y <- as.vector(D %*% beta_true) + qX + U

  # Tidy outputs
  df <- data.frame(X, D, Y = Y, q_oracle = qX)
  # add oracle mD and residual rD to df
  for (j in seq_len(p_d)) {
    df[[paste0("mD", j)]] <- mD[, j]
    df[[paste0("rD", j)]] <- D[, j] - mD[, j]
  }

  # Ground-truth closures
  mD_fn <- function(xx) {
    xx <- as.matrix(xx)
    xx %*% Gamma
  }

  list(
    data = df,
    x0 = as.numeric(x0),
    truth = list(beta_true = beta_true, q_fn = function(xx) q_fn(as.matrix(xx)), mD_fn = mD_fn),
    q_x0 = as.numeric(q_fn(matrix(x0, nrow = 1)))
  )
}

#' Orthogonal PLQR moment contributions for a given \eqn{\beta} (oracle nuisances)
#'
#' Computes, for each observation \eqn{i}, the orthogonal PLQR moment
#' \deqn{
#'   (\tau - \mathbf{1}\{\,Y_{i} \le D_{i}^{\top}\beta + q_{\tau}(X_{i})\,\})\,
#'   (D_{i} - m_{D}(X_{i})) \, ,
#' }
#' where \eqn{m_{D}(X) = \mathbb{E}[D \mid X]} is provided as an oracle nuisance.
#' This is the vector moment entering the Z-estimation.
#'
#' @param df A \code{data.frame} as returned by \code{plqr_demo_data()$data},
#'   containing columns \code{Y}, \code{D1..Dp_d}, oracle nuisances
#'   \code{mD1..mDp_d} (components of \eqn{m_{D}(X)}), and \code{q_oracle}
#'   (the value of \eqn{q_{\tau}(X)}).
#' @param beta Numeric length-\code{p_d} candidate coefficient vector \eqn{\beta}.
#' @param tau Quantile level \eqn{\tau} (should match the data-generating level).
#'
#' @return A numeric matrix of size \eqn{n \times p_d} with the moment
#'   contributions for each observation (rows) and each partially linear
#'   component (columns).
#'
#' @examples
#' demo <- plqr_demo_data(n = 200, p = 2, p_d = 2, tau = 0.5, seed = 1)
#' g <- plqr_moment(df = demo$data, beta = demo$truth$beta_true, tau = 0.5)
#' colMeans(g) # should be near zero at truth
#'
#' @export
plqr_moment <- function(df, beta, tau) {
  p_d <- sum(grepl("^D[0-9]+$", names(df)))
  stopifnot(length(beta) == p_d)
  Dmat <- as.matrix(df[paste0("D", seq_len(p_d))])
  mD <- as.matrix(df[paste0("mD", seq_len(p_d))])
  rD <- Dmat - mD
  thresh <- as.vector(Dmat %*% beta) + df$q_oracle
  ind <- as.numeric(df$Y <= thresh)
  resid <- (tau - ind)
  rD * resid
}

#' Quick helper: sample a batch of subsets (for U-stat / localization experiments)
#' @inheritParams cate_demo_subsets
#' @export
plqr_demo_subsets <- function(n, K = 1000L, s = 8L, seed = NULL) {
  stopifnot(s >= 1L, s <= n, K >= 1L)
  if (!is.null(seed)) set.seed(seed)
  t(replicate(K, sort(sample.int(n, s))))
}
