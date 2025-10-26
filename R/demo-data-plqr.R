# R/demo-data-plqr.R

# --- Internal: Asymmetric Laplace draws (location 0, scale = sigma, asymmetry = tau) ----
# Ensures Q_tau(U) = 0, i.e., the tau-quantile of U is zero.
.r_ald <- function(n, tau = 0.5, sigma = 1) {
  stopifnot(0 < tau, tau < 1, sigma > 0)
  u <- runif(n)
  out <- ifelse(
    u < tau,
    sigma * ( (1 - tau) * log(u / tau) ),
    -sigma * ( tau * log((1 - u) / (1 - tau)) )
  )
  out
}

#' Minimal PLQR demo data (oracle nuisances, no internal fitting)
#'
#' Generates data for a partially linear quantile regression model
#' \deqn{ Q_\tau(Y \mid X, D) = D^\top \beta_0(\tau) + q_\tau(X), }
#' with user-chosen design for \eqn{q_\tau}. Errors are Asymmetric Laplace
#' so the conditional \eqn{\tau}-quantile is exactly the right-hand side.
#'
#' Also returns oracle \eqn{m_D(X)=E[D\mid X]} to build the *orthogonal* PLQR score
#' \deqn{ \psi(W;\beta) = (\,\tau - 1\{Y \le D^\top\beta + q_\tau(X)\}\,)\,(D - m_D(X)). }
#'
#' @param n Sample size.
#' @param p Number of features in X (>= 2 recommended).
#' @param p_d Dimension of D (number of partially linear regressors).
#' @param tau Quantile level in (0,1).
#' @param sigma ALD scale for the regression error (larger means noisier Y).
#' @param design One of "linear", "sinusoidal", "bump" for \eqn{q_\tau(X)}.
#' @param seed Optional integer seed.
#' @param x0 Optional target point in R^p for localization (defaults to 0 vector).
#' @return A list with:
#' \itemize{
#'   \item \code{data}: data.frame with columns \code{X1..Xp}, \code{D1..Dp_d}, \code{Y},
#'         oracle nuisances \code{q_oracle} (== q_tau(X)) and \code{mD1..mDp_d} (== E[D|X]),
#'         and residuals \code{rD1..rDp_d} (== D - mD).
#'   \item \code{x0}: the target point used.
#'   \item \code{truth}: list with \code{beta_true} (length p_d),
#'         \code{q_fn(x)}, and \code{mD_fn(x)} closures to evaluate ground truth at new X.
#'   \item \code{q_x0}: the true \eqn{q_\tau(x_0)}.
#' }
#' @examples
#' demo <- plqr_demo_data(n = 300, p = 2, p_d = 1, tau = 0.5, design = "sinusoidal", seed = 1)
#' head(demo$data)
#' demo$truth$beta_true
#' demo$q_x0
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
  q_fn <- switch(
    design,
    linear     = function(xx) 0.5 * xx[, 1],
    sinusoidal = function(xx) 0.5 * sin(xx[, 1]) + 0.25 * (if (ncol(xx) >= 2) xx[, 2]^2 else 0),
    bump       = {
      center <- c(0.5, -0.5, rep(0, max(0, p - 2L)))
      function(xx) {
        xx <- as.matrix(xx)
        d2 <- rowSums( (xx - rep(center, each = nrow(xx)))^2 )
        0.2 + 0.8 * exp(-0.5 * d2)
      }
    }
  )
  qX <- q_fn(X)

  # Partially linear regressors D with oracle E[D|X] = mD(X)
  # We'll use a simple linear mD(X) with different loadings per column of D.
  Gamma <- matrix(0, nrow = p, ncol = p_d)
  for (j in seq_len(p_d)) {
    g <- c(0.6, -0.4, rep(0, max(0, p - 2L))) + (j - 1) * 0.2
    Gamma[, j] <- g[seq_len(p)]
  }
  mD <- X %*% Gamma                 # n x p_d matrix, oracle E[D|X]
  V  <- matrix(rnorm(n * p_d, sd = 0.6), n, p_d)
  D  <- mD + V
  colnames(D) <- paste0("D", seq_len(p_d))

  # Error with tau-quantile zero
  U <- .r_ald(n, tau = tau, sigma = sigma)

  # Outcome satisfying PLQR structure:
  # Y = D^T beta_true + q_tau(X) + U   with Q_tau(U | X, D) = 0
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
    xx <- as.matrix(xx); xx %*% Gamma
  }

  list(
    data = df,
    x0 = as.numeric(x0),
    truth = list(beta_true = beta_true, q_fn = function(xx) q_fn(as.matrix(xx)), mD_fn = mD_fn),
    q_x0 = as.numeric(q_fn(matrix(x0, nrow = 1)))
  )
}

#' Orthogonal PLQR moment contributions for a given beta (uses oracle nuisances)
#'
#' For each i, returns \eqn{ ( \tau - 1\{Y_i \le D_i^\top \beta + q_\tau(X_i)\} ) (D_i - m_D(X_i)) }.
#' This is the vector moment entering the Z-estimation.
#'
#' @param df A data.frame returned by \code{plqr_demo_data()$data}.
#' @param beta Numeric length-\code{p_d} candidate coefficient vector.
#' @param tau Quantile level (must match the data-generating tau).
#' @return An n x p_d matrix of moment contributions.
#' @examples
#' demo <- plqr_demo_data(n=200, p=2, p_d=2, tau=0.5, seed=1)
#' g <- plqr_moment(df = demo$data, beta = demo$truth$beta_true, tau = 0.5)
#' colMeans(g)  # should be near zero at truth
#' @export
plqr_moment <- function(df, beta, tau) {
  p_d <- sum(grepl("^D[0-9]+$", names(df)))
  stopifnot(length(beta) == p_d)
  Dmat <- as.matrix(df[ paste0("D", seq_len(p_d)) ])
  mD   <- as.matrix(df[ paste0("mD", seq_len(p_d)) ])
  rD   <- Dmat - mD
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
