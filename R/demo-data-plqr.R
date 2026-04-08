# PLQR demo data (oracle nuisances; no internal fitting)
# + test-friendly injection hooks
# + multi-scenario generator for unit tests
# -------------------------------------------------------------------

# --- Internal: Asymmetric Laplace draws (location 0, scale = sigma, asymmetry = tau) ----
# Ensures Q_tau(U) = 0, i.e., the tau-quantile of U is zero.
.r_ald <- function(n, tau = 0.5, sigma = 1) {
  n <- .validate_scalar_int(n, arg = "n", min = 1L)
  tau <- .validate_scalar_prob(tau, arg = "tau")
  if (!is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
    .abort_bad_arg("sigma", "must be a positive finite numeric scalar")
  }
  sigma <- as.numeric(sigma)
  u <- runif(n)
  ifelse(
    u < tau,
    sigma * ((1 - tau) * log(u / tau)),
    -sigma * (tau * log((1 - u) / (1 - tau)))
  )
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
#' @param inject Optional, one of \code{NULL}, \code{"na"}, \code{"inf"}; if provided,
#'   inject a single non-finite into \code{inject_col} at row \code{idx_inject}.
#'   Useful for testing runtime guards.
#' @param inject_col Character; which column to corrupt when \code{inject} is not \code{NULL}.
#'   One of \code{"Y"}, \code{"q_oracle"}, \code{"D1"}, \code{"mD1"}, \code{"rD1"} (extend as needed).
#' @param idx_inject Integer row index to corrupt (default 3).
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
                           x0 = NULL,
                           inject = NULL,
                           inject_col = c("Y", "q_oracle", "D1", "mD1", "rD1")[1],
                           idx_inject = 3L) {
  n <- .validate_scalar_int(n, arg = "n", min = 1L)
  p <- .validate_scalar_int(p, arg = "p", min = 1L)
  p_d <- .validate_scalar_int(p_d, arg = "p_d", min = 1L)
  tau <- .validate_scalar_prob(tau, arg = "tau")
  if (!is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
    .abort_bad_arg("sigma", "must be a positive finite numeric scalar")
  }
  sigma <- as.numeric(sigma)
  design <- match.arg(design)
  if (!is.null(seed)) set.seed(seed)
  if (is.null(x0)) {
    x0 <- rep(0, p)
  } else {
    x0 <- .validate_target_vector(x0, p = p, x0_arg = "x0", x_arg = "X", finite = TRUE)
  }

  # Covariates
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", seq_len(p))

  # True beta (partially linear effect at quantile tau)
  beta_true <- rep(0.5, p_d)

  # Nonparametric quantile component q_tau(X)
  q_fn <- switch(design,
    linear = function(xx) {
      xx <- as.matrix(xx)
      0.5 * xx[, 1]
    },
    sinusoidal = function(xx) {
      xx <- as.matrix(xx)
      0.5 * sin(xx[, 1]) + if (ncol(xx) >= 2) 0.25 * xx[, 2]^2 else 0
    },
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

  # Partially linear regressors D with oracle E(D | X) = mD(X):
  # simple linear mD with different loadings per column of D
  Gamma <- matrix(0, nrow = p, ncol = p_d)
  for (j in seq_len(p_d)) {
    g <- c(0.6, -0.4, rep(0, max(0, p - 2L))) + (j - 1) * 0.2
    Gamma[, j] <- g[seq_len(p)]
  }
  mD <- X %*% Gamma # n x p_d
  V <- matrix(rnorm(n * p_d, sd = 0.6), n, p_d)
  D <- mD + V
  colnames(D) <- paste0("D", seq_len(p_d))

  # Error with tau-quantile zero
  U <- .r_ald(n, tau = tau, sigma = sigma)

  # Outcome: Y = D^T beta_true + q_tau(X) + U  with Q_tau(U | X, D) = 0
  Y <- as.vector(D %*% beta_true) + qX + U

  # Tidy outputs
  df <- data.frame(X, D, Y = as.numeric(Y), q_oracle = as.numeric(qX))
  for (j in seq_len(p_d)) {
    df[[paste0("mD", j)]] <- mD[, j]
    df[[paste0("rD", j)]] <- D[, j] - mD[, j]
  }

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
      "Y", "q_oracle",
      paste0("D", seq_len(p_d)[1]),
      paste0("mD", seq_len(p_d)[1]),
      paste0("rD", seq_len(p_d)[1])
    )
    if (!inject_col %in% allowed_inject_cols) {
      stop(sprintf(
        "Column '%s' (specified by `inject_col`) is not a valid injection target.\nValid options: %s",
        inject_col, paste(allowed_inject_cols, collapse = ", ")
      ), call. = FALSE)
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
      beta_true = as.numeric(beta_true),
      q_fn = function(xx) q_fn(as.matrix(xx)),
      mD_fn = mD_fn
    ),
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
#'   \code{mD1..mDp_d}, and \code{q_oracle}.
#' @param beta Numeric length-\code{p_d} coefficient vector \eqn{\beta}.
#' @param tau Quantile level \eqn{\tau} (should match the data-generating level).
#'
#' @return A numeric matrix of size \eqn{n \times p_d}.
#' @examples
#' demo <- plqr_demo_data(n = 100, p = 2, p_d = 2, tau = 0.5, seed = 1)
#' plqr_moment(demo$data, beta = demo$truth$beta_true, tau = 0.5)
#' @export
plqr_moment <- function(df, beta, tau) {
  if (!is.data.frame(df)) {
    .abort_bad_arg("df", "must be a data.frame")
  }
  p_d <- sum(grepl("^D[0-9]+$", names(df)))
  if (!is.numeric(beta)) {
    .abort_bad_arg("beta", "must be numeric")
  }
  if (.has_na_only(beta)) {
    .abort_bad_arg("beta", "cannot contain `NA` values")
  }
  if (.has_nonfinite_non_na(beta)) {
    .abort_bad_arg("beta", "cannot contain non-finite values (`Inf`, `-Inf`, or `NaN`)")
  }
  if (length(beta) != p_d) {
    .abort_bad_arg(
      "beta",
      sprintf(
        "length must equal number of D-columns (got length(beta)=%d, p_d=%d)",
        length(beta), p_d
      )
    )
  }
  tau <- .validate_scalar_prob(tau, arg = "tau")
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
#' @examples
#' plqr_demo_subsets(n = 10, K = 4, s = 3, seed = 1)
#' @export
plqr_demo_subsets <- function(n, K = 1000L, s = 8L, seed = NULL) {
  n <- .validate_scalar_int(n, arg = "n", min = 1L)
  s <- .validate_scalar_int(s, arg = "s", min = 1L)
  K <- .validate_scalar_int(K, arg = "K", min = 1L)
  if (s > n) {
    .abort_bad_arg("s", sprintf("must satisfy s <= n (got s=%d, n=%d)", s, n))
  }
  if (!is.null(seed)) set.seed(seed)
  t(replicate(K, sort(sample.int(n, s))))
}

#' Ready-made multi-scenario generator for PLQR tests
#'
#' Produces a named list of scenarios for unit tests of PLQR scores/moments.
#' Scenarios include multiple designs and edge/failure cases (non-finite Y/q/mD/rD),
#' plus dimension/shape guards.
#'
#' Each scenario is a list with:
#' \itemize{
#'   \item \code{data}: \code{data.frame} with \code{X1..Xp}, \code{D1..Dp_d}, \code{Y},
#'         \code{q_oracle}, and \code{mD*}/\code{rD*}.
#'   \item \code{x0}: numeric length-\eqn{p}.
#'   \item \code{x_cols}: character names of X columns (or integer indices).
#'   \item \code{d_cols}: character names of D columns.
#'   \item \code{tau}: quantile level used in the DGP.
#'   \item \code{beta_true}: vector used in the DGP (length \eqn{p_d}).
#'   \item \code{s}, \code{k}: integers for uniform API.
#'   \item \code{subs_all}, \code{subs_empty}, \code{subs_bad_ncol}: subset matrices.
#'   \item \code{cols}: a list mapping standard column aliases: \code{Y}, \code{q}, \code{mD}, \code{rD}.
#' }
#'
#' @param n Integer; number of rows (>= 10 recommended).
#' @param p Integer; number of X columns (>= 1).
#' @param p_d Integer; number of D columns (>= 1).
#' @return Named list of scenarios.
#' @keywords internal
#' @noRd
plqr_demo_sets <- function(n = 120L, p = 2L, p_d = 2L) {
  n <- .validate_scalar_int(n, arg = "n", min = 10L)
  p <- .validate_scalar_int(p, arg = "p", min = 1L)
  p_d <- .validate_scalar_int(p_d, arg = "p_d", min = 1L)

  # Subset matrices
  subs_all <- matrix(seq_len(n), ncol = 1L)
  subs_empty <- matrix(integer(0), ncol = 1L)
  subs_bad_ncol <- cbind(seq_len(n), rev(seq_len(n)))

  X_names <- paste0("X", seq_len(p))
  D_names <- paste0("D", seq_len(p_d))
  mD_names <- paste0("mD", seq_len(p_d))
  rD_names <- paste0("rD", seq_len(p_d))

  cols <- list(Y = "Y", q = "q_oracle", mD = mD_names, rD = rD_names, D = D_names)

  make_case <- function(obj) {
    list(
      data = obj$data,
      x0 = obj$x0,
      x_cols = X_names,
      d_cols = D_names,
      tau = attr(obj, "tau"),
      beta_true = obj$truth$beta_true,
      s = 1L, k = 1L,
      subs_all = subs_all,
      subs_empty = subs_empty,
      subs_bad_ncol = subs_bad_ncol,
      cols = cols
    )
  }

  # --- SUCCESS CASES (seeded; include a noiseless design for exact checks) ---
  ok_lin <- plqr_demo_data(n, p, p_d, tau = 0.5, sigma = 0.5, design = "linear", seed = 1)
  attr(ok_lin, "tau") <- 0.5
  ok_sin <- plqr_demo_data(n, p, p_d, tau = 0.5, sigma = 0.4, design = "sinusoidal", seed = 2)
  attr(ok_sin, "tau") <- 0.5
  ok_bump <- plqr_demo_data(n, p, p_d, tau = 0.5, sigma = 0.6, design = "bump", seed = 3)
  attr(ok_bump, "tau") <- 0.5
  ok_clean <- plqr_demo_data(n, p, p_d, tau = 0.5, sigma = 1e-8, design = "linear", seed = 4)
  attr(ok_clean, "tau") <- 0.5

  # --- EDGE / FAILURE CASES (non-finites injected) ---
  nf_Y_na <- plqr_demo_data(n, p, p_d,
    tau = 0.5, sigma = 0.5, design = "linear",
    seed = 5, inject = "na", inject_col = "Y", idx_inject = 3L
  )
  attr(nf_Y_na, "tau") <- 0.5
  nf_q_inf <- plqr_demo_data(n, p, p_d,
    tau = 0.5, sigma = 0.5, design = "sinusoidal",
    seed = 6, inject = "inf", inject_col = "q_oracle", idx_inject = 2L
  )
  attr(nf_q_inf, "tau") <- 0.5
  nf_mD_na <- plqr_demo_data(n, p, p_d,
    tau = 0.5, sigma = 0.5, design = "bump",
    seed = 7, inject = "na", inject_col = "mD1", idx_inject = 4L
  )
  attr(nf_mD_na, "tau") <- 0.5
  nf_rD_inf <- plqr_demo_data(n, p, p_d,
    tau = 0.5, sigma = 0.5, design = "linear",
    seed = 8, inject = "inf", inject_col = "rD1", idx_inject = 5L
  )
  attr(nf_rD_inf, "tau") <- 0.5

  # Bad x0 dimension scenario (use ok_lin data)
  bad_x0 <- ok_lin
  bad_x0$x0 <- rep(0, p + 1L)

  list(
    ok_linear = make_case(ok_lin),
    ok_sinusoid = make_case(ok_sin),
    ok_bump = make_case(ok_bump),
    ok_clean = make_case(ok_clean),
    nonfinite_Y_na = make_case(nf_Y_na),
    nonfinite_q_inf = make_case(nf_q_inf),
    nonfinite_mD_na = make_case(nf_mD_na),
    nonfinite_rD_inf = make_case(nf_rD_inf),
    bad_x0_dim = list(
      data = bad_x0$data,
      x0 = bad_x0$x0,
      x_cols = X_names,
      d_cols = D_names,
      tau = attr(bad_x0, "tau"),
      beta_true = bad_x0$truth$beta_true,
      s = 1L, k = 1L,
      subs_all = subs_all,
      subs_empty = subs_empty,
      subs_bad_ncol = subs_bad_ncol,
      cols = cols
    ),
    bad_subs_ncol = make_case(ok_lin)
  )
}
