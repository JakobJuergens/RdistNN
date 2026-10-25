# Minimal CATE demo data (oracle nuisances + pseudo-outcome)
# + test-friendly injection hooks
# + multi-scenario generator for unit tests
# -------------------------------------------------------------------

#' Minimal CATE demo data (oracle nuisances + pseudo-outcome)
#'
#' Generates a small synthetic dataset for experimenting with the CATE score.
#' It returns \eqn{X}, \eqn{D}, \eqn{Y}, along with oracle nuisance quantities:
#' \deqn{
#'   m(x) \;=\; \mathbb{E}[Y \mid X=x]
#'          \;=\; m_{0}(x) \;+\; e(x)\,\tau(x) \, , \qquad
#'   e(x) \;=\; \Pr(D=1 \mid X=x) \, ,
#' }
#' where \eqn{m_{0}(x) = \mathbb{E}[Y \mid X=x, D=0]} and \eqn{\tau(x)} is the
#' treatment effect. It also returns the oracle AIPW pseudo-outcome.
#'
#' No model fitting is performed; all quantities come from the DGP.
#'
#' @param n Integer sample size.
#' @param p Integer feature dimension, at least 2.
#' @param design Character; one of \code{"simple"}, \code{"hetero"}, \code{"bump"} for \eqn{\tau(x)}.
#' @param sigma Noise standard deviation in \eqn{Y}.
#' @param clip Propensity clipping in \eqn{(0,\,0.5)}; applied to \eqn{e(x)}.
#' @param seed Optional integer seed, used locally without changing the caller's
#'   RNG state. `NULL` uses the current stream.
#' @param x0 Optional numeric length-\eqn{p} target point (defaults to the zero vector).
#' @param inject Optional, one of \code{NULL}, \code{"na"}, \code{"inf"}; if provided,
#'   inject a single non-finite value into one of the returned columns (see \code{inject_col})
#'   at row \code{idx_inject}; useful for testing runtime guards.
#' @param inject_col Character; which column to corrupt when \code{inject} is not \code{NULL}.
#'   One of \code{"Y"}, \code{"aipw_oracle"}, \code{"e_oracle"}, \code{"m_oracle"},
#'   \code{"m0_oracle"}, \code{"tau_oracle"}. Default \code{"aipw_oracle"}.
#' @param idx_inject Integer index to corrupt (default 3).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{data}: a \code{data.frame} with columns \code{X1..Xp}, \code{D}, \code{Y},
#'         and oracle nuisances \code{e_oracle} (=\eqn{e(X)}),
#'         \code{m_oracle} (=\eqn{m(X)}),
#'         \code{m0_oracle} (=\eqn{m_{0}(X)}),
#'         \code{m1_oracle} (treated outcome mean),
#'         \code{tau_oracle} (=\eqn{\tau(X)}),
#'         \code{aipw_oracle} (oracle AIPW signal).
#'   \item \code{x0}: the target point used (length \eqn{p}).
#'   \item \code{truth}: list of closures \code{tau_fn(x)}, \code{e_fn(x)},
#'         \code{m0_fn(x)}, \code{m_fn(x)} to evaluate ground truth at any \eqn{x}.
#'   \item \code{tau_x0}: the true \eqn{\tau(x_{0})}.
#' }
#'
#' @examples
#' demo <- cate_demo_data(n = 200, p = 2, design = "hetero", seed = 1)
#' str(demo$data)
#' demo$tau_x0
#'
#' @export
cate_demo_data <- function(n = 200,
                           p = 2,
                           design = c("simple", "hetero", "bump"),
                           sigma = 0.5,
                           clip = 1e-3,
                           seed = NULL,
                           x0 = NULL,
                           inject = NULL,
                           inject_col = "aipw_oracle",
                           idx_inject = 3L) {
  n <- .validate_scalar_int(n, arg = "n", min = 1L)
  p <- .validate_scalar_int(p, arg = "p", min = 2L)
  if (!is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma) || sigma < 0) {
    .abort_bad_arg("sigma", "must be a non-negative finite numeric scalar")
  }
  sigma <- as.numeric(sigma)
  if (!is.numeric(clip) || length(clip) != 1L || !is.finite(clip) || clip <= 0 || clip >= 0.5) {
    .abort_bad_arg("clip", "must be a finite numeric scalar in (0, 0.5)")
  }
  clip <- as.numeric(clip)
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

  # Draw covariates
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", seq_len(p))

  # Propensity: e(x) = logistic(beta^T x), clip to [clip, 1-clip]
  beta <- c(0.5, -0.5, rep(0, max(0, p - 2L)))
  lin_e <- as.vector(X %*% beta[seq_len(p)])
  e_true <- plogis(lin_e)
  e_true <- pmin(pmax(e_true, clip), 1 - clip)

  # Baseline outcome (under control)
  m0_true <- 1 + X[, 1] - X[, 2]

  # Treatment effect tau(x)
  tau_fn <- switch(design,
    simple = function(xx) {
      xx <- as.matrix(xx)
      0.5 + 0 * xx[, 1]
    },
    hetero = function(xx) {
      xx <- as.matrix(xx)
      0.5 + 0.5 * sin(xx[, 1])
    },
    bump = {
      # bump near c(0.5, -0.5), flat elsewhere
      function(xx) {
        xx <- as.matrix(xx)
        center <- c(0.5, -0.5, rep(0, max(0, ncol(xx) - 2L)))
        d2 <- rowSums((xx - rep(center, each = nrow(xx)))^2)
        0.2 + 0.8 * exp(-0.5 * d2) # bump width ~ 1
      }
    }
  )
  tau_true <- tau_fn(X)

  # Treatment and outcome
  D <- rbinom(n, size = 1, prob = e_true)
  eps <- rnorm(n, sd = sigma)
  Y <- m0_true + D * tau_true + eps

  # Unconditional regression function m(x) = E[Y|X=x]
  m_true <- m0_true + e_true * tau_true

  # AIPW pseudo-outcome (oracle nuisances)
  m1_true <- m0_true + tau_true
  rho <- cate_rho_aipw(Y, D, m1_true, m0_true, e_true, clip = clip)

  df <- data.frame(X,
    D = D, Y = as.numeric(Y),
    e_oracle = as.numeric(e_true),
    m_oracle = as.numeric(m_true),
    m0_oracle = as.numeric(m0_true),
    m1_oracle = as.numeric(m1_true),
    tau_oracle = as.numeric(tau_true),
    aipw_oracle = as.numeric(rho)
  )

  # Optional injection of non-finites for guard testing
  if (!is.null(inject)) {
    if (!inject %in% c("na", "inf")) {
      .abort_bad_arg("inject", "must be one of NULL, 'na', or 'inf'", quoted = FALSE)
    }
    idx_inject <- .validate_scalar_int(idx_inject, arg = "idx_inject", min = 1L)
    if (idx_inject > n) {
      .abort_bad_arg("idx_inject", sprintf("must be <= n (got idx_inject=%d, n=%d)", idx_inject, n))
    }
    allowed_inject_cols <- c("Y", "aipw_oracle", "e_oracle", "m_oracle", "m0_oracle", "tau_oracle")
    if (!inject_col %in% allowed_inject_cols) {
      .dnn_abort("argument_invalid", sprintf(
        "Column '%s' (specified by `inject_col`) is not a valid injection target.\nValid options: %s",
        inject_col, paste(allowed_inject_cols, collapse = ", ")
      ), argument = "inject_col")
    }
    if (identical(inject, "na")) df[[inject_col]][idx_inject] <- NA_real_
    if (identical(inject, "inf")) df[[inject_col]][idx_inject] <- Inf
  }

  # Ground-truth closures for convenience
  e_fn <- function(xx) {
    xx <- as.matrix(xx)
    lin <- as.vector(xx %*% beta[seq_len(ncol(xx))])
    ee <- plogis(lin)
    pmin(pmax(ee, clip), 1 - clip)
  }
  m0_fn <- function(xx) {
    xx <- as.matrix(xx)
    1 + xx[, 1] - xx[, 2]
  }
  m_fn <- function(xx) {
    xx <- as.matrix(xx)
    m0_fn(xx) + e_fn(xx) * tau_fn(xx)
  }

  list(
    data = df,
    x0 = as.numeric(x0),
    truth = list(
      tau_fn = function(xx) tau_fn(as.matrix(xx)),
      e_fn = e_fn, m0_fn = m0_fn, m_fn = m_fn
    ),
    tau_x0 = as.numeric(tau_fn(matrix(x0, nrow = 1)))
  )
}

#' Quick helper: sample a batch of subsets (for U-stat experiments)
#'
#' Draws \eqn{K} subsets of size \eqn{s} uniformly at random from
#' \eqn{\{1,\ldots,n\}}. Within each subset, elements are sampled **without**
#' replacement (all indices are distinct). Across draws, subsets are sampled
#' independently (the same index may reappear in different rows).
#'
#' @param n Integer; number of observations.
#' @param K Integer; number of subsets to sample.
#' @param s Integer; subset size (order of the kernel).
#' @param seed Optional integer seed, used locally without changing the caller's
#'   RNG state. `NULL` uses the current stream.
#'
#' @return An integer `K x s` matrix with one sorted subset per row, including
#'   when `K = 1` or `s = 1`. The result can be passed to [score_eval_batch()].
#'
#' @examples
#' subs <- cate_demo_subsets(n = 200, K = 1000, s = 8, seed = 1)
#' dim(subs)
#'
#' @export
cate_demo_subsets <- function(n, K = 1000L, s = 8L, seed = NULL) {
  .rdistnn_demo_subsets(n = n, K = K, s = s, seed = seed)
}
