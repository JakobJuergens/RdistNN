# Minimal nonparametric regression demo data (Y-only / singleton score)
# + ready-made multi-scenario generator for tests
# -------------------------------------------------------------------

#' Minimal nonparametric regression demo data (oracle m and Y only)
#'
#' Generates a simple nonparametric regression dataset with features \eqn{X}
#' and response \eqn{Y = m(X) + \epsilon}, returning the oracle regression
#' function \eqn{m(x)}. This is designed for experimenting with the Y-only
#' singleton score `score_np_y()`.
#'
#' @param n Integer; sample size.
#' @param p Integer; number of features in \eqn{X}.
#' @param design Character; one of \code{"linear"}, \code{"sinusoidal"}, \code{"bump"}, \code{"constant"}.
#' @param sigma Noise standard deviation in \eqn{Y}.
#' @param seed Optional integer seed, used locally without changing the caller's
#'   RNG state. `NULL` uses the current stream.
#' @param x0 Optional target point in \eqn{\mathbb{R}^p} (defaults to the zero vector).
#' @param inject Optional, one of \code{NULL}, \code{"na"}, \code{"inf"}; if provided,
#'   inject a single non-finite value into \code{Y} at \code{idx_inject} (useful for tests).
#' @param idx_inject Integer index at which to inject \code{NA}/\code{Inf} (default 3).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{data}: \code{data.frame} with columns \code{X1..Xp}, \code{Y}, and \code{m_oracle} (=\eqn{m(X)}).
#'   \item \code{x0}: the target point used.
#'   \item \code{truth}: list with closure \code{m_fn(x)} for evaluating \eqn{m(x)} at new \eqn{X}.
#'   \item \code{m_x0}: the true \eqn{m(x_0)}.
#' }
#'
#' @examples
#' demo <- np_demo_data(n = 200, p = 2, design = "sinusoidal", seed = 1)
#' head(demo$data)
#' demo$m_x0
#'
#' @export
np_demo_data <- function(n = 200,
                         p = 2,
                         design = c("linear", "sinusoidal", "bump", "constant"),
                         sigma = 0.5,
                         seed = NULL,
                         x0 = NULL,
                         inject = NULL,
                         idx_inject = 3L) {
  n <- .validate_scalar_int(n, arg = "n", min = 1L)
  p <- .validate_scalar_int(p, arg = "p", min = 1L)
  if (!is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma) || sigma < 0) {
    .abort_bad_arg("sigma", "must be a non-negative finite numeric scalar")
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

  # Oracle regression m(x)
  m_fn <- switch(design,
    linear = function(xx) {
      xx <- as.matrix(xx)
      0.5 * xx[, 1] + if (ncol(xx) >= 2) -0.25 * xx[, 2] else 0
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
    },
    constant = function(xx) {
      xx <- as.matrix(xx)
      rep(2, nrow(xx))
    }
  )

  mX <- m_fn(X)
  eps <- if (sigma > 0) rnorm(n, sd = sigma) else rep(0, n)
  Y <- mX + eps

  # Optional injection of non-finites (for guard tests)
  if (!is.null(inject)) {
    if (!inject %in% c("na", "inf")) {
      .abort_bad_arg("inject", "must be one of NULL, 'na', or 'inf'", quoted = FALSE)
    }
    idx_inject <- .validate_scalar_int(idx_inject, arg = "idx_inject", min = 1L)
    if (idx_inject > n) {
      .abort_bad_arg("idx_inject", sprintf("must be <= n (got idx_inject=%d, n=%d)", idx_inject, n))
    }
    if (identical(inject, "na")) Y[idx_inject] <- NA_real_
    if (identical(inject, "inf")) Y[idx_inject] <- Inf
  }

  df <- data.frame(X, Y = as.numeric(Y), m_oracle = as.numeric(mX))

  list(
    data = df,
    x0 = as.numeric(x0),
    truth = list(m_fn = function(xx) m_fn(as.matrix(xx))),
    m_x0 = as.numeric(m_fn(matrix(x0, nrow = 1)))
  )
}

#' Quick helper: sample a batch of subsets (singleton by default)
#' @inheritParams cate_demo_subsets
#' @inherit cate_demo_subsets return
#' @examples
#' draws <- np_demo_subsets(n = 10, K = 4, s = 1, seed = 1)
#' draws  # one singleton subset per row
#' @export
np_demo_subsets <- function(n, K = 1000L, s = 1L, seed = NULL) {
  .rdistnn_demo_subsets(n = n, K = K, s = s, seed = seed)
}
