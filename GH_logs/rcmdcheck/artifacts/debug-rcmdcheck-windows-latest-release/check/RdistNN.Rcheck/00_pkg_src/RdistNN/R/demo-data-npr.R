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
#' @param seed Optional integer seed for reproducibility.
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
  stopifnot(n >= 1L, p >= 1L, sigma >= 0)
  design <- match.arg(design)
  if (!is.null(seed)) set.seed(seed)
  if (is.null(x0)) x0 <- rep(0, p)

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
      center <- c(0.5, -0.5, rep(0, max(0, p - 2L)))
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
    stopifnot(idx_inject >= 1L, idx_inject <= n)
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
#' @examples
#' np_demo_subsets(n = 10, K = 4, s = 1, seed = 1)
#' @export
np_demo_subsets <- function(n, K = 1000L, s = 1L, seed = NULL) {
  stopifnot(s >= 1L, s <= n, K >= 1L)
  if (!is.null(seed)) set.seed(seed)
  t(replicate(K, sort(sample.int(n, s))))
}

#' Ready-made multi-scenario generator for `score_np_y` tests
#'
#' Produces a named list of scenarios used in unit tests of the Y-only singleton
#' score. Scenarios include successful cases (column-based Y, external Y, integer
#' \code{x_cols}, constant Y, integer Y) as well as failure/edge cases (non-finite
#' Y, bad \code{s}, bad \code{k}, dimension mismatch for \code{x0}, and a bad
#' \code{subs} matrix with wrong number of columns).
#'
#' @param n Integer; number of rows (>= 6 recommended).
#' @param p Integer; number of X columns (>= 1).
#' @return Named list of scenarios. Each scenario is a list with:
#'   \code{data}, \code{x0}, \code{x_cols}, \code{s}, \code{k}, \code{y_col}, \code{y},
#'   \code{subs_all}, \code{subs_empty}, \code{subs_bad_ncol}.
#' @keywords internal
#' @noRd
np_y_demo_sets <- function(n = 8L, p = 3L) {
  stopifnot(is.numeric(n), n >= 6L, is.numeric(p), p >= 1L)
  n <- as.integer(n)
  p <- as.integer(p)

  # Deterministic X for exact-equality expectations
  X <- as.data.frame(
    setNames(
      lapply(seq_len(p), function(j) j + 10 * j + seq_len(n)),
      paste0("X", seq_len(p))
    )
  )

  # A simple deterministic m(X) and Y (no noise)
  coeff <- (-1)^(seq_len(p)) * seq_len(p)
  mX <- as.numeric(as.matrix(X) %*% coeff) + seq_len(n)
  Y_const <- rep(5, n)
  Y_int <- as.integer(round(mX))

  # Subset matrices for batch eval tests
  subs_all <- matrix(seq_len(n), ncol = 1L)
  subs_empty <- matrix(integer(0), ncol = 1L)
  subs_bad_ncol <- cbind(seq_len(n), rev(seq_len(n)))

  # Helpers
  x0_ok <- rep(0, p)
  x0_bad <- rep(0, p + 1L)

  make_case <- function(data, x0, x_cols, s = 1L, k = 1L, y_col = "Y", y = NULL) {
    list(
      data = data,
      x0 = x0,
      x_cols = x_cols,
      s = as.integer(s),
      k = as.integer(k),
      y_col = y_col,
      y = if (!is.null(y)) as.numeric(y) else NULL,
      subs_all = subs_all,
      subs_empty = subs_empty,
      subs_bad_ncol = subs_bad_ncol
    )
  }

  list(
    # --- SUCCESS CASES ---
    ok_ycol = make_case(
      data = transform(X, Y = mX),
      x0 = x0_ok,
      x_cols = names(X),
      s = 1L, k = 1L,
      y_col = "Y", y = NULL
    ),
    ok_y_external = make_case(
      data = X,
      x0 = x0_ok,
      x_cols = names(X),
      s = 1L, k = 1L,
      y_col = NULL, y = mX
    ),
    ok_xcols_int = make_case(
      data = transform(X, Y = mX),
      x0 = x0_ok,
      x_cols = seq_len(p),
      s = 1L, k = 1L,
      y_col = "Y", y = NULL
    ),
    ok_constant_Y = make_case(
      data = transform(X, Y = Y_const),
      x0 = x0_ok,
      x_cols = names(X),
      s = 1L, k = 1L,
      y_col = "Y", y = NULL
    ),
    ok_integer_Y = make_case(
      data = transform(X, Y = Y_int),
      x0 = x0_ok,
      x_cols = names(X),
      s = 1L, k = 1L,
      y_col = "Y", y = NULL
    ),

    # --- EDGE / FAILURE CASES ---
    nonfinite_na = make_case(
      data = transform(X, Y = replace(mX, 3L, NA_real_)),
      x0 = x0_ok,
      x_cols = names(X),
      s = 1L, k = 1L,
      y_col = "Y", y = NULL
    ),
    nonfinite_inf = make_case(
      data = transform(X, Y = replace(mX, 2L, Inf)),
      x0 = x0_ok,
      x_cols = names(X),
      s = 1L, k = 1L,
      y_col = "Y", y = NULL
    ),
    bad_s_not1 = make_case(
      data = transform(X, Y = mX),
      x0 = x0_ok,
      x_cols = names(X),
      s = 2L, k = 1L,
      y_col = "Y", y = NULL
    ),
    bad_k_nonpos = make_case(
      data = transform(X, Y = mX),
      x0 = x0_ok,
      x_cols = names(X),
      s = 1L, k = 0L,
      y_col = "Y", y = NULL
    ),
    bad_x0_dim = make_case(
      data = transform(X, Y = mX),
      x0 = x0_bad,
      x_cols = names(X),
      s = 1L, k = 1L,
      y_col = "Y", y = NULL
    ),
    bad_subs_ncol = make_case(
      data = transform(X, Y = mX),
      x0 = x0_ok,
      x_cols = names(X),
      s = 1L, k = 1L,
      y_col = "Y", y = NULL
    )
  )
}
