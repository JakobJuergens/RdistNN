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
#' treatment effect. It also returns the R-learner–style pseudo-outcome
#' \deqn{
#'   \rho \;=\; \frac{(D - e)\,(Y - m)}{e\,(1 - e)} \, .
#' }
#'
#' No model fitting is performed; all quantities come from the DGP.
#'
#' @param n Integer sample size.
#' @param p Integer feature dimension (recommend \eqn{p \ge 2}).
#' @param design Character; one of \code{"simple"}, \code{"hetero"}, \code{"bump"} for \eqn{\tau(x)}.
#' @param sigma Noise standard deviation in \eqn{Y}.
#' @param clip Propensity clipping in \eqn{(0,\,0.5)}; applied to \eqn{e(x)}.
#' @param seed Optional integer seed for reproducibility.
#' @param x0 Optional numeric length-\eqn{p} target point (defaults to the zero vector).
#' @param inject Optional, one of \code{NULL}, \code{"na"}, \code{"inf"}; if provided,
#'   inject a single non-finite value into one of the returned columns (see \code{inject_col})
#'   at row \code{idx_inject}; useful for testing runtime guards.
#' @param inject_col Character; which column to corrupt when \code{inject} is not \code{NULL}.
#'   One of \code{"Y"}, \code{"rho_oracle"}, \code{"e_oracle"}, \code{"m_oracle"},
#'   \code{"m0_oracle"}, \code{"tau_oracle"}. Default \code{"rho_oracle"}.
#' @param idx_inject Integer index to corrupt (default 3).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{data}: a \code{data.frame} with columns \code{X1..Xp}, \code{D}, \code{Y},
#'         and oracle nuisances \code{e_oracle} (=\eqn{e(X)}),
#'         \code{m_oracle} (=\eqn{m(X)}),
#'         \code{m0_oracle} (=\eqn{m_{0}(X)}),
#'         \code{tau_oracle} (=\eqn{\tau(X)}),
#'         \code{rho_oracle} (=\eqn{\rho}).
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
                           inject_col = "rho_oracle",
                           idx_inject = 3L) {
  stopifnot(p >= 2L, clip > 0, clip < 0.5)
  design <- match.arg(design)
  if (!is.null(seed)) set.seed(seed)
  if (is.null(x0)) x0 <- rep(0, p)

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

  # R-learner pseudo-outcome (oracle nuisances)
  rho <- (D - e_true) * (Y - m_true) / (e_true * (1 - e_true))

  df <- data.frame(X,
    D = D, Y = as.numeric(Y),
    e_oracle = as.numeric(e_true),
    m_oracle = as.numeric(m_true),
    m0_oracle = as.numeric(m0_true),
    tau_oracle = as.numeric(tau_true),
    rho_oracle = as.numeric(rho)
  )

  # Optional injection of non-finites for guard testing
  if (!is.null(inject)) {
    stopifnot(idx_inject >= 1L, idx_inject <= n)
    stopifnot(inject_col %in% c("Y", "rho_oracle", "e_oracle", "m_oracle", "m0_oracle", "tau_oracle"))
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
#' replacement (all indices are distinct). Across rows, subsets are sampled
#' independently (the same index may reappear in different rows).
#'
#' @param n Integer; number of observations.
#' @param K Integer; number of subsets to sample.
#' @param s Integer; subset size (order of the kernel).
#' @param seed Optional integer seed for reproducibility.
#'
#' @return An integer matrix of size \eqn{K \times s}. Each row contains
#' \eqn{s} distinct integers in \eqn{\{1,\ldots,n\}} representing one subset.
#'
#' @examples
#' subs <- cate_demo_subsets(n = 200, K = 1000, s = 8, seed = 1)
#' dim(subs)
#'
#' @export
cate_demo_subsets <- function(n, K = 1000L, s = 8L, seed = NULL) {
  stopifnot(s >= 1L, s <= n, K >= 1L)
  if (!is.null(seed)) set.seed(seed)
  # Efficient in base R for moderate sizes
  t(replicate(K, sort(sample.int(n, s))))
}

#' Ready-made multi-scenario generator for CATE score tests
#'
#' Produces a named list of deterministic scenarios for unit tests of CATE score
#' factories and evaluators. Scenarios cover multiple \code{design}s and a set
#' of edge/failure cases (non-finite oracles, bad \code{x0} dimension, wrong
#' \code{subs} shape, empty \code{subs}).
#'
#' Each scenario is a list with:
#' \itemize{
#'   \item \code{data}: a \code{data.frame} containing \code{X1..Xp}, \code{D}, \code{Y},
#'         and oracle columns (\code{e_oracle}, \code{m_oracle}, \code{m0_oracle},
#'         \code{tau_oracle}, \code{rho_oracle}).
#'   \item \code{x0}: numeric length-\eqn{p} vector.
#'   \item \code{x_cols}: character names of the X columns (or integer indices).
#'   \item \code{s}, \code{k}: integers for score factories that accept these (uniform API).
#'   \item \code{subs_all}, \code{subs_empty}, \code{subs_bad_ncol}: subset matrices for batch eval tests.
#'   \item \code{cols}: a list with standard column names used by tests
#'         (\code{Y}, \code{D}, \code{e}, \code{m}, \code{m0}, \code{tau}, \code{rho}).
#' }
#' @param n Integer; number of rows (>= 6 recommended).
#' @param p Integer; number of X columns (>= 2).
#' @return Named list of scenarios.
#' @keywords internal
#' @noRd
cate_demo_sets <- function(n = 80L, p = 2L) {
  stopifnot(is.numeric(n), n >= 6L, is.numeric(p), p >= 2L)
  n <- as.integer(n)
  p <- as.integer(p)

  # Common subset matrices
  subs_all <- matrix(seq_len(n), ncol = 1L)
  subs_empty <- matrix(integer(0), ncol = 1L)
  subs_bad_ncol <- cbind(seq_len(n), rev(seq_len(n)))

  # Column map (standardized across tests)
  cols <- list(
    X   = paste0("X", seq_len(p)),
    Y   = "Y",
    D   = "D",
    e   = "e_oracle",
    m   = "m_oracle",
    m0  = "m0_oracle",
    tau = "tau_oracle",
    rho = "rho_oracle"
  )

  # Helper to assemble a scenario
  make_case <- function(dd, x0, x_cols, s = 1L, k = 1L) {
    list(
      data = dd,
      x0 = as.numeric(x0),
      x_cols = x_cols,
      s = as.integer(s),
      k = as.integer(k),
      subs_all = subs_all,
      subs_empty = subs_empty,
      subs_bad_ncol = subs_bad_ncol,
      cols = cols
    )
  }

  # --- SUCCESS CASES (seeded for determinism) ---
  ok_simple <- cate_demo_data(n = n, p = p, design = "simple", sigma = 0.4, clip = 1e-3, seed = 1)
  ok_hetero <- cate_demo_data(n = n, p = p, design = "hetero", sigma = 0.4, clip = 1e-3, seed = 2)
  ok_bump <- cate_demo_data(n = n, p = p, design = "bump", sigma = 0.4, clip = 1e-3, seed = 3)
  ok_clean <- cate_demo_data(n = n, p = p, design = "simple", sigma = 0.0, clip = 1e-3, seed = 4) # noiseless Y

  # --- EDGE / FAILURE CASES (non-finites injected) ---
  nf_rho_na <- cate_demo_data(
    n = n, p = p, design = "hetero", sigma = 0.3, clip = 1e-3, seed = 5,
    inject = "na", inject_col = "rho_oracle", idx_inject = 3L
  )
  nf_rho_inf <- cate_demo_data(
    n = n, p = p, design = "hetero", sigma = 0.3, clip = 1e-3, seed = 6,
    inject = "inf", inject_col = "rho_oracle", idx_inject = 2L
  )
  nf_Y_na <- cate_demo_data(
    n = n, p = p, design = "simple", sigma = 0.3, clip = 1e-3, seed = 7,
    inject = "na", inject_col = "Y", idx_inject = 4L
  )
  nf_e_inf <- cate_demo_data(
    n = n, p = p, design = "simple", sigma = 0.3, clip = 1e-3, seed = 8,
    inject = "inf", inject_col = "e_oracle", idx_inject = 5L
  )

  # Bad x0 dimension
  x0_ok <- rep(0, p)
  x0_bad <- rep(0, p + 1L)

  list(
    # ---- OK ----
    ok_simple = make_case(ok_simple$data, x0_ok, cols$X, s = 1L, k = 1L),
    ok_hetero = make_case(ok_hetero$data, x0_ok, cols$X, s = 1L, k = 1L),
    ok_bump = make_case(ok_bump$data, x0_ok, cols$X, s = 1L, k = 1L),
    ok_clean = make_case(ok_clean$data, x0_ok, cols$X, s = 1L, k = 1L),

    # ---- Non-finites (expect_error at eval/setup depending on your guards) ----
    nonfinite_rho_na = make_case(nf_rho_na$data, x0_ok, cols$X),
    nonfinite_rho_inf = make_case(nf_rho_inf$data, x0_ok, cols$X),
    nonfinite_Y_na = make_case(nf_Y_na$data, x0_ok, cols$X),
    nonfinite_e_inf = make_case(nf_e_inf$data, x0_ok, cols$X),

    # ---- Dimension / shape guards ----
    bad_x0_dim = make_case(ok_simple$data, x0_bad, cols$X),
    bad_subs_ncol = make_case(ok_simple$data, x0_ok, cols$X)
  )
}
