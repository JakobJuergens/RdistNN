# R/demo-data-cate.R

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
                           x0 = NULL) {
  stopifnot(p >= 2L)
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
    simple = function(xx) 0.5 + 0 * xx[, 1],
    hetero = function(xx) 0.5 + 0.5 * sin(xx[, 1]),
    bump = {
      # bump near c(0.5, -0.5), flat elsewhere
      center <- c(0.5, -0.5, rep(0, max(0, ncol(xx) - 2L)))
      function(xx) {
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
    D = D, Y = Y,
    e_oracle = e_true,
    m_oracle = m_true,
    m0_oracle = m0_true,
    tau_oracle = tau_true,
    rho_oracle = rho
  )

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
