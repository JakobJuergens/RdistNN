#' Bounded endogenous-treatment data for conditional LATE examples
#'
#' Generates a binary instrument/treatment model with known nuisance functions.
#' The latent compliance rank also affects untreated outcomes, so treatment is
#' endogenous conditional on the observed adjustment variables.
#'
#' @param n Positive integer number of observations.
#' @param seed Optional nonnegative integer seed, used locally without changing
#'   the caller's saved RNG state. `NULL` uses the current stream.
#' @param heterogeneous Whether effects also vary with an adjustment coordinate.
#'   This makes CLATE differ from CATE by 3/1024.
#' @param keep_latent Whether to return the simulated latent rank, potential
#'   treatments and potential outcomes for structural checks.
#'
#' @details X has two uniform reporting coordinates and Z has four independent
#'   uniform adjustment coordinates. With \eqn{S_j=\sin(2\pi Z_j)}, instrument
#'   propensity is \eqn{1/2+S_4/8}, always-taker probability is
#'   \eqn{a=1/4+S_1/32} and compliance probability is \eqn{c=1/3+S_2/32}.
#'   The independent latent rank H gives D(0)=1(H<=a), D(1)=1(H<=a+c).
#'   Outcomes are
#'   \deqn{Y(d)=g(X,Z)+\tau(X,Z)d+(H-1/2)/8+\varepsilon,}
#'   with equiprobable errors +/-1/8,
#'   \eqn{g=(S_1+S_4+X_1S_2+X_2S_3)/8} and
#'   \eqn{\tau=(X_1+X_2)/8}. The heterogeneous variant adds \eqn{S_2/16}
#'   to the treatment effect. The raw outcome is bounded by one in both models.
#'
#'   Core-model score variance is 13/(48 sqrt(15)); ratio local-limit variance
#'   is 39/(16 sqrt(15)). Those constants do not carry over to the heterogeneous
#'   variant: its sigma/omega metadata are explicitly NA.
#'
#' @return A list with data (observations and explicitly named oracle prediction
#'   columns), x0, truth (theta/cate/first_stage functions and core covariance),
#'   configuration and an optional latent audit frame. Truth functions accept a
#'   length-two vector or a matrix with two reporting columns in the unit square.
#'   Latent variables are excluded from the default adjustment columns.
#' @examples
#' z <- clate_demo_data(n = 50, seed = 1)
#' z$truth$theta(z$x0)
#' fit <- clate_dnn(
#'   data = z$data, x_cols = c("X1", "X2"), y_col = "Y", d_col = "D", a_col = "A",
#'   e_col = "e_oracle", m_y0_col = "m_y0_oracle", m_y1_col = "m_y1_oracle",
#'   m_d0_col = "m_d0_oracle", m_d1_col = "m_d1_oracle",
#'   x0 = z$x0, s = 8, training = "oracle"
#' )
#' summary(fit)
#' @export
clate_demo_data <- function(n = 200L, seed = NULL, heterogeneous = FALSE, keep_latent = FALSE) {
  n <- .affine_integer(n, "n")
  heterogeneous <- .validate_scalar_logical(heterogeneous, "heterogeneous")
  keep_latent <- .validate_scalar_logical(keep_latent, "keep_latent")
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != floor(seed)) {
      .abort_bad_arg("seed", "must be NULL or a nonnegative integer scalar")
    }
    seed <- as.integer(seed)
    restore_rng <- .rdistnn_rng_restore()
    on.exit(restore_rng(), add = TRUE)
    set.seed(seed)
  }
  X <- matrix(stats::runif(2L * n), n, 2L)
  Z <- matrix(stats::runif(4L * n), n, 4L)
  S <- sin(2 * pi * Z)
  e <- .5 + S[, 4L] / 8
  a <- .25 + S[, 1L] / 32
  c <- 1 / 3 + S[, 2L] / 32
  g <- (S[, 1L] + S[, 4L] + X[, 1L] * S[, 2L] + X[, 2L] * S[, 3L]) / 8
  tau <- rowSums(X) / 8 + if (heterogeneous) S[, 2L] / 16 else 0
  H <- stats::runif(n)
  epsilon <- sample(c(-1 / 8, 1 / 8), n, replace = TRUE)
  A <- stats::rbinom(n, 1L, e)
  D0 <- as.integer(H <= a)
  D1 <- as.integer(H <= a + c)
  D <- ifelse(A == 1L, D1, D0)
  Y0 <- g + (H - .5) / 8 + epsilon
  Y1 <- Y0 + tau
  Y <- Y0 + tau * D
  data <- data.frame(
    row_id = seq_len(n), X1 = X[, 1L], X2 = X[, 2L],
    Z1 = Z[, 1L], Z2 = Z[, 2L], Z3 = Z[, 3L], Z4 = Z[, 4L],
    A = A, D = D, Y = Y, e_oracle = e,
    m_y0_oracle = g + tau * a, m_y1_oracle = g + tau * (a + c),
    m_d0_oracle = a, m_d1_oracle = a + c
  )
  out <- list(
    data = data, x0 = c(.5, .5),
    truth = list(
      theta = .clate_effect_function(if (heterogeneous) 3 / 1024 else 0),
      cate = .clate_effect_function(0), first_stage = .clate_first_stage_function,
      sigma = if (heterogeneous) NA_real_ else 13 / (48 * sqrt(15)),
      omega = if (heterogeneous) NA_real_ else 39 / (16 * sqrt(15))
    ),
    configuration = list(n = n, seed = seed, heterogeneous = heterogeneous)
  )
  if (keep_latent) {
    out$latent <- data.frame(H = H, epsilon = epsilon, D0 = D0, D1 = D1, Y0 = Y0, Y1 = Y1)
  }
  out
}

.clate_truth_targets <- function(X) {
  if (is.null(dim(X))) X <- matrix(X, nrow = 1L)
  X <- .validate_matrix(X, "X", numeric = TRUE, finite = TRUE)
  if (ncol(X) != 2L || any(X < 0 | X > 1)) {
    .abort_bad_arg("X", "must have two reporting coordinates in [0, 1]")
  }
  X
}

# The factory captures the scalar offset, not the generator's observed/latent data.
.clate_effect_function <- function(offset) {
  force(offset)
  function(X) rowSums(.clate_truth_targets(X)) / 8 + offset
}

.clate_first_stage_function <- function(X) {
  rep(1 / 3, nrow(.clate_truth_targets(X)))
}
