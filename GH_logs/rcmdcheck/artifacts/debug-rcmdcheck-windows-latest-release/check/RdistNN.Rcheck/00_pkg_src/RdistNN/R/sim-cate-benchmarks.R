# sim-cate-benchmarks.R
# Benchmark estimator implementations for the CATE simulation study.
#
# All methods share a common internal interface:
#
#   .method_XXX(dgp, x_cols, k, s, alpha, ...)
#     dgp    -- list returned by sim_cate_data()
#     x_cols -- character vector of covariate column names
#     k, s   -- DNN tuning (only used by dnn_* methods)
#     alpha  -- nominal coverage level
#     ...    -- method-specific extra arguments (absorbed silently if unused)
#
#   Returns a named numeric vector:
#     tau_hat, tau0, bias, sq_err, se, cover, ci_len, neff,
#     bw_eff, n_folds, n_per_fold
#
#   tau_hat / tau0   -- raw estimate and truth; self-contained for downstream
#   neff             -- neff_x = 1/sum(w_i^2(x)); paper eq. (02b L.81); NA for
#                       global methods
#   bw_eff           -- empirical localization radius r_s(x): distance to the
#                       s-th nearest neighbour for DNN (paper 02b L.71-72;
#                       r_s(x) ≍ h_eff = s^{-1/k_x}); mean Silverman h for
#                       local_constant; NA for global methods. Informative about
#                       actual local density realizations across reps.
#   n_folds          -- cross-fitting folds; non-NA only for dnn_crossfit
#   n_per_fold       -- obs per fold (n %/% n_folds); non-NA only for dnn_crossfit
#
# Methods are registered in .SIM_CATE_METHODS and dispatched through
# .get_sim_method().  The public API is sim_cate_methods().
#
# GRF-based methods (causal_forest, llf) require the grf package listed
# under Suggests.  They return NA metrics with a warning when grf is absent.

# ── Shared metric constructor ─────────────────────────────────────────────────

.sim_metric <- function(tau_hat, se_hat, tau0, alpha, neff = NA_real_) {
  z       <- stats::qnorm(1 - alpha / 2)
  finite  <- is.finite(tau_hat) && is.finite(se_hat) && se_hat > 0
  covered <- if (finite) {
    (tau_hat - z * se_hat <= tau0) & (tau0 <= tau_hat + z * se_hat)
  } else {
    NA_real_
  }
  c(tau_hat = as.numeric(tau_hat),
    tau0    = as.numeric(tau0),
    bias    = as.numeric(tau_hat - tau0),
    sq_err  = as.numeric((tau_hat - tau0)^2),
    se      = as.numeric(se_hat),
    cover   = as.numeric(covered),
    ci_len  = as.numeric(2 * z * se_hat),
    neff    = as.numeric(neff))
}

# ── DNN oracle ────────────────────────────────────────────────────────────────
# Proposed method.  Uses oracle nuisances so that first-stage error does not
# contaminate the comparison.

.method_dnn_oracle <- function(dgp, x_cols, k, s, alpha, ...) {
  fit <- cate_dnn(
    data      = dgp$data,
    x_cols    = x_cols,
    y_col     = "Y",
    t_col     = "D",
    e_col     = "e_oracle",
    m0_col    = "m0_oracle",
    m1_col    = "m1_oracle",
    x0        = dgp$x0,
    k         = k,
    s         = s,
    se        = TRUE,
    jackknife = "fast",
    return    = "ess"
  )
  .sim_metric(fit$tau, fit$se, dgp$tau_x0, alpha, neff = fit$ess)
}

# ── DNN cross-fitted ──────────────────────────────────────────────────────────
# Proposed method with K-fold cross-fitted nuisances estimated via
# logistic regression (propensity) and OLS (outcome mean), both base-R.
# Compares same-sample inference against the cross-fitting alternative.

.method_dnn_crossfit <- function(dgp, x_cols, k, s, alpha,
                                  n_folds = 5L, clip = 1e-3, ...) {
  dat <- dgp$data
  n   <- nrow(dat)

  # Balanced fold assignment; uses the RNG state continuing from the DGP draw,
  # so the result is deterministic for a given rep seed.
  fold_ids <- sample(rep_len(seq_len(n_folds), n))

  e_cf <- numeric(n)
  m_cf <- numeric(n)

  e_fmla <- stats::as.formula(paste("D ~", paste(x_cols, collapse = " + ")))
  # Outcome model: E[Y | X] (marginalizing over D for Robinson residuals)
  m_fmla <- stats::as.formula(paste("Y ~", paste(x_cols, collapse = " + ")))

  for (fold_k in seq_len(n_folds)) {
    tr <- fold_ids != fold_k
    te <- fold_ids == fold_k

    e_fit    <- stats::glm(e_fmla, data = dat[tr, ], family = stats::binomial())
    e_raw    <- stats::predict(e_fit, newdata = dat[te, ], type = "response")
    e_cf[te] <- pmin(pmax(e_raw, clip), 1 - clip)

    m_fit    <- stats::lm(m_fmla, data = dat[tr, ])
    m_cf[te] <- stats::predict(m_fit, newdata = dat[te, ])
  }

  rho_cf <- (dat$D - e_cf) * (dat$Y - m_cf) / pmax(e_cf * (1 - e_cf), 1e-8)

  fit <- cate_dnn(
    X              = as.matrix(dat[, x_cols]),
    x0             = dgp$x0,
    pseudo_outcome = rho_cf,
    k              = k,
    s              = s,
    se             = TRUE,
    jackknife      = "fast",
    return         = "ess"
  )
  .sim_metric(fit$tau, fit$se, dgp$tau_x0, alpha, neff = fit$ess)
}

# ── DNN oracle with Mahalanobis distance ─────────────────────────────────────
# Proposed method with Mahalanobis distance; relevant for dgp5_correlated where
# the off-diagonal covariance makes L2 neighbourhoods geometrically misaligned.
# The sample covariance is estimated once from the training data (via the
# distance_mahalanobis() setup() mechanism) and reused for all evaluations.

.method_dnn_oracle_mahal <- function(dgp, x_cols, k, s, alpha, ...) {
  fit <- cate_dnn(
    data      = dgp$data,
    x_cols    = x_cols,
    y_col     = "Y",
    t_col     = "D",
    e_col     = "e_oracle",
    m0_col    = "m0_oracle",
    m1_col    = "m1_oracle",
    x0        = dgp$x0,
    k         = k,
    s         = s,
    distance  = "mahalanobis",
    se        = TRUE,
    jackknife = "fast",
    return    = "ess"
  )
  .sim_metric(fit$tau, fit$se, dgp$tau_x0, alpha, neff = fit$ess)
}

# ── Local constant (Nadaraya-Watson) ──────────────────────────────────────────
# Product Gaussian kernel applied to oracle pseudo-outcomes.  Per-variable
# Silverman bandwidth.  SE via the vectorized delete-1 jackknife.
# Serves as the simplest nonparametric fixed-point comparator.

.method_local_constant <- function(dgp, x_cols, alpha, bw_factor = 1.06, ...) {
  dat <- dgp$data
  X   <- as.matrix(dat[, x_cols])
  rho <- dat$rho_oracle
  x0  <- dgp$x0
  n   <- nrow(X)
  kx  <- ncol(X)

  # Per-variable Silverman bandwidth: h_j = 1.06 * sd_j * n^{-1/(kx+4)}
  h_vec <- bw_factor * apply(X, 2, stats::sd) * n^(-1 / (kx + 4))
  h_vec <- pmax(h_vec, .Machine$double.eps)

  # Proportional product Gaussian kernel (constants cancel in NW ratio)
  Z      <- sweep(X, 2, x0) / rep(h_vec, each = n)
  K_vals <- exp(-0.5 * rowSums(Z^2))

  K_sum     <- sum(K_vals)
  K_rho_sum <- sum(K_vals * rho)
  tau_hat   <- K_rho_sum / K_sum

  # Vectorized LOO: tau_loo[i] = (K_rho_sum - K[i]*rho[i]) / (K_sum - K[i])
  denom_loo  <- K_sum - K_vals
  safe       <- denom_loo > .Machine$double.eps
  tau_loo    <- numeric(n)
  tau_loo[safe]  <- (K_rho_sum - K_vals[safe] * rho[safe]) / denom_loo[safe]
  tau_loo[!safe] <- tau_hat  # isolated dominant observation; leave-one-out = full

  se_hat <- sqrt((n - 1L) / n * sum((tau_loo - tau_hat)^2))

  neff <- K_sum^2 / sum(K_vals^2)

  .sim_metric(tau_hat, se_hat, dgp$tau_x0, alpha, neff = neff)
}

# ── OLS with treatment-covariate interactions ─────────────────────────────────
# Fits Y ~ D * (X1 + ... + Xkx), extracts tau(x0) via the plug-in formula,
# and computes the SE by the delta method.  Correctly specified only for DGP 1;
# tests global linear misspecification for DGPs 2 and 3.

.method_ols_interact <- function(dgp, x_cols, alpha, ...) {
  dat <- dgp$data
  x0  <- dgp$x0

  fmla <- stats::as.formula(paste(
    "Y ~ D * (", paste(x_cols, collapse = " + "), ")"
  ))
  fit <- stats::lm(fmla, data = dat)

  b       <- stats::coef(fit)
  int_nms <- paste0("D:", x_cols)

  # tau(x0) = beta_D + sum_j beta_{D:Xj} * x0[j]
  tau_hat <- as.numeric(b["D"]) + sum(b[int_nms] * x0, na.rm = TRUE)

  # Delta method SE: gradient g = (1, x0[1], ..., x0[kx]) for present coefs
  all_nms <- c("D", int_nms)
  present <- intersect(all_nms, rownames(stats::vcov(fit)))
  g_raw   <- stats::setNames(c(1, x0), all_nms)
  g       <- g_raw[present]
  se_hat  <- as.numeric(sqrt(t(g) %*% stats::vcov(fit)[present, present] %*% g))

  .sim_metric(tau_hat, se_hat, dgp$tau_x0, alpha)
}

# ── GRF causal forest ─────────────────────────────────────────────────────────
# Uses oracle nuisances (Y.hat = m_oracle, W.hat = e_oracle) to isolate
# the estimator's adaptation capability from first-stage quality.
# Primary comparison for the paper's self-normalisation / power claim.

.method_causal_forest <- function(dgp, x_cols, alpha, ...) {
  if (!requireNamespace("grf", quietly = TRUE)) {
    warning("grf not available; causal_forest returns NA.")
    return(.sim_metric(NA_real_, NA_real_, dgp$tau_x0, alpha))
  }
  dat    <- dgp$data
  X      <- as.matrix(dat[, x_cols])
  x0_mat <- matrix(dgp$x0, nrow = 1L)

  cf <- grf::causal_forest(
    X     = X,
    Y     = dat$Y,
    W     = dat$D,
    Y.hat = dat$m_oracle,
    W.hat = dat$e_oracle
  )
  pred    <- stats::predict(cf, newdata = x0_mat, estimate.variance = TRUE)
  tau_hat <- pred$predictions[1L]
  se_hat  <- sqrt(pred$variance.estimates[1L])

  .sim_metric(tau_hat, se_hat, dgp$tau_x0, alpha)
}

# ── GRF local linear forest ───────────────────────────────────────────────────
# Causal forest with local-linear corrections applied at prediction time via
# linear.correction.variables.  Oracle nuisances for first-stage isolation.

.method_llf <- function(dgp, x_cols, alpha, ll_lambda = 0.1, ...) {
  if (!requireNamespace("grf", quietly = TRUE)) {
    warning("grf not available; llf returns NA.")
    return(.sim_metric(NA_real_, NA_real_, dgp$tau_x0, alpha))
  }
  dat    <- dgp$data
  X      <- as.matrix(dat[, x_cols])
  kx     <- ncol(X)
  x0_mat <- matrix(dgp$x0, nrow = 1L)

  # ll.lambda / ll.weight.penalty are predict()-time options in grf >= 2.x;
  # passing them to causal_forest() breaks variance estimation.
  cf <- grf::causal_forest(
    X     = X,
    Y     = dat$Y,
    W     = dat$D,
    Y.hat = dat$m_oracle,
    W.hat = dat$e_oracle
  )
  pred <- stats::predict(
    cf,
    newdata                    = x0_mat,
    ll.lambda                  = ll_lambda,
    estimate.variance          = TRUE,
    linear.correction.variables = seq_len(kx)
  )
  tau_hat <- pred$predictions[1L]
  se_hat  <- sqrt(pred$variance.estimates[1L])

  .sim_metric(tau_hat, se_hat, dgp$tau_x0, alpha)
}

# ── Registry ──────────────────────────────────────────────────────────────────

.SIM_CATE_METHODS <- list(
  dnn_oracle       = .method_dnn_oracle,
  dnn_oracle_mahal = .method_dnn_oracle_mahal,
  dnn_crossfit     = .method_dnn_crossfit,
  local_constant   = .method_local_constant,
  ols_interact     = .method_ols_interact,
  causal_forest    = .method_causal_forest,
  llf              = .method_llf
)

.get_sim_method <- function(name) {
  fn <- .SIM_CATE_METHODS[[name]]
  if (is.null(fn))
    stop("Unknown simulation method: '", name, "'.  ",
         "Use sim_cate_methods() to see available methods.")
  fn
}

#' List available CATE simulation methods
#'
#' Returns the names of all estimator methods registered for use with
#' [sim_cate_one_rep()] and [sim_cate_run()].
#'
#' | Name                | Description |
#' |---------------------|-------------|
#' | `dnn_oracle`        | Proposed DNN Z-estimator, L2 distance, oracle nuisances |
#' | `dnn_oracle_mahal`  | Proposed DNN Z-estimator, Mahalanobis distance, oracle nuisances (`dgp5_correlated`) |
#' | `dnn_crossfit`      | Proposed DNN Z-estimator, L2 distance, K-fold cross-fitted nuisances |
#' | `local_constant`    | Nadaraya-Watson smoother on oracle pseudo-outcome |
#' | `ols_interact`      | OLS with D*X interactions, delta-method SE at x0 |
#' | `causal_forest`     | GRF causal forest, oracle nuisances (requires **grf**) |
#' | `llf`               | GRF causal forest with local-linear corrections (requires **grf**) |
#'
#' @return Character vector of registered method names.
#'
#' @seealso [sim_cate_run()], [sim_cate_one_rep()]
#'
#' @examples
#' sim_cate_methods()
#'
#' @export
sim_cate_methods <- function() names(.SIM_CATE_METHODS)
