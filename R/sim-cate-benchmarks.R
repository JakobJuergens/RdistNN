# sim-cate-benchmarks.R
# Benchmark estimator implementations for the CATE simulation study.
#
# All methods share a common internal interface:
#
#   .method_XXX(dgp, x_cols, k, s, tuning, tuning_mode, alpha, ...)
#     dgp    -- list returned by sim_cate_data()
#     x_cols -- character vector of covariate column names
#     k, s   -- scalar DNN tuning (only used by dnn_* methods in default mode)
#     tuning, tuning_mode -- optional DNN candidate grid controls
#     alpha  -- nominal coverage level
#     ...    -- method-specific extra arguments (absorbed silently if unused)
#
#   Returns a data.frame with one or more rows:
#     candidate_id, kappa, s, tau_hat, tau0, bias, sq_err, se, cover, ci_len,
#     neff, status
#
#   tau_hat / tau0   -- raw estimate and truth; self-contained for downstream
#   neff             -- neff_x = 1/sum(w_i^2(x)); paper eq. (02b L.81); NA for
#                       global methods
# Methods are registered in .SIM_CATE_METHODS and dispatched through
# .get_sim_method().  The public API is sim_cate_methods().
#
# GRF-based methods (causal_forest, llf) require the grf package listed
# under Suggests.  They return NA metrics with a warning when grf is absent.

# ── Shared metric constructors ────────────────────────────────────────────────

.sim_cate_metric_names <- function() {
  c("tau_hat", "tau0", "bias", "sq_err", "se", "cover", "ci_len", "neff")
}

.sim_cate_row_names <- function() {
  c("candidate_id", "kappa", "s", .sim_cate_metric_names(), "status")
}

.sim_recycle <- function(x, n) {
  if (length(x) == 0L) {
    return(rep(NA, n))
  }
  if (length(x) == n) {
    return(x)
  }
  rep(x, length.out = n)
}

.sim_metric <- function(
  tau_hat,
  se_hat,
  tau0,
  alpha,
  neff = NA_real_,
  candidate_id = NA_character_,
  kappa = NA_real_,
  s = NA_integer_,
  status = "ok"
) {
  n_out <- max(
    length(tau_hat), length(se_hat), length(tau0), length(neff),
    length(candidate_id), length(kappa), length(s), length(status)
  )
  tau_hat <- as.numeric(.sim_recycle(tau_hat, n_out))
  se_hat <- as.numeric(.sim_recycle(se_hat, n_out))
  tau0 <- as.numeric(.sim_recycle(tau0, n_out))
  neff <- as.numeric(.sim_recycle(neff, n_out))
  candidate_id <- as.character(.sim_recycle(candidate_id, n_out))
  kappa <- as.numeric(.sim_recycle(kappa, n_out))
  s <- as.integer(.sim_recycle(s, n_out))
  status <- as.character(.sim_recycle(status, n_out))

  z <- stats::qnorm(1 - alpha / 2)
  finite <- is.finite(tau_hat) & is.finite(se_hat) & se_hat > 0
  covered <- rep(NA_real_, n_out)
  if (any(finite)) {
    covered[finite] <- as.numeric(
      (tau_hat[finite] - z * se_hat[finite] <= tau0[finite]) &
        (tau0[finite] <= tau_hat[finite] + z * se_hat[finite])
    )
  }

  data.frame(
    candidate_id = candidate_id,
    kappa = kappa,
    s = s,
    tau_hat = tau_hat,
    tau0 = tau0,
    bias = as.numeric(tau_hat - tau0),
    sq_err = as.numeric((tau_hat - tau0)^2),
    se = se_hat,
    cover = covered,
    ci_len = as.numeric(2 * z * se_hat),
    neff = neff,
    status = status,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

.sim_cate_default_kappa <- function(tuning, n, kx, k) {
  if (is.null(tuning)) {
    return(NA_real_)
  }
  default <- tryCatch(
    tuning_eval(
      tuning = tuning,
      n = n,
      kx = kx,
      k = k,
      se = TRUE,
      jackknife = "fast",
      mode = "default"
    ),
    error = function(e) NULL
  )
  if (is.null(default) || is.null(default$meta$kappa)) {
    return(NA_real_)
  }
  out <- suppressWarnings(as.numeric(default$meta$kappa))
  if (length(out) != 1L || !is.finite(out)) NA_real_ else out
}

.sim_cate_dnn_grid_metrics <- function(fit, tau0, alpha) {
  grid <- fit$tuning_grid
  candidates <- grid$candidates
  est <- grid$estimates
  diag <- grid$diagnostics
  est_idx <- match(candidates$candidate_id, est$candidate_id)
  diag_idx <- match(candidates$candidate_id, diag$candidate_id)

  .sim_metric(
    tau_hat = est$estimate[est_idx],
    se_hat = est$se[est_idx],
    tau0 = tau0,
    alpha = alpha,
    neff = diag$ess[diag_idx],
    candidate_id = candidates$candidate_id,
    kappa = candidates$kappa,
    s = candidates$s,
    status = candidates$status
  )
}

# ── DNN oracle ────────────────────────────────────────────────────────────────
# Proposed method.  Uses the oracle AIPW score from the active coupled-score
# CATE route so that first-stage error does not contaminate the comparison and
# the oracle/crossfit comparison differs only in nuisance quality, not in score
# construction.

.method_dnn_oracle <- function(dgp, x_cols, k, s, tuning = NULL,
                               tuning_mode = "default", alpha, ...) {
  if (identical(tuning_mode, "grid")) {
    fit <- cate_dnn(
      X              = as.matrix(dgp$data[, x_cols]),
      x0             = dgp$x0,
      pseudo_outcome = dgp$data$aipw_oracle,
      k              = k,
      tuning         = tuning,
      tuning_mode    = "grid",
      se             = TRUE,
      jackknife      = "fast",
      return         = "ess"
    )
    return(.sim_cate_dnn_grid_metrics(fit, dgp$tau_x0, alpha))
  }

  fit <- cate_dnn(
    X              = as.matrix(dgp$data[, x_cols]),
    x0             = dgp$x0,
    pseudo_outcome = dgp$data$aipw_oracle,
    k              = k,
    s              = s,
    se             = TRUE,
    jackknife      = "fast",
    return         = "ess"
  )
  .sim_metric(
    fit$tau, fit$se, dgp$tau_x0, alpha,
    neff = fit$ess,
    kappa = .sim_cate_default_kappa(tuning, nrow(dgp$data), length(x_cols), k),
    s = s
  )
}

# ── DNN with estimated nuisances ──────────────────────────────────────────────
# Same nuisance model families are used for both same-sample and cross-fitted
# variants: logistic regression for the propensity score and arm-specific OLS
# for the potential outcome means.  The comparison between dnn_fullsample and
# dnn_crossfit then isolates the cost of sample splitting rather than the score
# construction.

.fit_sim_outcome_arm <- function(dat, x_cols, arm, train = rep(TRUE, nrow(dat))) {
  m_fmla <- stats::as.formula(paste("Y ~", paste(x_cols, collapse = " + ")))
  train_dat <- dat[train & dat$D == arm, , drop = FALSE]
  if (nrow(train_dat) == 0L) {
    .abort_bad_state(sprintf(
      "Cannot fit AIPW outcome nuisance: no D=%d observations in training sample.",
      arm
    ))
  }
  fit <- stats::lm(m_fmla, data = train_dat)
  as.numeric(stats::predict(fit, newdata = dat))
}

.fit_sim_aipw_nuisances <- function(dat, x_cols, train = rep(TRUE, nrow(dat)),
                                    clip = 1e-3) {
  e_fmla <- stats::as.formula(paste("D ~", paste(x_cols, collapse = " + ")))
  train_dat <- dat[train, , drop = FALSE]

  e_fit <- stats::glm(e_fmla, data = train_dat, family = stats::binomial())
  e_hat <- stats::predict(e_fit, newdata = dat, type = "response")
  e_hat <- pmin(pmax(e_hat, clip), 1 - clip)

  list(
    e_hat = e_hat,
    m0_hat = .fit_sim_outcome_arm(dat, x_cols = x_cols, arm = 0L, train = train),
    m1_hat = .fit_sim_outcome_arm(dat, x_cols = x_cols, arm = 1L, train = train)
  )
}

# ── DNN full-sample ───────────────────────────────────────────────────────────
# Same nuisance estimators as dnn_crossfit, but fit once on the full sample and
# reused on that same sample. This is the direct non-cross-fitted counterpart of
# the realistic cross-fit implementation.

.method_dnn_fullsample <- function(dgp, x_cols, k, s, tuning = NULL,
                                   tuning_mode = "default", alpha,
                                   clip = 1e-3, ...) {
  dat <- dgp$data
  nuis <- .fit_sim_aipw_nuisances(dat, x_cols = x_cols, clip = clip)
  aipw_fs <- cate_rho_aipw(
    y = dat$Y,
    d = dat$D,
    m1hat = nuis$m1_hat,
    m0hat = nuis$m0_hat,
    ehat = nuis$e_hat
  )

  if (identical(tuning_mode, "grid")) {
    fit <- cate_dnn(
      X              = as.matrix(dat[, x_cols]),
      x0             = dgp$x0,
      pseudo_outcome = aipw_fs,
      k              = k,
      tuning         = tuning,
      tuning_mode    = "grid",
      se             = TRUE,
      jackknife      = "fast",
      return         = "ess"
    )
    return(.sim_cate_dnn_grid_metrics(fit, dgp$tau_x0, alpha))
  }

  fit <- cate_dnn(
    X              = as.matrix(dat[, x_cols]),
    x0             = dgp$x0,
    pseudo_outcome = aipw_fs,
    k              = k,
    s              = s,
    se             = TRUE,
    jackknife      = "fast",
    return         = "ess"
  )
  .sim_metric(
    fit$tau, fit$se, dgp$tau_x0, alpha,
    neff = fit$ess,
    kappa = .sim_cate_default_kappa(tuning, nrow(dat), length(x_cols), k),
    s = s
  )
}

# ── DNN cross-fitted ──────────────────────────────────────────────────────────
# Proposed method with K-fold cross-fitted nuisances estimated via logistic
# regression (propensity) and arm-specific OLS (potential outcome means), both
# base-R.  Compares same-sample inference against the cross-fitting alternative.

.method_dnn_crossfit <- function(dgp, x_cols, k, s, tuning = NULL,
                                 tuning_mode = "default", alpha,
                                 n_folds = 5L, clip = 1e-3, ...) {
  dat <- dgp$data
  n   <- nrow(dat)

  # Balanced fold assignment; uses the RNG state continuing from the DGP draw,
  # so the result is deterministic for a given rep seed.
  fold_ids <- sample(rep_len(seq_len(n_folds), n))

  e_cf <- numeric(n)
  m0_cf <- numeric(n)
  m1_cf <- numeric(n)

  for (fold_k in seq_len(n_folds)) {
    te <- fold_ids == fold_k
    nuis <- .fit_sim_aipw_nuisances(
      dat = dat,
      x_cols = x_cols,
      train = !te,
      clip = clip
    )
    e_cf[te] <- nuis$e_hat[te]
    m0_cf[te] <- nuis$m0_hat[te]
    m1_cf[te] <- nuis$m1_hat[te]
  }

  aipw_cf <- cate_rho_aipw(
    y = dat$Y,
    d = dat$D,
    m1hat = m1_cf,
    m0hat = m0_cf,
    ehat = e_cf
  )

  if (identical(tuning_mode, "grid")) {
    fit <- cate_dnn(
      X              = as.matrix(dat[, x_cols]),
      x0             = dgp$x0,
      pseudo_outcome = aipw_cf,
      k              = k,
      tuning         = tuning,
      tuning_mode    = "grid",
      se             = TRUE,
      jackknife      = "fast",
      return         = "ess"
    )
    return(.sim_cate_dnn_grid_metrics(fit, dgp$tau_x0, alpha))
  }

  fit <- cate_dnn(
    X              = as.matrix(dat[, x_cols]),
    x0             = dgp$x0,
    pseudo_outcome = aipw_cf,
    k              = k,
    s              = s,
    se             = TRUE,
    jackknife      = "fast",
    return         = "ess"
  )
  .sim_metric(
    fit$tau, fit$se, dgp$tau_x0, alpha,
    neff = fit$ess,
    kappa = .sim_cate_default_kappa(tuning, nrow(dat), length(x_cols), k),
    s = s
  )
}

# ── DNN oracle with Mahalanobis distance ─────────────────────────────────────
# Proposed method with Mahalanobis distance and the oracle AIPW
# pseudo-outcome; relevant for dgp5_correlated where the off-diagonal covariance
# makes L2 neighbourhoods geometrically misaligned. The sample covariance is
# estimated once from the training data (via the distance_mahalanobis() setup()
# mechanism) and reused for all evaluations.

.method_dnn_oracle_mahal <- function(dgp, x_cols, k, s, tuning = NULL,
                                     tuning_mode = "default", alpha, ...) {
  if (identical(tuning_mode, "grid")) {
    fit <- cate_dnn(
      X              = as.matrix(dgp$data[, x_cols]),
      x0             = dgp$x0,
      pseudo_outcome = dgp$data$aipw_oracle,
      k              = k,
      tuning         = tuning,
      tuning_mode    = "grid",
      distance       = "mahalanobis",
      se             = TRUE,
      jackknife      = "fast",
      return         = "ess"
    )
    return(.sim_cate_dnn_grid_metrics(fit, dgp$tau_x0, alpha))
  }

  fit <- cate_dnn(
    X              = as.matrix(dgp$data[, x_cols]),
    x0             = dgp$x0,
    pseudo_outcome = dgp$data$aipw_oracle,
    k              = k,
    s              = s,
    distance       = "mahalanobis",
    se             = TRUE,
    jackknife      = "fast",
    return         = "ess"
  )
  .sim_metric(
    fit$tau, fit$se, dgp$tau_x0, alpha,
    neff = fit$ess,
    kappa = .sim_cate_default_kappa(tuning, nrow(dgp$data), length(x_cols), k),
    s = s
  )
}

# ── Local constant (Nadaraya-Watson) ──────────────────────────────────────────
# Product Gaussian kernel applied to oracle AIPW pseudo-outcomes.  Per-variable
# Silverman bandwidth.  SE via the vectorized delete-1 jackknife.
# Serves as the simplest nonparametric fixed-point comparator.

.method_local_constant <- function(dgp, x_cols, alpha, bw_factor = 1.06, ...) {
  dat <- dgp$data
  X   <- as.matrix(dat[, x_cols])
  aipw <- dat$aipw_oracle
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
  K_aipw_sum <- sum(K_vals * aipw)
  tau_hat    <- K_aipw_sum / K_sum

  # Vectorized LOO: tau_loo[i] = (K_aipw_sum - K[i]*aipw[i]) / (K_sum - K[i])
  denom_loo  <- K_sum - K_vals
  safe       <- denom_loo > .Machine$double.eps
  tau_loo    <- numeric(n)
  tau_loo[safe]  <- (K_aipw_sum - K_vals[safe] * aipw[safe]) / denom_loo[safe]
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
# Use one GRF thread per Monte Carlo worker: the outer simulation harness
# parallelizes over replications, so letting each forest fit use all cores would
# oversubscribe the machine badly on local multisession runs.

.method_causal_forest <- function(dgp, x_cols, alpha, ...) {
  if (!requireNamespace("grf", quietly = TRUE)) {
    warning("grf not available; causal_forest returns NA.", call. = FALSE)
    return(.sim_metric(NA_real_, NA_real_, dgp$tau_x0, alpha, status = "failed"))
  }
  dat    <- dgp$data
  X      <- as.matrix(dat[, x_cols])
  x0_mat <- matrix(dgp$x0, nrow = 1L)

  cf <- grf::causal_forest(
    X     = X,
    Y     = dat$Y,
    W     = dat$D,
    Y.hat = dat$m_oracle,
    W.hat = dat$e_oracle,
    num.threads = 1L
  )
  pred    <- stats::predict(
    cf,
    newdata = x0_mat,
    estimate.variance = TRUE,
    num.threads = 1L
  )
  tau_hat <- pred$predictions[1L]
  se_hat  <- sqrt(pred$variance.estimates[1L])

  .sim_metric(tau_hat, se_hat, dgp$tau_x0, alpha)
}

# ── GRF local linear forest ───────────────────────────────────────────────────
# Causal forest with local-linear corrections applied at prediction time via
# linear.correction.variables.  Oracle nuisances for first-stage isolation.
# Keep GRF itself single-threaded here for the same reason as above: rep-level
# parallelism lives outside this method.

.method_llf <- function(dgp, x_cols, alpha, ll_lambda = 0.1, ...) {
  if (!requireNamespace("grf", quietly = TRUE)) {
    warning("grf not available; llf returns NA.", call. = FALSE)
    return(.sim_metric(NA_real_, NA_real_, dgp$tau_x0, alpha, status = "failed"))
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
    W.hat = dat$e_oracle,
    num.threads = 1L
  )
  pred <- stats::predict(
    cf,
    newdata                    = x0_mat,
    ll.lambda                  = ll_lambda,
    estimate.variance          = TRUE,
    linear.correction.variables = seq_len(kx),
    num.threads                = 1L
  )
  tau_hat <- pred$predictions[1L]
  se_hat  <- sqrt(pred$variance.estimates[1L])

  .sim_metric(tau_hat, se_hat, dgp$tau_x0, alpha)
}

# ── Registry ──────────────────────────────────────────────────────────────────

.SIM_CATE_METHODS <- list(
  dnn_oracle       = .method_dnn_oracle,
  dnn_oracle_mahal = .method_dnn_oracle_mahal,
  dnn_fullsample   = .method_dnn_fullsample,
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
         "Use sim_cate_methods() to see available methods.", call. = FALSE)
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
#' | `dnn_fullsample`    | Proposed DNN Z-estimator, L2 distance, full-sample estimated AIPW nuisances |
#' | `dnn_crossfit`      | Proposed DNN Z-estimator, L2 distance, K-fold cross-fitted AIPW nuisances |
#' | `local_constant`    | Nadaraya-Watson smoother on the oracle AIPW score |
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
