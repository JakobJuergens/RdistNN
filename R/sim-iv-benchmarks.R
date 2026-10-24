# sim-iv-benchmarks.R
# Oracle and classical benchmark helpers for IV-GMM simulations.
#
# Methods are registered in .SIM_IV_METHODS and dispatched through
# sim_iv_one_rep()/sim_iv_run().

.sim_metric_iv <- function(theta_hat, se_hat, theta0, alpha, neff = NA_real_,
                           J = NA_real_, df = NA_integer_) {
  z <- .dnn_normal_critical_alpha(alpha)
  finite_ci <- is.finite(theta_hat) && is.finite(se_hat) && se_hat > 0
  covered <- if (finite_ci) {
    (theta_hat - z * se_hat <= theta0) & (theta0 <= theta_hat + z * se_hat)
  } else {
    NA_real_
  }

  valid_j <- is.finite(J) && !is.na(df) && df > 0L
  pval_J <- if (valid_j) stats::pchisq(J, df = df, lower.tail = FALSE) else NA_real_
  reject_J <- if (is.finite(pval_J)) as.numeric(pval_J < alpha) else NA_real_

  c(
    theta_hat = as.numeric(theta_hat),
    theta0 = as.numeric(theta0),
    bias = as.numeric(theta_hat - theta0),
    sq_err = as.numeric((theta_hat - theta0)^2),
    se = as.numeric(se_hat),
    cover = as.numeric(covered),
    ci_len = as.numeric(2 * z * se_hat),
    neff = as.numeric(neff),
    J = as.numeric(J),
    df = as.numeric(if (!is.na(df) && df > 0L) df else NA_integer_),
    pval_J = as.numeric(pval_J),
    reject_J = as.numeric(reject_J)
  )
}

.sim_iv_colnames <- function(dgp, prefix) {
  grep(paste0("^", prefix, "[0-9]+"), names(dgp$data), value = TRUE)
}

.sim_iv_na_metric <- function(theta0, alpha, reason) {
  out <- .sim_metric_iv(NA_real_, NA_real_, theta0, alpha)
  attr(out, "method_failures") <- data.frame(method = NA_character_,
    reason = .dnn_brief(reason, 240L), stringsAsFactors = FALSE)
  out
}

.sim_iv_safe_solve <- function(A, b = NULL) {
  # Match solve.default's existing estimated 1-norm reciprocal-condition
  # cutoff. Keep host handling of nonfinite matrices and unexpected failures.
  if (all(is.finite(A)) && rcond(A, norm = "1") < .Machine$double.eps) return(NULL)
  if (is.null(b)) solve(A) else solve(A, b)
}

.sim_iv_fit_lm <- function(formula, train, test) {
  fit <- stats::lm(formula, data = train)
  as.numeric(stats::predict(fit, newdata = test))
}

.sim_iv_fit_glm_prob <- function(formula, train, test, clip) {
  fit <- stats::glm(formula, data = train, family = stats::binomial())
  pred <- as.numeric(stats::predict(fit, newdata = test, type = "response"))
  pmin(pmax(pred, clip), 1 - clip)
}

.sim_iv_instrument_types <- function(dgp) {
  types <- dgp$truth$instrument_type
  if (is.null(types)) {
    types <- rep("continuous", dgp$q)
  }
  rep_len(types, dgp$q)
}

.fit_sim_iv_nuisances_fullsample <- function(dat, x_cols, q, instrument_type,
                                             clip = 1e-3) {
  z_cols <- paste0("Z", seq_len(q))
  pi_hat <- matrix(NA_real_, nrow(dat), q)

  for (j in seq_len(q)) {
    fmla <- stats::reformulate(x_cols, response = z_cols[j])
    pi_hat[, j] <- if (identical(instrument_type[j], "binary")) {
      .sim_iv_fit_glm_prob(fmla, train = dat, test = dat, clip = clip)
    } else {
      .sim_iv_fit_lm(fmla, train = dat, test = dat)
    }
  }

  muY <- .sim_iv_fit_lm(stats::reformulate(x_cols, response = "Y"), dat, dat)
  muD <- .sim_iv_fit_glm_prob(stats::reformulate(x_cols, response = "D"), dat, dat, clip)

  list(pi = pi_hat, muY = muY, muD = muD)
}

.fit_sim_iv_nuisances_crossfit <- function(dat, x_cols, q, instrument_type,
                                           n_folds = 5L, clip = 1e-3) {
  n <- nrow(dat)
  n_folds <- .validate_scalar_int(n_folds, arg = "n_folds", min = 2L)
  if (n_folds > n) {
    .abort_bad_arg(
      "n_folds",
      "must be <= nrow(data)",
      expected = sprintf("<= %d", n),
      actual = sprintf("n_folds=%d", n_folds)
    )
  }

  z_cols <- paste0("Z", seq_len(q))
  fold_ids <- sample(rep_len(seq_len(n_folds), n))
  pi_hat <- matrix(NA_real_, n, q)
  muY <- muD <- rep(NA_real_, n)

  for (fold_k in seq_len(n_folds)) {
    tr <- fold_ids != fold_k
    te <- fold_ids == fold_k
    train <- dat[tr, , drop = FALSE]
    test <- dat[te, , drop = FALSE]

    for (j in seq_len(q)) {
      fmla <- stats::reformulate(x_cols, response = z_cols[j])
      pi_hat[te, j] <- if (identical(instrument_type[j], "binary")) {
        .sim_iv_fit_glm_prob(fmla, train = train, test = test, clip = clip)
      } else {
        .sim_iv_fit_lm(fmla, train = train, test = test)
      }
    }
    muY[te] <- .sim_iv_fit_lm(stats::reformulate(x_cols, response = "Y"), train, test)
    muD[te] <- .sim_iv_fit_glm_prob(stats::reformulate(x_cols, response = "D"), train, test, clip)
  }

  list(pi = pi_hat, muY = muY, muD = muD)
}

.sim_iv_add_nuisance_cols <- function(dat, nuis, prefix = "hat") {
  out <- dat
  q <- ncol(nuis$pi)
  pi_cols <- paste0("piZ", seq_len(q), "_", prefix)
  for (j in seq_len(q)) {
    out[[pi_cols[j]]] <- nuis$pi[, j]
  }
  out[[paste0("muY_", prefix)]] <- nuis$muY
  out[[paste0("muD_", prefix)]] <- nuis$muD
  list(
    data = out,
    pi_cols = pi_cols,
    muY_col = paste0("muY_", prefix),
    muD_col = paste0("muD_", prefix)
  )
}

.sim_iv_gmm_metric <- function(dgp, x_cols, alpha, pi_cols, muY_col, muD_col,
                               k, s, distance = "l2", method = "exact") {
  z_cols <- paste0("Z", seq_len(dgp$q))
  model <- gmm_model_dr_iv(
    y_col = "Y",
    d_col = "D",
    z_cols = z_cols,
    pi_cols = pi_cols,
    muY_col = muY_col,
    muD_col = muD_col
  )
  fit <- gmm_dnn(
    data = dgp$data,
    x_cols = x_cols,
    x0 = dgp$x0,
    model = model,
    distance = distance,
    method = method,
    k = k,
    s = s,
    se = TRUE,
    return = "ess"
  )
  .sim_metric_iv(
    theta_hat = fit$theta[1L, 1L],
    se_hat = fit$se[1L, 1L],
    theta0 = dgp$theta_x0,
    alpha = alpha,
    neff = fit$neff[1L],
    J = fit$J[1L],
    df = fit$df[1L]
  )
}

.method_iv_oracle <- function(dgp, x_cols, k, s, alpha,
                              distance = "l2", method = "exact", ...) {
  dat <- dgp$data
  X <- as.matrix(dat[, x_cols, drop = FALSE])
  z_cols <- paste0("Z", seq_len(dgp$q))
  pi_cols <- paste0("piZ", seq_len(dgp$q), "_oracle")

  wt <- weights_dnn(
    X = X,
    x0 = dgp$x0,
    distance = distance,
    method = method,
    k = k,
    s = s,
    return = "W"
  )
  w <- as.numeric(wt$W[, 1L])
  neff <- 1 / sum(w^2)

  Z_resid <- as.matrix(dat[, z_cols, drop = FALSE]) -
    as.matrix(dat[, pi_cols, drop = FALSE])
  y_resid <- dat$Y - dat$muY_oracle
  d_resid <- dat$D - dat$muD_oracle

  A_obs <- Z_resid * y_resid
  B_obs <- Z_resid * d_resid
  a_bar <- as.numeric(crossprod(A_obs, w))
  b_bar <- as.numeric(crossprod(B_obs, w))

  denom_identity <- sum(b_bar^2)
  if (!is.finite(denom_identity) || denom_identity <= .Machine$double.eps) {
    return(.sim_iv_na_metric(dgp$theta_x0, alpha,
      "The oracle IV denominator is nonfinite or no larger than the existing numerical tolerance."))
  }
  theta_1 <- sum(b_bar * a_bar) / denom_identity
  g1 <- A_obs - B_obs * theta_1
  Omega <- crossprod(g1, g1 * w)
  Omega <- (Omega + t(Omega)) / 2
  Omega_inv <- .sim_iv_safe_solve(Omega)
  if (is.null(Omega_inv) || any(!is.finite(Omega_inv))) {
    return(.sim_iv_na_metric(dgp$theta_x0, alpha,
      "The oracle IV score covariance is singular or its inverse is nonfinite."))
  }

  info <- as.numeric(t(b_bar) %*% Omega_inv %*% b_bar)
  if (!is.finite(info) || info <= .Machine$double.eps) {
    return(.sim_iv_na_metric(dgp$theta_x0, alpha,
      "The oracle IV information is nonfinite or no larger than the existing numerical tolerance."))
  }
  theta_hat <- as.numeric(t(b_bar) %*% Omega_inv %*% a_bar / info)
  se_hat <- sqrt(1 / (neff * info))

  J <- NA_real_
  df <- NA_integer_
  if (dgp$q > 1L) {
    gbar <- a_bar - b_bar * theta_hat
    J <- as.numeric(neff * t(gbar) %*% Omega_inv %*% gbar)
    df <- as.integer(dgp$q - 1L)
  }

  .sim_metric_iv(theta_hat, se_hat, dgp$theta_x0, alpha, neff = neff, J = J, df = df)
}

.method_dnn_iv_oracle <- function(dgp, x_cols, k, s, alpha,
                                  distance = "l2", method = "exact", ...) {
  .sim_iv_gmm_metric(
    dgp = dgp,
    x_cols = x_cols,
    alpha = alpha,
    pi_cols = paste0("piZ", seq_len(dgp$q), "_oracle"),
    muY_col = "muY_oracle",
    muD_col = "muD_oracle",
    k = k,
    s = s,
    distance = distance,
    method = method
  )
}

.method_dnn_iv_oracle_mahal <- function(dgp, x_cols, k, s, alpha, ...) {
  .method_dnn_iv_oracle(
    dgp = dgp,
    x_cols = x_cols,
    k = k,
    s = s,
    alpha = alpha,
    distance = "mahalanobis",
    ...
  )
}

.method_dnn_iv_fullsample <- function(dgp, x_cols, k, s, alpha,
                                      clip = 1e-3, ...) {
  instrument_type <- .sim_iv_instrument_types(dgp)
  nuis <- .fit_sim_iv_nuisances_fullsample(
      dat = dgp$data,
      x_cols = x_cols,
      q = dgp$q,
      instrument_type = instrument_type,
      clip = clip
  )
  if (is.null(nuis)) return(.sim_iv_na_metric(dgp$theta_x0, alpha,
      "Full-sample nuisance estimation returned no fitted predictions."))

  augmented <- .sim_iv_add_nuisance_cols(dgp$data, nuis, prefix = "fs")
  dgp_aug <- dgp
  dgp_aug$data <- augmented$data
  .sim_iv_gmm_metric(
    dgp = dgp_aug,
    x_cols = x_cols,
    alpha = alpha,
    pi_cols = augmented$pi_cols,
    muY_col = augmented$muY_col,
    muD_col = augmented$muD_col,
    k = k,
    s = s
  )
}

.method_dnn_iv_crossfit <- function(dgp, x_cols, k, s, alpha,
                                    n_folds = 5L, clip = 1e-3, ...) {
  instrument_type <- .sim_iv_instrument_types(dgp)
  nuis <- .fit_sim_iv_nuisances_crossfit(
      dat = dgp$data,
      x_cols = x_cols,
      q = dgp$q,
      instrument_type = instrument_type,
      n_folds = n_folds,
      clip = clip
  )
  if (is.null(nuis)) return(.sim_iv_na_metric(dgp$theta_x0, alpha,
      "Cross-fitted nuisance estimation returned no fitted predictions."))

  augmented <- .sim_iv_add_nuisance_cols(dgp$data, nuis, prefix = "cf")
  dgp_aug <- dgp
  dgp_aug$data <- augmented$data
  .sim_iv_gmm_metric(
    dgp = dgp_aug,
    x_cols = x_cols,
    alpha = alpha,
    pi_cols = augmented$pi_cols,
    muY_col = augmented$muY_col,
    muD_col = augmented$muD_col,
    k = k,
    s = s
  )
}

.method_2sls_global <- function(dgp, x_cols, alpha, ...) {
  dat <- dgp$data
  z_cols <- paste0("Z", seq_len(dgp$q))
  fs_formula <- stats::reformulate(c(z_cols, x_cols), response = "D")
  fs <- stats::lm(fs_formula, data = dat)
  if (is.null(fs)) return(.sim_iv_na_metric(dgp$theta_x0, alpha,
      "The first-stage regression returned no fitted model."))

  df2 <- dat
  df2$Dhat <- stats::predict(fs, newdata = dat)
  ss_formula <- stats::reformulate(c("Dhat", x_cols), response = "Y")
  ss <- stats::lm(ss_formula, data = df2)
  if (is.null(ss)) return(.sim_iv_na_metric(dgp$theta_x0, alpha,
      "The second-stage regression returned no fitted model."))

  co <- stats::coef(ss)
  vc <- stats::vcov(ss)
  if (!"Dhat" %in% names(co) || !is.finite(co[["Dhat"]]) ||
      !"Dhat" %in% rownames(vc)) {
    return(.sim_iv_na_metric(dgp$theta_x0, alpha,
      "The second-stage treatment coefficient or covariance is unavailable."))
  }

  .sim_metric_iv(
    theta_hat = co[["Dhat"]],
    se_hat = sqrt(vc["Dhat", "Dhat"]),
    theta0 = dgp$theta_x0,
    alpha = alpha
  )
}

.method_2sls_interact <- function(dgp, x_cols, alpha, ...) {
  dat <- dgp$data
  z_cols <- paste0("Z", seq_len(dgp$q))

  stage_df <- dat[, c("Y", "D", x_cols, z_cols), drop = FALSE]
  endog_cols <- c("D", paste0("D_", x_cols))
  stage_df[endog_cols[1L]] <- dat$D
  for (j in seq_along(x_cols)) {
    stage_df[[endog_cols[j + 1L]]] <- dat$D * dat[[x_cols[j]]]
  }

  inst_cols <- c(x_cols, z_cols)
  for (z in z_cols) {
    for (x in x_cols) {
      nm <- paste0(z, "_", x)
      stage_df[[nm]] <- dat[[z]] * dat[[x]]
      inst_cols <- c(inst_cols, nm)
    }
  }

  fit_cols <- paste0("hat_", endog_cols)
  for (j in seq_along(endog_cols)) {
    fs_formula <- stats::reformulate(inst_cols, response = endog_cols[j])
    fs <- stats::lm(fs_formula, data = stage_df)
    if (is.null(fs)) return(.sim_iv_na_metric(dgp$theta_x0, alpha,
      "An interacted first-stage regression returned no fitted model."))
    stage_df[[fit_cols[j]]] <- stats::predict(fs, newdata = stage_df)
  }

  ss_formula <- stats::reformulate(c(x_cols, fit_cols), response = "Y")
  ss <- stats::lm(ss_formula, data = stage_df)
  if (is.null(ss)) return(.sim_iv_na_metric(dgp$theta_x0, alpha,
      "The interacted second-stage regression returned no fitted model."))

  co <- stats::coef(ss)
  vc <- stats::vcov(ss)
  needed <- fit_cols
  if (!all(needed %in% names(co)) || any(!is.finite(co[needed]))) {
    return(.sim_iv_na_metric(dgp$theta_x0, alpha,
      "A required interacted second-stage coefficient is unavailable."))
  }

  grad <- stats::setNames(c(1, dgp$x0), needed)
  theta_hat <- sum(grad * co[needed])
  vc_sub <- vc[needed, needed, drop = FALSE]
  se_hat <- as.numeric(sqrt(t(grad) %*% vc_sub %*% grad))

  .sim_metric_iv(theta_hat, se_hat, dgp$theta_x0, alpha)
}

.method_local_2sls <- function(dgp, x_cols, alpha, bw_factor = 1.06, ...) {
  dat <- dgp$data
  X <- as.matrix(dat[, x_cols, drop = FALSE])
  x0 <- dgp$x0
  n <- nrow(X)
  kx <- ncol(X)
  z_cols <- paste0("Z", seq_len(dgp$q))

  h_vec <- bw_factor * apply(X, 2, stats::sd) * n^(-1 / (kx + 4))
  h_vec <- pmax(h_vec, .Machine$double.eps)
  Zloc <- sweep(X, 2, x0) / rep(h_vec, each = n)
  K_vals <- exp(-0.5 * rowSums(Zloc^2))
  if (!is.finite(sum(K_vals)) || sum(K_vals) <= .Machine$double.eps) {
    return(.sim_iv_na_metric(dgp$theta_x0, alpha,
      "The local IV kernel mass is nonfinite or no larger than the existing numerical tolerance."))
  }

  fs_formula <- stats::reformulate(z_cols, response = "D")
  fs <- stats::lm(fs_formula, data = dat, weights = K_vals)
  if (is.null(fs)) return(.sim_iv_na_metric(dgp$theta_x0, alpha,
      "The local first-stage regression returned no fitted model."))

  df2 <- dat
  df2$Dhat <- stats::predict(fs, newdata = dat)
  ss <- stats::lm(Y ~ Dhat, data = df2, weights = K_vals)
  if (is.null(ss)) return(.sim_iv_na_metric(dgp$theta_x0, alpha,
      "The local second-stage regression returned no fitted model."))

  co <- stats::coef(ss)
  vc <- stats::vcov(ss)
  if (!"Dhat" %in% names(co) || !is.finite(co[["Dhat"]]) ||
      !"Dhat" %in% rownames(vc)) {
    return(.sim_iv_na_metric(dgp$theta_x0, alpha,
      "The local second-stage treatment coefficient or covariance is unavailable."))
  }

  neff <- sum(K_vals)^2 / sum(K_vals^2)
  .sim_metric_iv(co[["Dhat"]], sqrt(vc["Dhat", "Dhat"]), dgp$theta_x0, alpha, neff = neff)
}

.SIM_IV_METHODS <- list(
  dnn_iv_oracle = .method_dnn_iv_oracle,
  dnn_iv_oracle_mahal = .method_dnn_iv_oracle_mahal,
  dnn_iv_fullsample = .method_dnn_iv_fullsample,
  dnn_iv_crossfit = .method_dnn_iv_crossfit,
  iv_oracle = .method_iv_oracle,
  `2sls_global` = .method_2sls_global,
  `2sls_interact` = .method_2sls_interact,
  local_2sls = .method_local_2sls
)

.get_sim_iv_method <- function(name) {
  fn <- .SIM_IV_METHODS[[name]]
  if (is.null(fn)) {
    .dnn_abort("argument_invalid", paste0(
      "Unknown IV simulation method: '", name, "'. Valid options: ",
      paste(names(.SIM_IV_METHODS), collapse = ", "), "."), argument = "name")
  }
  fn
}

#' List available IV simulation methods
#'
#' Returns the names of all estimator methods registered for use with
#' [sim_iv_one_rep()] and [sim_iv_run()].
#'
#' @return Character vector of registered IV simulation method names.
#' @examples
#' sim_iv_methods()
#' @export
sim_iv_methods <- function() names(.SIM_IV_METHODS)
