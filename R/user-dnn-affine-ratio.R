#' DNN inference for a ratio of two supplied signal averages
#'
#' Localizes two aligned, externally supplied signals with common DNN
#' weights. This estimator uses paired frozen-score covariance;
#' it never divides individual signals or deleted component averages.
#'
#' @param data Optional data frame or tibble. Do not combine with X or direct
#'   sample vectors.
#' @param x_cols Reporting-covariate columns in data, in target-coordinate order.
#' @param numerator_col,denominator_col Signal columns for the data-frame route.
#' @param X Numeric reporting matrix, observations in rows.
#' @param X0 Numeric target matrix, targets in rows.
#' @param x0 Numeric vector for a single target. Supply exactly one of X0,x0.
#' @param numerator,denominator Finite numeric signal vectors for the matrix route.
#' @param distance Existing distance specification, registered name or function.
#'   Geometry is frozen across deletions; input position resolves distance ties.
#' @param method One of "exact" or "asymp", as in [cate_dnn()] and [np_dnn()].
#' @param k Number of neighbors averaged within each size-s subsample, 1 <= k <= s.
#' @param s Deterministic integer order. Covariance requires s < n.
#' @param tuning Optional existing tuning specification, mutually exclusive with s.
#'   If neither is provided, the package paper rule resolves the order.
#' @param tuning_mode "default" retains a single fit; "grid" also evaluates
#'   feasible tuning candidates using the same frozen signals and ranks.
#' @param normalize Logical normalization control passed to the asymptotic
#'   weight rule, as in [weights_dnn()]. Exact weights are normalized by construction.
#' @param se Whether to compute covariance.
#' @param jackknife "fast" uses the existing prefix/suffix engine; "loo" uses
#'   explicit frozen-rank deletions; "none" disables covariance.
#' @param denominator_tol Nonnegative numerical cancellation tolerance. Changes
#'   a diagnostic flag only, never the ratio or covariance.
#' @param row_ids Optional aligned unique character/numeric identities.
#' @param ctx Optional context passed to the existing distance machinery.
#' @param return Optional extras: W, ranks, ord, dx, ess, signals, deletions,
#'   pseudo_values, neighbors or tuning. Unknown values warn and are ignored.
#'
#' @details Both full and deleted averages use the same order s. The target
#' and nuisance predictions remain frozen. If V is the paired estimator-scale
#' component covariance, the ratio variance is
#' \deqn{(V_{YY}+\hat\theta^2 V_{DD}-2\hat\theta V_{YD})/\hat\Delta_D^2.}
#' The computation uses centered deletion scores to retain cross-covariance.
#' theta_vcov is already on the estimator scale; theta_omega multiplies it by
#' neff. Do not divide theta_vcov by neff again.
#'
#' Signed signals and finite negative denominators are retained. At exact zero
#' denominator the components remain available and the ratio is NA. Computational
#' statuses do not verify causal assumptions, learner stability, point-target
#' undersmoothing or weak-instrument-uniform validity.
#' If nonzero component variation underflows its covariance, uncertainty is
#' unavailable with numerical_failure status; it is not treated as deterministic.
#'
#' Prediction retains the resolved distance specification and fitted context;
#' re-registering the distance name does not change an existing fit.
#'
#' @return An affine_ratio_result with theta, se, m-by-2 components,
#'   lists component_vcov, theta_vcov and theta_omega (NULL without covariance),
#'   x0, neff, per-target status, row_ids, meta and retained frozen model state.
#'   Optional deletions/pseudo_values have dimensions n-by-2-by-m. Covariance
#'   matrices are marginal across targets, not cross-target covariances.
#' @examples
#' fit <- affine_ratio_dnn(
#'   X = matrix(1:4, ncol = 1), x0 = 0,
#'   numerator = c(2, 1, 4, 3), denominator = c(1, -1, 2, 1), s = 2
#' )
#' fit$theta
#' fit$theta_vcov
#' d <- data.frame(x = 1:4, y = c(2, 1, 4, 3), d = c(1, -1, 2, 1))
#' affine_ratio_dnn(
#'   data = d, x_cols = "x", numerator_col = "y",
#'   denominator_col = "d", x0 = 0, s = 2, se = FALSE
#' )$theta
#' @export
affine_ratio_dnn <- function(
  data = NULL, x_cols = NULL,
  numerator_col = NULL, denominator_col = NULL,
  X = NULL, X0 = NULL, x0 = NULL,
  numerator = NULL, denominator = NULL,
  distance = "l2", method = c("exact", "asymp"), k = 1L,
  s = NULL, tuning = NULL, tuning_mode = c("default", "grid"),
  normalize = TRUE,
  se = TRUE, jackknife = c("fast", "loo", "none"),
  denominator_tol = sqrt(.Machine$double.eps),
  row_ids = NULL, ctx = NULL, return = c("W")
) {
  method <- .affine_weight_method(method)
  tuning_mode <- match.arg(tuning_mode)
  if (tuning_mode == "grid" && !is.null(s)) {
    .abort_bad_arg("tuning_mode", "cannot be 'grid' with explicit s; omit s and use tuning")
  }
  k <- .affine_integer(k, "k")
  normalize <- .validate_scalar_logical(normalize, "normalize")
  se <- .validate_scalar_logical(se, "se")
  jackknife <- match.arg(jackknife)
  denominator_tol <- .affine_nonnegative_scalar(denominator_tol, "denominator_tol")
  valid_returns <- c(
    "W", "ranks", "ord", "dx", "ess", "signals",
    "deletions", "pseudo_values", "neighbors", "tuning"
  )
  .warn_unknown_return(return, valid = valid_returns)
  return <- intersect(return, valid_returns)
  resolved_x_cols <- NULL
  if (!is.null(data)) {
    if (!is.null(X) || !is.null(numerator) || !is.null(denominator)) {
      .abort_bad_arg("data/X", "cannot combine data with X or direct signal vectors", quoted = FALSE)
    }
    if (!is.data.frame(data)) .abort_bad_arg("data", "must be a data.frame")
    n <- .affine_sample_size(data, "data")
    resolved_x_cols <- .dnn_resolve_data_cols(x_cols, data, "x_cols")
    x_mat <- .dnn_col_mat(resolved_x_cols, data, "x_cols", out_arg = "X")
    num <- .affine_column(numerator_col, data, n, "numerator_col")
    den <- .affine_column(denominator_col, data, n, "denominator_col")
  } else {
    if (!is.null(x_cols) || !is.null(numerator_col) || !is.null(denominator_col)) {
      .abort_bad_arg("data", "must be supplied when column selectors are used")
    }
    if (is.null(X)) .abort_bad_arg("X", "must be supplied for the matrix interface")
    x_mat <- .validate_matrix_nonempty(X, "X")
    n <- .affine_sample_size(x_mat, "X")
    num <- .affine_vector(numerator, n, "numerator")
    den <- .affine_vector(denominator, n, "denominator")
  }
  .abort_if_na_or_nonfinite(x_mat, "X", quoted = TRUE)
  ids <- .affine_row_ids(row_ids, n)
  x0_mat <- .dnn_resolve_targets(x_mat, X0, x0)
  active_covariance <- se && jackknife != "none"
  if (!is.null(s)) {
    s <- .affine_integer(s, "s")
    if (s > n) .abort_bad_arg("s", sprintf("must satisfy s <= n=%d", n))
  }
  resolved <- .dnn_resolve_s(s, tuning, n, ncol(x_mat), k, se, jackknife, ctx)
  s <- resolved$s
  if (s > n) .abort_bad_arg("s", sprintf("must satisfy s <= n=%d", n))
  if (active_covariance && (n < 2L || s >= n)) {
    .abort_bad_arg("s", sprintf("must satisfy s < n=%d for delete-1 covariance", n))
  }
  geometry <- .dnn_fitted_geometry(distance, x_mat, ctx)
  # Resolve the geometry once, including for the explicit deletion reference.
  weights <- .dnn_weights_with_extras(
    x_mat, x0_mat, geometry$distance, method, k, s, normalize, geometry$ctx,
    return = c("ranks", "ord", "dx")
  )
  W <- weights$W
  wt <- weights$wt
  m <- ncol(W)
  signal <- cbind(numerator = num, denominator = den)
  model_weights <- list(k = k, s = s, method = method, normalize = normalize)
  out <- .affine_ratio_evaluate(signal, W, wt$ranks, .affine_rank_weight(model_weights),
                                active_covariance, jackknife, denominator_tol)
  loo <- out$deletions
  out$deletions <- NULL
  components <- out$components
  out$x0 <- .dnn_x0_with_names(x0_mat, x_mat)
  out$row_ids <- ids
  out$meta <- list(
    n = n, m = m, p = ncol(x_mat), k = k, s = s, method = method, normalize = normalize,
    distance = .dnn_distance_name(distance), jackknife = if (active_covariance) jackknife else "none",
    covariance = if (active_covariance) "paired_frozen_score" else "none",
    component_roles = c(numerator = "numerator", denominator = "denominator"),
    tie_rule = "distance_then_input_position", tuning = resolved$tuning, tuning_mode = tuning_mode
  )
  out <- .dnn_attach_common_extras(out, return, wt, W, k, resolved$tuning)
  if ("signals" %in% return) out$signals <- signal
  if ("deletions" %in% return && !is.null(loo)) out$deletions <- loo
  if ("pseudo_values" %in% return && !is.null(loo)) {
    pv <- loo
    for (j in seq_len(m)) {
      full <- matrix(components[j, ], n, 2L, byrow = TRUE)
      pv[, , j] <- full + (n - 1) * (full - loo[, , j])
    }
    out$pseudo_values <- pv
  }
  out$.model <- c(list(
    x_mat = x_mat, signal = signal, x_cols = resolved_x_cols,
    k = k, s = s, method = method, normalize = normalize,
    denominator_tol = denominator_tol, row_ids = ids,
    geometry = list(ranks = wt$ranks, ord = wt$ord, dx = wt$dx)
  ), geometry)
  if (tuning_mode == "grid") {
    out$tuning_grid <- .affine_ratio_tuning_grid(out, signal, wt$ranks, model_weights,
      active_covariance, jackknife, denominator_tol, resolved$tuning)
  }
  class(out) <- "affine_ratio_result"
  out
}

# Common paired calculation for the default order and cached grid candidates.
.affine_ratio_evaluate <- function(signal, W, ranks, rw, active, jackknife, tol) {
  n <- nrow(signal)
  components <- t(.compensated_crossprod(signal, W))
  colnames(components) <- c("numerator", "denominator")
  rownames(components) <- NULL
  neff <- as.numeric(1 / colSums(W^2))
  theta <- components[, 1L] / components[, 2L]
  status <- .affine_ratio_status(components, theta, signal[, 2L], tol)
  theta[status$ratio_status != "ok"] <- NA_real_
  loo <- NULL
  covariance <- list(se = NULL, component_vcov = NULL, theta_vcov = NULL, theta_omega = NULL)
  if (active) {
    jk <- .dnn_loo_centered(signal, ranks, .rank_weight_table(rw, n - 1L),
      tab_full = .rank_weight_table(rw, n), rw = if (jackknife == "fast") NULL else rw)
    loo <- jk$loo
    covariance <- .affine_ratio_covariance(loo, components, theta, neff)
    status$covariance_status <- covariance$status
    status$wald_available <- is.finite(theta) & is.finite(covariance$se) & covariance$se > 0
    loo <- .dnn_restore_loo(loo, jk$offset)
  }
  c(list(theta = as.numeric(theta)), covariance[c("se", "component_vcov", "theta_vcov", "theta_omega")],
    list(components = components, neff = neff, status = status, deletions = loo))
}

.affine_ratio_status <- function(components, theta, denominator, tol) {
  d <- components[, 2L]
  finite_components <- apply(components, 1L, function(z) all(is.finite(z)))
  ratio_status <- rep("ok", nrow(components))
  ratio_status[!finite_components] <- "component_overflow"
  ratio_status[finite_components & d == 0] <- "zero_denominator"
  ratio_status[finite_components & d != 0 & !is.finite(theta)] <- "ratio_overflow"
  den_scale <- max(abs(denominator))
  near <- if (tol > 0 && den_scale > 0) {
    is.finite(d) & d != 0 & abs(d) / den_scale <= tol
  } else {
    rep(FALSE, length(d))
  }
  data.frame(
    target_id = seq_len(nrow(components)), ratio_status = ratio_status,
    negative_denominator = is.finite(d) & d < 0,
    near_zero_denominator = near, covariance_status = "not_computed",
    wald_available = FALSE, stringsAsFactors = FALSE
  )
}

.affine_ratio_covariance <- function(loo, components, theta, neff) {
  n <- dim(loo)[1L]
  m <- nrow(components)
  paired <- .dnn_jackknife_vcov_from_loo(loo, q = 2L, n = n, m = m,
    allow_failure = TRUE)$vcov
  theta_vcov <- theta_omega <- vector("list", m)
  se <- rep(NA_real_, m)
  status <- rep("ok", m)
  for (j in seq_len(m)) {
    dimnames(paired[[j]]) <- list(colnames(components), colnames(components))
    theta_vcov[[j]] <- theta_omega[[j]] <- matrix(NA_real_, 1L, 1L, dimnames = list("theta", "theta"))
    centered <- sweep(loo[, , j], 2L, colMeans(loo[, , j]), "-")
    varying <- colSums(centered != 0) > 0
    if (any(!is.finite(paired[[j]])) || any(!is.finite(centered)) ||
      any(varying & diag(paired[[j]]) == 0)) {
      paired[[j]][, ] <- NA_real_
      status[j] <- "numerical_failure"
      next
    }
    if (!is.finite(theta[j])) next
    # Divide each centered component before combining to avoid denominator squaring.
    g <- theta[j] * (centered[, 2L] / components[j, 2L]) -
      centered[, 1L] / components[j, 2L]
    scale <- max(abs(g))
    if (!is.finite(scale)) {
      status[j] <- "numerical_failure"
      next
    }
    v <- .covariance_matrix(matrix(g), (n - 1) / n, center = FALSE,
      allow_failure = TRUE)[1L, 1L]
    sd <- sqrt(v)
    omega <- neff[j] * v
    if (!is.finite(v) || !is.finite(omega) || (scale > 0 && v == 0)) {
      status[j] <- "numerical_failure"
      next
    }
    theta_vcov[[j]][1L, 1L] <- v
    theta_omega[[j]][1L, 1L] <- omega
    se[j] <- sd
    if (scale == 0) status[j] <- "zero_score_variance"
  }
  list(
    se = se, component_vcov = paired, theta_vcov = theta_vcov,
    theta_omega = theta_omega, status = status
  )
}
