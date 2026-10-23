# DNN-based local GMM wrapper.

.gmm_prepare_return <- function(return) {
  valid <- c(
    "W", "ranks", "ord", "dx", "g_obs", "omega", "jacobian",
    "ess", "neighbors", "tuning"
  )
  if (is.null(return)) {
    return <- character(0)
  } else {
    return <- unique(as.character(return))
  }
  .warn_unknown_return(return, valid = valid)
  intersect(return, valid)
}

.gmm_moment_matrix <- function(g, n, arg = "moment_obs") {
  if (is.null(dim(g))) {
    g <- .validate_vector_nrow(g, n, arg = arg, ref = "X", finite = TRUE)
    return(matrix(g, ncol = 1L))
  }
  g_mat <- as.matrix(g)
  if (nrow(g_mat) != n) {
    .abort_bad_arg(
      arg,
      "must return one row per observation",
      expected = sprintf("nrow=%d", n),
      actual = sprintf("nrow=%d", nrow(g_mat)),
      quoted = FALSE
    )
  }
  if (!is.numeric(g_mat)) {
    .abort_bad_arg(arg, "must be numeric", quoted = FALSE)
  }
  .abort_if_na_or_nonfinite(g_mat, arg = arg, quoted = FALSE)
  g_mat
}

.gmm_eval_moment <- function(model, sample, theta, ctx, n) {
  g <- model$moment_obs(sample = sample, theta = theta, ctx = ctx)
  .gmm_moment_matrix(g, n = n, arg = "model$moment_obs()")
}

.gmm_gbar <- function(model, sample, theta, ctx, w, n) {
  g <- .gmm_eval_moment(model, sample, theta, ctx, n)
  as.numeric(.compensated_crossprod(g, w))
}

.gmm_objective <- function(theta, model, sample, ctx, w, n, omega_inv = NULL) {
  val <- tryCatch({
    gbar <- .gmm_gbar(model, sample, theta, ctx, w, n)
    if (is.null(omega_inv)) {
      sum(gbar^2)
    } else {
      as.numeric(t(gbar) %*% omega_inv %*% gbar)
    }
  }, error = function(e) Inf)
  if (!is.finite(val)) .Machine$double.xmax else as.numeric(val)
}

.gmm_safe_solve <- function(A, b = NULL) {
  if (any(!is.finite(A)) || (!is.null(b) && any(!is.finite(b)))) {
    return(NULL)
  }
  tryCatch(
    if (is.null(b)) solve(A) else solve(A, b),
    error = function(e) NULL
  )
}

.gmm_condition_number <- function(A) {
  vals <- tryCatch(eigen(A, symmetric = TRUE, only.values = TRUE)$values,
                   error = function(e) NA_real_)
  if (length(vals) == 1L && is.na(vals)) {
    return(NA_real_)
  }
  vals <- abs(vals)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) {
    return(NA_real_)
  }
  max(vals) / max(min(vals), .Machine$double.eps)
}

.gmm_regularize_omega <- function(Omega, omega_ridge) {
  Omega <- .covariance_symmetrize(Omega)
  scale <- max(1, mean(abs(diag(Omega))))
  ridge <- as.numeric(omega_ridge) * scale
  if (ridge > 0) {
    Omega <- Omega + diag(ridge, nrow(Omega))
  }
  list(Omega = Omega, ridge = ridge, condition = .gmm_condition_number(Omega))
}

.gmm_eval_jacobian <- function(model, sample, theta, ctx, w, n, q, p,
                               mode, fd_step) {
  use_analytic <- mode %in% c("auto", "analytic") && !is.null(model$jacobian_obs)
  if (use_analytic) {
    jac <- model$jacobian_obs(sample = sample, theta = theta, ctx = ctx)
    if (p == 1L && is.null(dim(jac))) {
      jac <- .validate_vector_nrow(jac, n, arg = "model$jacobian_obs()", ref = "X", finite = TRUE)
      jac <- matrix(jac, ncol = 1L)
    }
    if (p == 1L && length(dim(jac)) == 2L) {
      jac_mat <- as.matrix(jac)
      if (nrow(jac_mat) != n || ncol(jac_mat) != q) {
        .abort_bad_arg(
          "model$jacobian_obs()",
          "must return an n x q matrix for scalar theta",
          expected = sprintf("%d x %d", n, q),
          actual = sprintf("%d x %d", nrow(jac_mat), ncol(jac_mat)),
          quoted = FALSE
        )
      }
      if (!is.numeric(jac_mat)) {
        .abort_bad_arg("model$jacobian_obs()", "must be numeric", quoted = FALSE)
      }
      .abort_if_na_or_nonfinite(jac_mat, arg = "model$jacobian_obs()", quoted = FALSE)
      return(matrix(as.numeric(.compensated_crossprod(jac_mat, w)), nrow = q, ncol = 1L))
    }
    if (length(dim(jac)) != 3L) {
      .abort_bad_arg(
        "model$jacobian_obs()",
        "must return an n x q x p array",
        quoted = FALSE
      )
    }
    if (dim(jac)[1L] != n || dim(jac)[2L] != q || dim(jac)[3L] != p) {
      .abort_bad_arg(
        "model$jacobian_obs()",
        "must return an n x q x p array",
        expected = sprintf("%d x %d x %d", n, q, p),
        actual = paste(dim(jac), collapse = " x "),
        quoted = FALSE
      )
    }
    if (!is.numeric(jac)) {
      .abort_bad_arg("model$jacobian_obs()", "must be numeric", quoted = FALSE)
    }
    .abort_if_na_or_nonfinite(jac, arg = "model$jacobian_obs()", quoted = FALSE)
    D <- matrix(NA_real_, nrow = q, ncol = p)
    for (j in seq_len(p)) {
      jac_j <- matrix(jac[, , j], nrow = n, ncol = q)
      D[, j] <- as.numeric(.compensated_crossprod(jac_j, w))
    }
    return(D)
  }

  if (identical(mode, "analytic")) {
    .abort_bad_arg("jacobian", "`analytic` requested but model has no `jacobian_obs()`")
  }

  D <- matrix(NA_real_, nrow = q, ncol = p)
  theta <- as.numeric(theta)
  for (j in seq_len(p)) {
    h <- fd_step * max(1, abs(theta[j]))
    plus <- theta
    minus <- theta
    plus[j] <- plus[j] + h
    minus[j] <- minus[j] - h
    D[, j] <- (.gmm_gbar(model, sample, plus, ctx, w, n) -
      .gmm_gbar(model, sample, minus, ctx, w, n)) / (2 * h)
  }
  D
}

.gmm_validate_start_vector <- function(start, p, arg = "start") {
  if (!is.numeric(start)) {
    .abort_bad_arg(arg, "must be numeric")
  }
  start_vec <- as.numeric(start)
  if (length(start_vec) != p || anyNA(start_vec) || any(!is.finite(start_vec))) {
    .abort_bad_arg(
      arg,
      "must resolve to a finite numeric vector with length equal to length(coef_names)",
      expected = sprintf("length=%d", p),
      actual = sprintf("length=%d", length(start_vec))
    )
  }
  start_vec
}

.gmm_start_for_target <- function(model, sample, ctx, w, x0, p,
                                  user_start = NULL) {
  if (!is.null(user_start)) {
    return(.gmm_validate_start_vector(user_start, p = p, arg = "start"))
  }
  if (!is.null(model$start)) {
    start <- model$start(sample = sample, ctx = ctx, w = w, x0 = x0)
    return(.gmm_validate_start_vector(start, p = p, arg = "model$start()"))
  }
  rep(0, p)
}

.gmm_fit_one_target <- function(model, sample, ctx, w, x0, n, q, p,
                                start, optimizer, jacobian, omega_ridge,
                                fd_step, se, need_jacobian,
                                moment_names, coef_names) {
  neff <- 1 / sum(w^2)
  theta_start <- .gmm_start_for_target(model, sample, ctx, w, x0, p, start)
  optim_control <- if (identical(optimizer, "L-BFGS-B")) {
    list(maxit = 1000L, factr = 1e-12 / .Machine$double.eps)
  } else {
    list(maxit = 1000L, reltol = 1e-12)
  }

  opt1 <- stats::optim(
    par = theta_start,
    fn = .gmm_objective,
    model = model,
    sample = sample,
    ctx = ctx,
    w = w,
    n = n,
    omega_inv = NULL,
    method = optimizer,
    control = optim_control
  )
  theta_first <- as.numeric(opt1$par)
  g1 <- .gmm_eval_moment(model, sample, theta_first, ctx, n)
  Omega_raw <- .covariance_matrix(g1 * sqrt(w), center = FALSE)
  omega_reg <- .gmm_regularize_omega(Omega_raw, omega_ridge = omega_ridge)
  Omega <- omega_reg$Omega
  rownames(Omega) <- colnames(Omega) <- moment_names
  Omega_inv <- .gmm_safe_solve(Omega)
  if (is.null(Omega_inv) || any(!is.finite(Omega_inv))) {
    .abort_bad_state("Could not invert the local GMM moment covariance matrix.")
  }

  opt2 <- stats::optim(
    par = theta_first,
    fn = .gmm_objective,
    model = model,
    sample = sample,
    ctx = ctx,
    w = w,
    n = n,
    omega_inv = Omega_inv,
    method = optimizer,
    control = optim_control
  )
  theta_hat <- as.numeric(opt2$par)
  g2 <- .gmm_eval_moment(model, sample, theta_hat, ctx, n)
  gbar2 <- as.numeric(.compensated_crossprod(g2, w))

  D_hat <- if (need_jacobian) {
    .gmm_eval_jacobian(
      model = model, sample = sample, theta = theta_hat, ctx = ctx,
      w = w, n = n, q = q, p = p, mode = jacobian, fd_step = fd_step
    )
  } else {
    NULL
  }
  if (!is.null(D_hat)) {
    rownames(D_hat) <- moment_names
    colnames(D_hat) <- coef_names
  }

  vc <- NULL
  se_hat <- rep(NA_real_, p)
  if (isTRUE(se)) {
    info <- t(D_hat) %*% Omega_inv %*% D_hat
    rownames(info) <- colnames(info) <- coef_names
    info_inv <- .gmm_safe_solve(info)
    outside_range <- any(!is.finite(info)) ||
      (!is.null(info_inv) && any(!is.finite(info_inv)))
    if (outside_range) {
      warning(paste0(
        "gmm_dnn: local GMM information matrix or its inverse is outside the numeric range; ",
        "coefficient SEs set to NA. Rescale model moments or the parameterization."
      ), call. = FALSE)
    } else if (!is.null(info_inv)) {
      vc <- info_inv / neff
      rownames(vc) <- colnames(vc) <- coef_names
      se_hat <- sqrt(pmax(diag(vc), 0))
    }
  }

  df <- as.integer(q - p)
  J <- pval <- NA_real_
  if (df > 0L) {
    J <- as.numeric(neff * t(gbar2) %*% Omega_inv %*% gbar2)
    pval <- stats::pchisq(J, df = df, lower.tail = FALSE)
  } else {
    df <- NA_integer_
  }

  list(
    theta = theta_hat,
    theta_first = theta_first,
    se = se_hat,
    vcov = vc,
    omega = Omega,
    jacobian = D_hat,
    g_obs = g2,
    J = J,
    df = df,
    pval_J = pval,
    neff = neff,
    convergence = list(
      first = opt1$convergence,
      second = opt2$convergence,
      value_first = opt1$value,
      value_second = opt2$value,
      omega_ridge = omega_reg$ridge,
      omega_condition = omega_reg$condition,
      message = NA_character_
    )
  )
}

.gmm_validate_start <- function(start, p, m) {
  if (is.null(start)) {
    return(NULL)
  }
  if (is.null(dim(start))) {
    start_vec <- .gmm_validate_start_vector(start, p = p, arg = "start")
    return(matrix(start_vec, nrow = p, ncol = m))
  }

  start_mat <- as.matrix(start)
  if (!is.numeric(start_mat)) {
    .abort_bad_arg("start", "must be numeric")
  }
  if (nrow(start_mat) != p || ncol(start_mat) != m) {
    .abort_bad_arg(
      "start",
      "matrix form must have one row per coefficient and one column per target",
      expected = sprintf("%d x %d", p, m),
      actual = sprintf("%d x %d", nrow(start_mat), ncol(start_mat))
    )
  }
  .abort_if_na_or_nonfinite(start_mat, arg = "start")
  start_mat
}

.gmm_probe_moments <- function(model, sample, ctx, n, theta) {
  .gmm_eval_moment(model, sample, theta, ctx, n)
}

#' Local GMM estimation via DNN rank weights
#'
#' Estimates a local parameter by solving a two-step GMM problem after
#' localizing observations with DNN rank weights.  The first step uses an
#' identity moment weighting matrix; the second step uses the local moment
#' covariance estimated at the first-step solution.
#' This local IV-GMM interface is experimental and may change.
#'
#' @details
#' The documented sampling framework uses independent observations; cluster
#' and survey-design inference are not supplied. Observations with exactly
#' equal computed distance to a target share their block's average rank weight,
#' even when their covariates differ. Groups are target-specific and use exact
#' represented equality, without a tolerance. Ordering uses distance then input
#' row index; outcomes and nuisances do not define groups. Every jackknife
#' deletion averages the surviving blocks with the sample-size-n-1 rule.
#' A warning is emitted once per public call with distance ties: the accompanying
#' theory assumes a continuous distance distribution and does not establish
#' inference under persistent ties. Earlier retained fit schemas require refitting.
#'
#' Two interfaces are supported.
#'
#' 1. Data-frame interface: pass `data`, `x_cols`, and a `gmm_model`.
#' 2. Matrix interface: pass `X`, one of `X0` or `x0`, a `gmm_model`, and
#'    `model_data` when the model needs observation-level variables.
#'
#' Nuisance quantities are not fit internally.  They must be supplied to the
#' model through `data` or `model_data`.
#'
#' The resolved distance specification and its fitted training context are
#' retained for weight reconstruction. Replacing a registered metric after
#' fitting does not change the fit's geometry. Custom callbacks with mutable
#' external state retain their own reproducibility obligations.
#'
#' Exact DNN localization uses the full rank permutation at each target rather
#' than an `s`-nearest-neighbor truncation. For exact DNN fits, the
#' unaveraged rank-table support edge is rank `n - s + k`, not rank `s`.
#'
#' Averaging can extend nonzero observation weights through the end of a
#' distance block that straddles the nominal rank-table edge. Effective sample
#' sizes and weight diagnostics use the resulting observation weights.
#'
#' @param data Optional data frame containing covariates and model variables.
#'   Column names must be non-missing, non-empty, and unique.
#' @param x_cols Names or indices of covariate columns in `data`.
#' @param X Numeric `n x p` covariate matrix for the matrix interface.
#' @param X0 Numeric matrix of target points; mutually exclusive with `x0`.
#' @param x0 Numeric vector for one target point.
#' @param model A [new_gmm_model()] object.
#' @param model_data Optional model data used with the matrix interface.
#' @param distance Distance name or specification.
#' @param method One of `"exact"` or `"asymp"` for DNN weights.
#' @param k Integer DNN neighbour count inside each size-`s` subsample. The
#'   paper theory uses `k = 1` (the default).
#' @param s Integer kernel order (undersmoothing parameter), or `NULL`. When
#'   both `s` and `tuning` are `NULL` (the default), `s` is resolved from the
#'   paper undersmoothing rule [tuning_dnn_paper()],
#'   \eqn{s = \lfloor 0.4 n^{(k_x+1)/(k_x+2)} \rfloor}. Note that small fixed
#'   values such as `s = 1` (or `k = s`) yield uniform weights and hence no
#'   localization; a warning is emitted in that case.
#' @param tuning Optional kernel-order tuning specification.  When supplied,
#'   `s` must be omitted.
#' @param normalize Logical; normalize weight columns to sum to one.
#' @param start Optional numeric vector of length `p_theta`, or `p_theta x m`
#'   matrix, used as optimizer starts. When omitted, the model's `start()`
#'   callback is evaluated separately at every target, or zero is used when
#'   the model has no callback. Targets do not inherit previous targets' starts
#'   or solutions. An explicit vector is reused at all targets; matrix column
#'   `j` supplies target `j`. Either explicit form takes precedence over the
#'   model callback.
#' @param optimizer Optimizer method passed to [stats::optim()]. Both stages
#'   allow 1000 iterations. BFGS and Nelder-Mead use `reltol = 1e-12`;
#'   L-BFGS-B uses `factr = 1e-12 / .Machine$double.eps` and the default
#'   projected-gradient tolerance. Convergence statuses remain available in
#'   the returned object.
#' @param jacobian One of `"auto"`, `"analytic"`, or `"fd"`.
#' @param omega_ridge Non-negative ridge multiplier for local moment covariance.
#' @param fd_step Positive finite finite-difference step multiplier.
#' @param se Logical; compute sandwich standard errors.
#' @param ctx Optional context passed to distance and ranking components.
#' @param return Optional extras:
#'   `c("W", "ranks", "ord", "dx", "g_obs", "omega", "jacobian", "ess",
#'   "neighbors", "tuning")`. `neighbors` is a `k x m` matrix of nearest
#'   training-observation indices, not the full DNN weight support.
#'
#' @return A `gmm_result` object with estimates, diagnostics, metadata, the
#'   `m x p` evaluation-point matrix `x0` (with covariate column names
#'   propagated from the training data when available), and requested extras.
#' @examples
#' dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 120, kx = 2, q = 2, seed = 1)
#' model <- gmm_model_dr_iv(
#'   z_cols = c("Z1", "Z2"),
#'   pi_cols = c("piZ1_oracle", "piZ2_oracle"),
#'   muY_col = "muY_oracle",
#'   muD_col = "muD_oracle"
#' )
#' fit <- gmm_dnn(
#'   data = dgp$data,
#'   x_cols = c("X1", "X2"),
#'   x0 = dgp$x0,
#'   model = model,
#'   k = 1, s = 10,
#'   return = "ess"
#' )
#' fit$theta
#' @export
gmm_dnn <- function(
  data = NULL,
  x_cols = NULL,
  X = NULL,
  X0 = NULL,
  x0 = NULL,
  model,
  model_data = NULL,
  distance = "l2",
  method = c("exact", "asymp"),
  k = 1L,
  s = NULL,
  tuning = NULL,
  normalize = TRUE,
  start = NULL,
  optimizer = c("BFGS", "Nelder-Mead", "L-BFGS-B"),
  jacobian = c("auto", "analytic", "fd"),
  omega_ridge = 1e-8,
  fd_step = 1e-5,
  se = TRUE,
  ctx = NULL,
  return = c("W")
) {
  method <- match.arg(method)
  optimizer <- match.arg(optimizer)
  jacobian <- match.arg(jacobian)
  s_arg <- s
  k <- .validate_scalar_int(k, "k")
  return <- .gmm_prepare_return(return)

  model <- .gmm_validate_model_object(model)
  if (!is.logical(normalize) || length(normalize) != 1L || is.na(normalize)) {
    .abort_bad_arg("normalize", "must be TRUE or FALSE")
  }
  if (!is.logical(se) || length(se) != 1L || is.na(se)) {
    .abort_bad_arg("se", "must be TRUE or FALSE")
  }
  if (!is.numeric(omega_ridge) || length(omega_ridge) != 1L ||
      !is.finite(omega_ridge) || omega_ridge < 0) {
    .abort_bad_arg("omega_ridge", "must be a non-negative finite numeric scalar")
  }
  if (!is.numeric(fd_step) || length(fd_step) != 1L ||
      !is.finite(fd_step) || fd_step <= 0) {
    .abort_bad_arg("fd_step", "must be a positive finite numeric scalar")
  }

  if (!is.null(data)) {
    if (!is.null(X) || !is.null(model_data)) {
      .abort_bad_arg(
        "data/X/model_data",
        "must use either the data-frame interface or the matrix interface, not both",
        quoted = FALSE
      )
    }
    if (!is.data.frame(data)) {
      .abort_bad_arg("data", "must be a data.frame")
    }
    x_mat <- .dnn_col_mat(x_cols, data = data, arg = "x_cols", out_arg = "X")
    sample <- data
  } else {
    if (!is.null(x_cols)) {
      .abort_bad_arg("data", "must be supplied when column selectors are used")
    }
    if (is.null(X)) {
      .abort_bad_arg(
        "X",
        "must be provided for the matrix interface when `data` is NULL",
        expected = "non-NULL `X`, or use the `data=` interface"
      )
    }
    x_mat <- .validate_matrix_nonempty(X, arg = "X")
    if (is.null(model_data)) {
      sample <- data.frame(.row_id = seq_len(nrow(x_mat)))
    } else {
      sample <- model_data
      if (is.data.frame(sample) && nrow(sample) != nrow(x_mat)) {
        .abort_bad_arg(
          "nrow(model_data)",
          "must equal nrow(X)",
          expected = sprintf("nrow(X)=%d", nrow(x_mat)),
          actual = sprintf("nrow(model_data)=%d", nrow(sample)),
          quoted = FALSE
        )
      }
    }
  }

  n <- nrow(x_mat)
  p_x <- ncol(x_mat)
  p_theta <- length(model$coef_names)

  s_resolved <- .dnn_resolve_s(
    s = s_arg,
    tuning = tuning,
    n = n,
    p = p_x,
    k = k,
    se = FALSE,
    jackknife = "none",
    ctx = ctx
  )
  s <- s_resolved$s
  tuning_payload <- s_resolved$tuning
  tuning_used <- isTRUE(s_resolved$tuning_used)

  x0_mat <- .dnn_resolve_targets(x_mat, X0 = X0, x0 = x0)
  m_out <- nrow(x0_mat)
  start_mat <- .gmm_validate_start(start, p = p_theta, m = m_out)

  model_ctx <- model$setup(sample = sample, n = n)
  weight_return <- return
  if ("neighbors" %in% return && !any(c("ranks", "ord") %in% return)) {
    weight_return <- unique(c(weight_return, "ranks"))
  }
  .dnn_distance_scope()
  geometry <- .dnn_fitted_geometry(distance, x_mat, ctx)
  weight_fit <- .dnn_weights_with_extras(
    x_mat = x_mat,
    x0_mat = x0_mat,
    distance = geometry$distance,
    method = method,
    k = k,
    s = s,
    normalize = normalize,
    ctx = geometry$ctx,
    return = weight_return
  )
  wt <- weight_fit$wt
  W <- weight_fit$W

  # Find moment dimensions at an actual target start, including models whose
  # moments are not defined at zero. Cache every attempted start or failure;
  # an earlier target must not prevent a later valid target from supplying q.
  cached_starts <- vector("list", m_out)
  attempted_start <- rep(FALSE, m_out)
  probe <- NULL
  for (j in seq_len(m_out)) {
    attempted_start[j] <- TRUE
    start_j <- tryCatch(.gmm_start_for_target(
      model, sample, model_ctx, as.numeric(W[, j]), x0_mat[j, ], p_theta,
      if (is.null(start_mat)) NULL else start_mat[, j]), error = function(e) e)
    cached_starts[j] <- list(start_j)
    if (inherits(start_j, "error")) next
    candidate_probe <- tryCatch(.gmm_probe_moments(
      model, sample, model_ctx, n = n, theta = start_j), error = function(e) e)
    if (inherits(candidate_probe, "error")) {
      cached_starts[j] <- list(candidate_probe)
    } else {
      probe <- candidate_probe
      break
    }
  }
  # Preserve all-target warning/NA results when every start fails but the
  # former zero probe can still establish the model's moment dimensions.
  if (is.null(probe)) {
    probe <- .gmm_probe_moments(model, sample, model_ctx, n = n, theta = rep(0, p_theta))
  }
  q <- ncol(probe)
  moment_names <- model$moment_names
  if (is.null(moment_names)) {
    moment_names <- paste0("moment_", seq_len(q))
  }
  if (length(moment_names) != q) {
    .abort_bad_arg(
      "moment_names",
      "length must match the number of model moments",
      expected = sprintf("length=%d", q),
      actual = sprintf("length=%d", length(moment_names))
    )
  }

  theta <- matrix(NA_real_, nrow = p_theta, ncol = m_out)
  theta_first <- matrix(NA_real_, nrow = p_theta, ncol = m_out)
  se_mat <- if (isTRUE(se)) matrix(NA_real_, nrow = p_theta, ncol = m_out) else NULL
  J <- pval_J <- neff <- rep(NA_real_, m_out)
  df <- rep(NA_integer_, m_out)
  vcov_list <- vector("list", m_out)
  omega_list <- vector("list", m_out)
  jac_list <- vector("list", m_out)
  g_list <- vector("list", m_out)
  conv <- vector("list", m_out)

  for (j in seq_len(m_out)) {
    w <- as.numeric(W[, j])
    start_j <- if (attempted_start[j]) cached_starts[[j]] else if (!is.null(start_mat)) start_mat[, j] else NULL

    fit_j <- tryCatch({
      if (inherits(start_j, "error")) stop(start_j)
      .gmm_fit_one_target(
        model = model,
        sample = sample,
        ctx = model_ctx,
        w = w,
        x0 = x0_mat[j, ],
        n = n,
        q = q,
        p = p_theta,
        start = start_j,
        optimizer = optimizer,
        jacobian = jacobian,
        omega_ridge = omega_ridge,
        fd_step = fd_step,
        se = se,
        need_jacobian = isTRUE(se) || "jacobian" %in% return,
        moment_names = moment_names,
        coef_names = model$coef_names
      )
    },
      error = function(e) {
        warning(sprintf(
          "gmm_dnn: target %d failed: %s",
          j, conditionMessage(e)
        ), call. = FALSE)
        NULL
      }
    )

    if (is.null(fit_j)) {
      conv[[j]] <- list(
        first = NA_integer_,
        second = NA_integer_,
        value_first = NA_real_,
        value_second = NA_real_,
        omega_ridge = NA_real_,
        omega_condition = NA_real_,
        message = "target failed"
      )
      next
    }

    theta[, j] <- fit_j$theta
    theta_first[, j] <- fit_j$theta_first
    if (isTRUE(se)) se_mat[, j] <- fit_j$se
    J[j] <- fit_j$J
    df[j] <- fit_j$df
    pval_J[j] <- fit_j$pval_J
    neff[j] <- fit_j$neff
    vcov_list[j] <- list(fit_j$vcov)
    omega_list[j] <- list(fit_j$omega)
    jac_list[j] <- list(fit_j$jacobian)
    g_list[j] <- list(fit_j$g_obs)
    conv[[j]] <- fit_j$convergence
  }

  target_names <- paste0("target_", seq_len(m_out))
  rownames(theta) <- rownames(theta_first) <- model$coef_names
  colnames(theta) <- colnames(theta_first) <- target_names
  if (isTRUE(se)) {
    rownames(se_mat) <- model$coef_names
    colnames(se_mat) <- target_names
  }
  names(J) <- target_names
  names(df) <- target_names
  names(pval_J) <- target_names
  names(neff) <- target_names
  names(vcov_list) <- target_names
  names(omega_list) <- target_names
  names(jac_list) <- target_names
  names(g_list) <- target_names

  convergence <- do.call(rbind, lapply(seq_along(conv), function(j) {
    cj <- conv[[j]]
    data.frame(
      target = target_names[j],
      convergence_first = cj$first,
      convergence_second = cj$second,
      value_first = cj$value_first,
      value_second = cj$value_second,
      omega_ridge = cj$omega_ridge,
      omega_condition = cj$omega_condition,
      message = cj$message,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))

  out <- list(
    theta = theta,
    theta_first = theta_first,
    se = se_mat,
    J = J,
    df = df,
    pval_J = pval_J,
    neff = neff,
    convergence = convergence,
    x0 = .dnn_x0_with_names(x0_mat, x_mat),
    meta = list(
      n = n, m = m_out, p = p_x, p_theta = p_theta, q = q,
      k = as.integer(k), s = as.integer(s),
      method = method, distance = .dnn_distance_name(distance),
      tie_rule = .dnn_tie_policy(),
      model = model$name,
      jacobian = jacobian,
      optimizer = optimizer,
      se = isTRUE(se),
      coef_names = model$coef_names,
      moment_names = moment_names
    )
  )
  if (isTRUE(se)) out$vcov <- vcov_list
  if ("omega" %in% return) out$omega <- omega_list
  if ("jacobian" %in% return) out$jacobian <- jac_list
  if ("g_obs" %in% return) out$g_obs <- g_list
  if (tuning_used) out$meta$tuning <- tuning_payload
  out <- .dnn_attach_common_extras(
    out,
    return = return,
    wt = wt,
    W = W,
    k = k,
    tuning_payload = tuning_payload
  )
  out$.model <- c(list(
    x_mat = x_mat,
    x_cols = if (!is.null(data)) x_cols else NULL,
    method = method,
    k = k,
    s = s,
    normalize = normalize, ranks = wt$ranks, groups = wt$groups, dx = wt$dx
  ), geometry)

  class(out) <- "gmm_result"
  out
}
