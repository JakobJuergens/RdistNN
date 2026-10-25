# PLQR score arithmetic, density-backed inference and retained prediction state.

.plqr_dnn_score_matrix <- function(d_mat, mD_mat, g_vec, y_vec, theta_eval, tau,
                                    below_threshold = NULL) {
  if (is.null(below_threshold)) {
    below_threshold <- .plqr_below_threshold(d_mat, theta_eval, g_vec, y_vec)
  }
  resid <- as.numeric(tau) - as.numeric(below_threshold)
  score <- (d_mat - mD_mat) * resid
  if (any(!is.finite(score))) {
    .dnn_abort("numerical_range",
      "Computed PLQR score is outside the numeric range; rescale regressor and nuisance inputs."
    )
  }
  score
}

.plqr_dnn_jackknife_se_loo <- function(
  x_mat, x0_mat, score_mat, distance, method, k, s, normalize, ctx,
  jackknife = c("fast", "loo"),
  precomputed_ranks = NULL, precomputed_groups = NULL
) {
  jackknife <- .dnn_match_arg(jackknife)
  .dnn_jackknife_se_loo_matrix(
    x_mat = x_mat,
    x0_mat = x0_mat,
    signal = score_mat,
    distance = distance,
    method = method,
    k = k,
    s = s,
    normalize = normalize,
    ctx = ctx,
    signal_arg = "score",
    jackknife = jackknife,
    precomputed_ranks = precomputed_ranks,
    precomputed_groups = precomputed_groups
  )
}

.plqr_validate_density <- function(density, n, arg = "density") {
  density_vec <- .validate_vector_nrow(
    density,
    n,
    arg = arg,
    ref = "X",
    strict = TRUE,
    finite = TRUE
  )
  if (any(density_vec <= 0)) {
    .abort_bad_arg(arg, "must contain positive density values", quoted = FALSE)
  }
  density_vec
}

.plqr_safe_solve <- function(A) {
  if (any(!is.finite(A))) return(NULL)
  power <- 0
  if (!is.finite(norm(A, type = "1"))) {
    # Finite entries can have an overflowing norm inside LAPACK's condition
    # estimate. Binary normalization leaves the ordinary solver route intact.
    power <- .numeric_binary_parts(max(abs(A)))$power
    A <- A / 2^power
  }
  if (rcond(A) < .Machine$double.eps) return(NULL)
  if (power == 0) solve(A) else solve(A) / 2^power
}

.plqr_check_report_model <- function(object, n, m, invalid, require_fields,
                                    check_array, check_covariance) {
  model <- object$.model
  require_fields(c("d_mat", "mD_mat", "d_resid_mat", "below_threshold", "density",
    "theta_eval", "tau", "coef_names"), model, ".model$")
  .dnn_validate_retained(.plqr_validate_prediction_model(model))
  if (!is.matrix(model$d_mat) || ncol(model$d_mat) < 1L) {
    invalid(".model$d_mat", "must retain the canonical PLQR coefficient matrix")
  }
  q <- ncol(model$d_mat)
  for (field in c("d_mat", "mD_mat", "d_resid_mat")) {
    check_array(model[[field]], paste0(".model$", field), c(n, q))
  }
  check_array(model$theta_eval, ".model$theta_eval", q)
  check_array(model$below_threshold, ".model$below_threshold", n)
  if (any(!model$below_threshold %in% c(0, 1))) {
    invalid(".model$below_threshold", "must retain the fitted binary threshold indicators")
  }
  if (any(model$d_resid_mat != model$d_mat - model$mD_mat)) {
    invalid(".model$d_resid_mat", "is inconsistent with the retained regressor and nuisance matrices")
  }
  if (!is.null(dim(model$coef_names))) invalid(".model$coef_names", "must be a character vector")
  meta <- object$meta
  require_fields(c("p_d", "coef_names", "theta_source", "theta_se_method",
    "theta_meat_method", "moment_se_method", "jackknife"), meta, "meta$")
  if (!is.numeric(meta$p_d) || length(meta$p_d) != 1L || !is.finite(meta$p_d) || meta$p_d != q ||
      !identical(meta$coef_names, model$coef_names) || !identical(meta$theta_source, "supplied")) {
    invalid("meta", "must retain the fitted PLQR coefficient dimension, names and supplied evaluation source")
  }
  choices <- list(theta_se_method = c("none", "direct_sandwich"),
    theta_meat_method = c("none", "analytic", "jackknife"),
    moment_se_method = c("none", "jackknife"), jackknife = c("none", "fast", "loo"))
  for (field in names(choices)) {
    value <- meta[[field]]
    if (!is.character(value) || length(value) != 1L || is.na(value) || !value %in% choices[[field]]) {
      invalid(paste0("meta$", field), "must retain the fitted PLQR uncertainty method")
    }
  }
  coefficient_se <- meta$theta_se_method == "direct_sandwich"
  moment_se <- meta$moment_se_method == "jackknife"
  inference_fields <- c("theta_vcov", "theta_omega", "jacobian", "neff",
    "score_covariance", "moment_vcov")
  if (coefficient_se != !is.null(object$theta_se) || moment_se != !is.null(object$moment_se) ||
      (coefficient_se && (is.null(model$density) || meta$theta_meat_method == "none")) ||
      (!coefficient_se && meta$theta_meat_method != "none") ||
      (!coefficient_se && any(vapply(object[inference_fields], Negate(is.null), logical(1L)))) ||
      (moment_se && meta$jackknife == "none")) {
    invalid("meta/theta_se/moment_se", "contains inconsistent PLQR uncertainty availability")
  }
  require_fields("tau")
  check_array(object$tau, "tau", 1L)
  if (object$tau != model$tau || !identical(as.numeric(object$theta_eval), as.numeric(model$theta_eval))) {
    invalid("tau/theta_eval", "is inconsistent with the retained PLQR evaluation values")
  }
  if (coefficient_se) {
    require_fields(c("theta_vcov", "theta_omega", "jacobian", "neff"))
    check_covariance(object$theta_omega, "theta_omega", q)
    check_array(object$neff, "neff", m, nonnegative = TRUE)
    if (any(object$neff <= 0)) invalid("neff", "must retain positive effective sample sizes")
    if (!is.list(object$jacobian) || length(object$jacobian) != m) {
      invalid("jacobian", "must retain one PLQR Jacobian matrix per target")
    }
    for (j in seq_len(m)) {
      J <- object$jacobian[[j]]
      if (!is.matrix(J) || !is.numeric(J) || is.complex(J) || !identical(dim(J), c(q, q)) ||
          (any(!is.finite(J)) && any(!is.na(object$theta_se[, j])))) {
        invalid("jacobian", paste0("must retain aligned numeric bread matrices, ",
          "with unavailable uncertainty outside the numeric range"))
      }
    }
    meat <- if (meta$theta_meat_method == "analytic") "score_covariance" else "moment_vcov"
    require_fields(meat)
    check_covariance(object[[meat]], meat, q)
  }
  invisible(NULL)
}

.plqr_dimnames_matrix <- function(x, coef_names, target_names) {
  rownames(x) <- coef_names
  colnames(x) <- target_names
  x
}

# The direct-sandwich theorem convention assumes weights that sum to one per
# target; unnormalized weights would make coefficient SEs scale-dependent.
.plqr_validate_normalized_weights <- function(W, tol = 1e-8) {
  col_sums <- colSums(as.matrix(W))
  if (any(!is.finite(col_sums)) || any(abs(col_sums - 1) > tol)) {
    .dnn_abort("internal_invariant",
      "Direct PLQR sandwich requires normalized weight columns (each summing to 1)."
    , hint = "Report a reproducible example of this internal inconsistency.")
  }
  invisible(NULL)
}

# Density-backed bread of the PLQR coefficient sandwich: per-target Jacobian
# J(x) = sum_i w_i f_i tilde_D_i D_i' plus effective sample sizes.
.plqr_direct_bread <- function(d_mat, mD_mat, density_vec, W,
                               coef_names, target_names) {
  d_resid_mat <- d_mat - mD_mat
  m <- ncol(W)
  jacobian <- vector("list", m)
  for (j in seq_len(m)) {
    J <- crossprod(d_resid_mat, d_mat * as.numeric(density_vec * W[, j]))
    dimnames(J) <- list(coef_names, coef_names)
    jacobian[[j]] <- J
  }
  names(jacobian) <- target_names
  list(
    jacobian = jacobian,
    neff = stats::setNames(as.numeric(1 / colSums(W^2)), target_names)
  )
}

# Sandwich J^-1 meat J^-T per target. `meat_is_vcov = FALSE` treats the meat
# as an asymptotic score covariance (analytic route): omega = J^-1 S J^-T and
# vcov = omega / neff. `meat_is_vcov = TRUE` treats the meat as a direct
# estimate of Var(moment) (jackknife route): vcov = J^-1 meat J^-T and
# omega = vcov * neff.
.plqr_sandwich_from_meat <- function(jacobian, meat, neff,
                                     coef_names, target_names,
                                     meat_is_vcov) {
  m <- length(jacobian)
  p_d <- length(coef_names)
  theta_omega <- vector("list", m)
  theta_vcov <- vector("list", m)
  theta_se <- matrix(NA_real_, nrow = p_d, ncol = m)

  for (j in seq_len(m)) {
    outside_range <- any(!is.finite(jacobian[[j]]))
    J_inv <- if (outside_range) NULL else .plqr_safe_solve(jacobian[[j]])
    outside_range <- outside_range || (!is.null(J_inv) && any(!is.finite(J_inv)))
    covariance_outside_range <- FALSE
    if (!outside_range && !is.null(J_inv)) {
      sandwiched <- J_inv %*% meat[[j]] %*% t(J_inv)
      sandwiched <- .covariance_symmetrize(sandwiched)
      if (isTRUE(meat_is_vcov)) {
        vcov_j <- sandwiched
        omega_j <- sandwiched * neff[[j]]
      } else {
        omega_j <- sandwiched
        vcov_j <- sandwiched / neff[[j]]
      }
      covariance_outside_range <- any(!is.finite(omega_j)) || any(!is.finite(vcov_j))
      outside_range <- covariance_outside_range
    }
    if (outside_range || is.null(J_inv)) {
      message <- if (outside_range) {
        component <- if (covariance_outside_range) "sandwich covariance" else "Jacobian or its inverse"
        paste0(
          "plqr_dnn: direct PLQR ", component, " is outside the numeric range ",
          "at target %d; coefficient SEs set to NA. Rescale regressor and density inputs."
        )
      } else {
        "plqr_dnn: direct PLQR Jacobian is singular at target %d; coefficient SEs set to NA."
      }
      warning(sprintf(message, j), call. = FALSE)
      omega_j <- matrix(NA_real_, nrow = p_d, ncol = p_d)
      vcov_j <- matrix(NA_real_, nrow = p_d, ncol = p_d)
    } else {
      theta_se[, j] <- sqrt(pmax(diag(vcov_j), 0))
    }
    dimnames(omega_j) <- list(coef_names, coef_names)
    dimnames(vcov_j) <- list(coef_names, coef_names)
    theta_omega[[j]] <- omega_j
    theta_vcov[[j]] <- vcov_j
  }

  names(theta_omega) <- names(theta_vcov) <- target_names
  list(
    theta_se = .plqr_dimnames_matrix(theta_se, coef_names, target_names),
    theta_vcov = theta_vcov,
    theta_omega = theta_omega
  )
}

.plqr_direct_sandwich <- function(d_mat, mD_mat, density_vec, W, tau,
                                  coef_names, target_names) {
  .plqr_validate_normalized_weights(W)
  d_resid_mat <- d_mat - mD_mat
  m <- ncol(W)

  bread <- .plqr_direct_bread(
    d_mat = d_mat,
    mD_mat = mD_mat,
    density_vec = density_vec,
    W = W,
    coef_names = coef_names,
    target_names = target_names
  )

  score_covariance <- vector("list", m)
  for (j in seq_len(m)) {
    S <- .covariance_matrix(d_resid_mat * sqrt(as.numeric(W[, j])),
      factor = as.numeric(tau) * (1 - as.numeric(tau)), center = FALSE)
    dimnames(S) <- list(coef_names, coef_names)
    score_covariance[[j]] <- S
  }
  names(score_covariance) <- target_names

  sandwich <- .plqr_sandwich_from_meat(
    jacobian = bread$jacobian,
    meat = score_covariance,
    neff = bread$neff,
    coef_names = coef_names,
    target_names = target_names,
    meat_is_vcov = FALSE
  )

  c(
    sandwich,
    list(
      jacobian = bread$jacobian,
      score_covariance = score_covariance,
      neff = bread$neff
    )
  )
}

# Jackknife-meat coefficient sandwich: full delete-1 covariance of the
# localized vector moment as the meat, density-backed bread unchanged.
.plqr_jackknife_meat_sandwich <- function(
  x_mat, x0_mat, score_mat, d_mat, mD_mat, density_vec, W,
  distance, method, k, s, normalize, ctx, jackknife,
  precomputed_ranks, coef_names, target_names, precomputed_groups = NULL
) {
  .plqr_validate_normalized_weights(W)

  loo <- .dnn_jackknife_loo_array(
    x_mat = x_mat,
    x0_mat = x0_mat,
    signal = score_mat,
    distance = distance,
    method = method,
    k = k,
    s = s,
    normalize = normalize,
    ctx = ctx,
    signal_arg = "score",
    jackknife = jackknife,
    precomputed_ranks = precomputed_ranks,
    precomputed_groups = precomputed_groups
  )
  jk <- .dnn_jackknife_vcov_from_loo(
    loo = loo,
    q = ncol(score_mat),
    n = nrow(x_mat),
    m = ncol(W)
  )
  moment_vcov <- lapply(jk$vcov, function(v) {
    dimnames(v) <- list(coef_names, coef_names)
    v
  })
  names(moment_vcov) <- target_names
  moment_se <- .plqr_dimnames_matrix(jk$se, coef_names, target_names)

  bread <- .plqr_direct_bread(
    d_mat = d_mat,
    mD_mat = mD_mat,
    density_vec = density_vec,
    W = W,
    coef_names = coef_names,
    target_names = target_names
  )
  sandwich <- .plqr_sandwich_from_meat(
    jacobian = bread$jacobian,
    meat = moment_vcov,
    neff = bread$neff,
    coef_names = coef_names,
    target_names = target_names,
    meat_is_vcov = TRUE
  )

  c(
    sandwich,
    list(
      jacobian = bread$jacobian,
      neff = bread$neff,
      moment_vcov = moment_vcov,
      moment_se = moment_se
    )
  )
}

.plqr_validate_prediction_targets <- function(x, arg) {
  .validate_matrix(
    x,
    arg = arg,
    min_rows = 1L,
    min_cols = 1L,
    numeric = TRUE,
    finite = TRUE
  )
}

.plqr_validate_prediction_model <- function(model) {
  .dnn_check_tie_model(model)
  if (!is.list(model) || is.null(model$x_mat) || is.null(model$d_mat) ||
      is.null(model$mD_mat) || is.null(model$theta_eval) || is.null(model$tau)) {
    .dnn_abort("state_invalid", paste0(
      "This object does not contain fitted model state for prediction. ",
      "Refit with the current plqr_dnn() implementation."
    ))
  }

  x_mat <- .validate_matrix(
    model$x_mat,
    arg = "object$.model$x_mat",
    min_rows = 1L,
    min_cols = 1L,
    numeric = TRUE,
    finite = TRUE
  )
  n <- nrow(x_mat)
  d_mat <- .validate_matrix(
    model$d_mat,
    arg = "object$.model$d_mat",
    min_rows = 1L,
    min_cols = 1L,
    numeric = TRUE,
    finite = TRUE
  )
  mD_mat <- .validate_matrix(
    model$mD_mat,
    arg = "object$.model$mD_mat",
    min_rows = 1L,
    min_cols = 1L,
    numeric = TRUE,
    finite = TRUE
  )
  if (nrow(d_mat) != n || nrow(mD_mat) != n) {
    .abort_bad_arg(
      "object$.model$d_mat/object$.model$mD_mat",
      "must have row counts equal to nrow(object$.model$x_mat)",
      expected = sprintf("nrow(object$.model$x_mat)=%d", n),
      actual = sprintf(
        "nrow(object$.model$d_mat)=%d, nrow(object$.model$mD_mat)=%d",
        nrow(d_mat), nrow(mD_mat)
      ),
      quoted = FALSE
    )
  }
  if (ncol(d_mat) != ncol(mD_mat)) {
    .abort_bad_arg(
      "object$.model$mD_mat",
      "must have the same number of columns as object$.model$d_mat",
      expected = sprintf("ncol(object$.model$d_mat)=%d", ncol(d_mat)),
      actual = sprintf("ncol(object$.model$mD_mat)=%d", ncol(mD_mat)),
      quoted = FALSE
    )
  }

  theta_eval <- .validate_finite_numeric_vector(
    model$theta_eval,
    arg = "object$.model$theta_eval",
    min_length = 1L
  )
  if (length(theta_eval) != ncol(d_mat)) {
    .abort_bad_arg(
      "object$.model$theta_eval",
      "length must equal ncol(object$.model$d_mat)",
      expected = sprintf("ncol(object$.model$d_mat)=%d", ncol(d_mat)),
      actual = sprintf("length(object$.model$theta_eval)=%d", length(theta_eval)),
      quoted = FALSE
    )
  }
  tau <- .validate_scalar_prob(model$tau, arg = "object$.model$tau")
  k <- .validate_scalar_int(model$k, arg = "object$.model$k", min = 1L)
  s <- .validate_scalar_int(model$s, arg = "object$.model$s", min = 1L)
  if (k > s) {
    .abort_bad_arg(
      "object$.model$k",
      sprintf("must satisfy k <= s (got k=%d, s=%d)", k, s),
      quoted = FALSE
    )
  }
  if (s > n) {
    .abort_bad_arg(
      "object$.model$s",
      "must satisfy s <= n",
      expected = sprintf("s <= n=%d", n),
      actual = sprintf("s=%d", s),
      quoted = FALSE
    )
  }
  if (!is.character(model$method) || length(model$method) != 1L ||
      is.na(model$method) || !(model$method %in% c("exact", "asymp"))) {
    .abort_bad_arg("object$.model$method", "must be one of 'exact' or 'asymp'", quoted = FALSE)
  }
  if (!is.logical(model$normalize) || length(model$normalize) != 1L || is.na(model$normalize)) {
    .abort_bad_arg("object$.model$normalize", "must be TRUE or FALSE", quoted = FALSE)
  }

  density <- if (is.null(model$density)) {
    NULL
  } else {
    .plqr_validate_density(model$density, n, arg = "object$.model$density")
  }
  coef_names <- model$coef_names
  if (!is.character(coef_names) || length(coef_names) != ncol(d_mat) ||
      anyNA(coef_names) || any(!nzchar(coef_names))) {
    .abort_bad_arg(
      "object$.model$coef_names",
      "must contain one non-empty coefficient name per PLQR coefficient",
      expected = sprintf("length %d", ncol(d_mat)),
      quoted = FALSE
    )
  }

  model$x_mat <- x_mat
  model$d_mat <- d_mat
  model$mD_mat <- mD_mat
  model$density <- density
  model$theta_eval <- theta_eval
  model$tau <- tau
  model$k <- k
  model$s <- s
  model$method <- model$method[[1L]]
  model$normalize <- model$normalize[[1L]]
  model$coef_names <- coef_names
  model
}

.plqr_resolve_targets <- function(model, newdata = NULL, X0 = NULL, x0 = NULL) {
  if (!is.null(newdata) && (!is.null(X0) || !is.null(x0))) {
    .abort_bad_arg(
      "newdata/X0/x0",
      "must provide either `newdata` or `X0`/`x0`, not both",
      quoted = FALSE
    )
  }

  if (!is.null(newdata)) {
    if (is.data.frame(newdata) && !is.null(model$x_cols)) {
      x0_mat <- .dnn_col_mat(model$x_cols, data = newdata, arg = "x_cols", out_arg = "newdata")
    } else {
      x0_mat <- .plqr_validate_prediction_targets(newdata, arg = "newdata")
    }
    x0_mat <- .plqr_validate_prediction_targets(x0_mat, arg = "newdata")
    return(.validate_ncol_relation(model$x_mat, x0_mat, y_arg = "newdata", x_arg = "X", relation = "match"))
  }

  if (!is.null(X0) && !is.null(x0)) {
    .abort_bad_arg("X0/x0", "must provide either `X0` or `x0`, not both", quoted = FALSE)
  }
  if (is.null(X0)) {
    if (is.null(x0)) {
      .abort_bad_arg(
        "newdata/X0/x0",
        "must provide `newdata`, `X0`, or a single `x0`",
        quoted = FALSE
      )
    }
    x0_vec <- .validate_target_vector(
      x0,
      ncol(model$x_mat),
      x0_arg = "x0",
      x_arg = "X",
      strict = TRUE,
      finite = TRUE
    )
    return(matrix(x0_vec, nrow = 1L))
  }
  x0_mat <- .plqr_validate_prediction_targets(X0, arg = "X0")
  .validate_ncol_relation(model$x_mat, x0_mat, y_arg = "X0", x_arg = "X", relation = "equal")
}

.plqr_prediction_theta <- function(theta_eval, m, coef_names, target_names) {
  theta <- matrix(rep(theta_eval, times = m), nrow = length(theta_eval), ncol = m)
  .plqr_dimnames_matrix(theta, coef_names, target_names)
}
