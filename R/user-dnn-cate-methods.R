# Methods for cate_result objects.

#' Joint jackknife covariance across CATE reporting profiles
#'
#' @param object A fitted `cate_result` from [cate_dnn()].
#' @param ... Not used.
#' @return An `m` by `m` covariance matrix, with axes in `coef(object)` order.
#' @details Uses common deleted observations with frozen supplied signals.
#'   The diagonal agrees with the fit's squared standard errors. The `fast`
#'   route freezes fitted geometry; `loo` retains deletion-specific geometry
#'   refitting. Covariance is calculated on request from retained model state,
#'   independently of optional `return` extras, including after serialization.
#'   Fitted uncertainty and an active jackknife are required. Singular matrices
#'   and exact zero covariance are retained. Joint covariance supplies the
#'   arithmetic for prespecified contrasts; valid contrast inference requires
#'   the corresponding joint assumptions and does not supply simultaneous bands.
#' @examples
#' fit <- cate_dnn(X = matrix(1:7, ncol = 1),
#'   X0 = matrix(c(2, 5), ncol = 1),
#'   pseudo_outcome = c(2, 1, 4, 3, 6, 5, 7), s = 3,
#'   return = character())
#' vcov(fit)
#' @export
vcov.cate_result <- function(object, ...) {
  jackknife <- object$meta$jackknife
  if (is.null(object$se) || is.null(jackknife) ||
      !(jackknife %in% c("fast", "loo"))) {
    .abort_bad_state("Covariance is not available. Refit with `se = TRUE` and an active jackknife.")
  }
  model <- .cate_validate_prediction_model(object$.model)
  required <- if (identical(jackknife, "loo")) {
    c("inference_distance", "inference_ctx")
  } else {
    c("distance", "ctx")
  }
  if (!all(required %in% names(model))) {
    .abort_bad_state("This object lacks fitted geometry for covariance. Refit with the current cate_dnn().")
  }
  x0_mat <- .dnn_resolve_targets(model$x_mat, X0 = object$x0)
  loo <- .dnn_jackknife_loo_array(
    x_mat = model$x_mat, x0_mat = x0_mat, signal = model$psi,
    distance = model[[required[[1L]]]], ctx = model[[required[[2L]]]],
    method = model$method, k = model$k, s = model$s,
    normalize = model$normalize, signal_arg = "pseudo_outcome",
    jackknife = jackknife, precomputed_ranks = model$ranks
  )
  n <- nrow(model$x_mat)
  V <- .covariance_matrix(matrix(loo, nrow = n), (n - 1) / n)
  dimnames(V) <- rep(list(names(stats::coef(object))), 2L)
  V
}

.cate_validate_prediction_model <- function(model) {
  if (!is.list(model) || is.null(model$x_mat) || is.null(model$psi)) {
    .abort_bad_state(paste0(
      "This object does not contain fitted model state for prediction. ",
      "Refit with the current cate_dnn() implementation."
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
  psi <- .validate_vector_nrow(
    model$psi,
    n,
    arg = "object$.model$psi",
    ref = "object$.model$x_mat",
    strict = TRUE,
    finite = TRUE
  )
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

  model$x_mat <- x_mat
  model$psi <- psi
  model$k <- k
  model$s <- s
  model$method <- model$method[[1L]]
  model$normalize <- model$normalize[[1L]]
  model
}

#' Summarize a `cate_result` object
#'
#' @param object A `cate_result` object from [cate_dnn()].
#' @param level Confidence level used for interval columns.
#' @param ... Not used.
#'
#' @return An object of class `summary.cate_result` with a coefficient table.
#' @examples
#' set.seed(1)
#' n <- 30
#' x_mat <- matrix(rnorm(n * 2), n, 2)
#' out <- cate_dnn(X = x_mat, x0 = c(0, 0), pseudo_outcome = rnorm(n))
#' summary(out)
#' @export
summary.cate_result <- function(object, level = 0.95, ...) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    .abort_bad_arg("level", "must be a numeric scalar in (0,1)")
  }

  est <- as.numeric(object$tau)
  m <- length(est)
  se <- object$se
  has_se <- is.numeric(se) && length(se) == m && all(is.finite(se))

  z <- pval <- lwr <- upr <- rep(NA_real_, m)
  se_col <- rep(NA_real_, m)
  if (has_se) {
    se_col <- as.numeric(se)
    z <- est / se_col
    pval <- 2 * stats::pnorm(-abs(z))
    crit <- stats::qnorm(0.5 + level / 2)
    lwr <- est - crit * se_col
    upr <- est + crit * se_col
  }

  coef_tab <- data.frame(
    estimate = est,
    std.error = se_col,
    statistic = z,
    p.value = pval,
    conf.low = lwr,
    conf.high = upr,
    row.names = paste0("target_", seq_len(m))
  )

  out <- list(
    coefficients = coef_tab,
    level = as.numeric(level),
    has_se = has_se,
    meta = object$meta
  )
  class(out) <- "summary.cate_result"
  out
}

#' Print a `summary.cate_result` object
#'
#' @param x A `summary.cate_result` object.
#' @param digits Number of significant digits for printing.
#' @param ... Not used.
#'
#' @return `x`, invisibly.
#' @export
print.summary.cate_result <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("<summary.cate_result>\n", sep = "")
  .dnn_print_localization(x$meta)
  cat("  targets: ", nrow(x$coefficients), "\n", sep = "")
  if (!isTRUE(x$has_se)) {
    cat("  note   : standard errors are not available in this object.\n", sep = "")
  }
  print(signif(x$coefficients, digits = digits))
  invisible(x)
}

#' Confidence intervals for `cate_result`
#'
#' @param object A `cate_result` object from [cate_dnn()].
#' @param parm Optional integer vector selecting target indices.
#' @param level Confidence level.
#' @param ... Not used.
#'
#' @return A matrix with lower and upper confidence bounds.
#' @examples
#' set.seed(2)
#' n <- 25
#' x_mat <- matrix(rnorm(n * 2), n, 2)
#' out <- cate_dnn(X = x_mat, x0 = c(0, 0), pseudo_outcome = rnorm(n), se = TRUE)
#' confint(out)
#' @export
confint.cate_result <- function(object, parm = NULL, level = 0.95, ...) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    .abort_bad_arg("level", "must be a numeric scalar in (0,1)")
  }
  if (is.null(object$se)) {
    .abort_bad_state("Standard errors are not available. Refit with `se = TRUE`.")
  }

  m <- length(object$tau)
  idx <- if (is.null(parm)) {
    seq_len(m)
  } else {
    .validate_index_vector(
      parm,
      upper = m,
      arg = "parm",
      label = "target indices"
    )
  }

  est <- as.numeric(object$tau[idx])
  se <- as.numeric(object$se[idx])
  crit <- stats::qnorm(0.5 + level / 2)
  ci <- cbind(
    conf.low = est - crit * se,
    conf.high = est + crit * se
  )
  rownames(ci) <- paste0("target_", idx)
  ci
}

#' Predict from a `cate_result` object at new targets
#'
#' @param object A `cate_result` object from [cate_dnn()].
#' @param newdata Optional matrix/data frame of target covariates.
#' @param X0 Optional numeric matrix of targets (rows).
#' @param x0 Optional numeric vector for a single target.
#' @param se.fit Logical; if `TRUE`, also return jackknife standard errors.
#' @param interval One of `"none"` or `"confidence"`.
#' @param level Confidence level for interval output.
#' @param ... Not used.
#'
#' @return If `interval = "none"` and `se.fit = FALSE`, a numeric vector of
#'   predictions. Otherwise a matrix/list analogous to `predict.lm()`.
#' @examples
#' set.seed(3)
#' n <- 40
#' x_mat <- matrix(rnorm(n * 2), n, 2)
#' fit <- cate_dnn(X = x_mat, x0 = c(0, 0), pseudo_outcome = rnorm(n), se = TRUE)
#' predict(fit, x0 = c(0.2, -0.1))
#' @export
predict.cate_result <- function(
    object,
    newdata = NULL,
    X0 = NULL,
    x0 = NULL,
    se.fit = FALSE,
    interval = c("none", "confidence"),
    level = 0.95,
    ...) {
  interval <- match.arg(interval)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    .abort_bad_arg("level", "must be a numeric scalar in (0,1)")
  }
  if (!is.logical(se.fit) || length(se.fit) != 1L || is.na(se.fit)) {
    .abort_bad_arg("se.fit", "must be TRUE or FALSE")
  }

  model <- .cate_validate_prediction_model(object$.model)

  x0_mat <- .dnn_resolve_prediction_targets(model, newdata = newdata, X0 = X0, x0 = x0)

  need_se <- se.fit || interval == "confidence"
  wt_return <- if (need_se) c("W", "ranks") else "W"
  wt <- weights_dnn(
    X = model$x_mat,
    X0 = x0_mat,
    distance = model$distance,
    method = model$method,
    k = model$k,
    s = model$s,
    normalize = model$normalize,
    ctx = model$ctx,
    return = wt_return
  )
  fit <- as.numeric(.compensated_crossprod(model$psi, wt$W))

  if (!need_se) {
    return(fit)
  }

  if (as.integer(model$s) >= nrow(model$x_mat)) {
    .abort_bad_arg(
      "model$s",
      "must satisfy s < n when jackknife-based prediction intervals are requested",
      expected = sprintf("s < n=%d", nrow(model$x_mat)),
      actual = sprintf("s=%d", as.integer(model$s)),
      quoted = FALSE
    )
  }
  se_hat <- .cate_dnn_jackknife_se_loo(
    x_mat = model$x_mat,
    x0_mat = x0_mat,
    psi = model$psi,
    distance = model$distance,
    method = model$method,
    k = as.integer(model$k),
    s = as.integer(model$s),
    normalize = isTRUE(model$normalize),
    ctx = model$ctx,
    jackknife = "fast",
    precomputed_ranks = wt$ranks
  )

  if (interval == "none") {
    if (!se.fit) {
      return(fit)
    }
    return(list(fit = fit, se.fit = se_hat))
  }

  crit <- stats::qnorm(0.5 + level / 2)
  pred_mat <- cbind(
    fit = fit,
    lwr = fit - crit * se_hat,
    upr = fit + crit * se_hat
  )
  if (!se.fit) {
    return(pred_mat)
  }
  list(fit = pred_mat, se.fit = se_hat)
}

#' Extract coefficients from a `cate_result`
#'
#' @param object A `cate_result` object.
#' @param ... Not used.
#' @return Named numeric vector of point estimates (`object$tau`), one entry
#'   per target point.
#' @export
coef.cate_result <- function(object, ...) {
  stats::setNames(as.numeric(object$tau), paste0("target_", seq_along(object$tau)))
}
