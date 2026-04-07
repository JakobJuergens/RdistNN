# Methods for cate_result objects.

.cate_resolve_targets <- function(model, newdata = NULL, X0 = NULL, x0 = NULL) {
  if (!is.null(newdata) && (!is.null(X0) || !is.null(x0))) {
    stop("Provide either `newdata` or `X0`/`x0`, not both.", call. = FALSE)
  }

  if (!is.null(newdata)) {
    if (is.data.frame(newdata) && !is.null(model$x_cols)) {
      x0_mat <- as.matrix(newdata[, model$x_cols, drop = FALSE])
    } else {
      x0_mat <- as.matrix(newdata)
    }
    return(.validate_ncol_relation(model$x_mat, x0_mat, y_arg = "newdata", x_arg = "X", relation = "match"))
  }

  if (!is.null(X0) && !is.null(x0)) {
    stop("Provide either `X0` or `x0`, not both.", call. = FALSE)
  }
  if (is.null(X0)) {
    if (is.null(x0)) {
      stop("Provide `newdata`, `X0`, or a single `x0`.", call. = FALSE)
    }
    x0_vec <- .validate_target_vector(x0, ncol(model$x_mat), x0_arg = "x0", x_arg = "X")
    return(matrix(x0_vec, nrow = 1L))
  }
  .validate_ncol_relation(model$x_mat, X0, y_arg = "X0", x_arg = "X", relation = "equal")
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
    stop("Standard errors are not available. Refit with `se = TRUE`.", call. = FALSE)
  }

  m <- length(object$tau)
  idx <- if (is.null(parm)) {
    seq_len(m)
  } else {
    as.integer(parm)
  }
  if (length(idx) == 0L || anyNA(idx) || any(idx < 1L | idx > m)) {
    .abort_bad_arg("parm", "must contain valid target indices in 1..length(tau)")
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

  model <- object$.model
  if (!is.list(model) || is.null(model$x_mat) || is.null(model$psi)) {
    stop(
      "This object does not contain fitted model state for prediction. ",
      "Refit with the current cate_dnn() implementation.",
      call. = FALSE
    )
  }

  x0_mat <- .cate_resolve_targets(model, newdata = newdata, X0 = X0, x0 = x0)

  need_se <- isTRUE(se.fit) || interval == "confidence"
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
  fit <- as.numeric(crossprod(model$psi, wt$W))

  if (!need_se) {
    return(fit)
  }

  if (as.integer(model$s) >= nrow(model$x_mat)) {
    stop(sprintf(
      "Delete-1 jackknife requires s < n (got s=%d, n=%d). Use `se.fit = FALSE` or `interval = 'none'`.",
      as.integer(model$s), nrow(model$x_mat)
    ), call. = FALSE)
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
    if (!isTRUE(se.fit)) {
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
  if (!isTRUE(se.fit)) {
    return(pred_mat)
  }
  list(fit = pred_mat, se.fit = se_hat)
}
