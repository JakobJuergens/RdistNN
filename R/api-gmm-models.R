# GMM model objects for local DNN-GMM estimators.

.gmm_validate_nonempty_chr <- function(x, arg) {
  if (!is.character(x) || length(x) < 1L || anyNA(x) || any(!nzchar(x))) {
    .abort_bad_arg(arg, "must be a non-empty character vector")
  }
  as.character(x)
}

#' Create a local GMM model object
#'
#' Defines a vector-moment model for use with [gmm_dnn()].  A model supplies
#' observation-level moments and, optionally, observation-level Jacobians and
#' target-specific starting values.
#'
#' The model functions are called with named arguments:
#' `setup(sample, n, ...)`, `moment_obs(sample, theta, ctx)`,
#' `jacobian_obs(sample, theta, ctx)`, and `start(sample, ctx, w, x0)`.
#' The `setup()` return value is passed back as `ctx`.
#'
#' @param name Non-empty model name.
#' @param setup Function that validates and prepares model data.
#' @param moment_obs Function returning an `n x q` numeric moment matrix.
#' @param jacobian_obs Optional function returning observation-level Jacobians.
#'   For scalar `theta`, an `n x q` matrix is accepted; otherwise use an
#'   `n x q x p` array.
#' @param start Optional function `(sample, ctx, w, x0)` returning a finite
#'   numeric starting vector, evaluated independently for each target.
#'   An explicit `start` passed to [gmm_dnn()] takes precedence.
#' @param coef_names Non-empty character vector naming the parameter entries.
#' @param moment_names Optional character vector naming the moment components.
#'
#' @return An object of class `rdistnn_gmm_model`.
#' @examples
#' mod <- new_gmm_model(
#'   name = "mean_zero",
#'   setup = function(sample, n, ...) list(g = sample$g),
#'   moment_obs = function(sample, theta, ctx) matrix(ctx$g - theta[1], ncol = 1),
#'   jacobian_obs = function(sample, theta, ctx) matrix(-1, length(ctx$g), 1),
#'   start = function(sample, ctx, w, x0) sum(w * ctx$g),
#'   coef_names = "theta",
#'   moment_names = "mean"
#' )
#' is_gmm_model(mod)
#' @export
new_gmm_model <- function(name,
                          setup,
                          moment_obs,
                          jacobian_obs = NULL,
                          start = NULL,
                          coef_names,
                          moment_names = NULL) {
  name <- .gmm_validate_nonempty_chr(name, "name")
  if (length(name) != 1L) {
    .abort_bad_arg("name", "must be a character scalar")
  }
  if (!is.function(setup)) {
    .abort_bad_arg("setup", "must be a function")
  }
  if (!is.function(moment_obs)) {
    .abort_bad_arg("moment_obs", "must be a function")
  }
  if (!is.null(jacobian_obs) && !is.function(jacobian_obs)) {
    .abort_bad_arg("jacobian_obs", "must be NULL or a function")
  }
  if (!is.null(start) && !is.function(start)) {
    .abort_bad_arg("start", "must be NULL or a function")
  }

  coef_names <- .gmm_validate_nonempty_chr(coef_names, "coef_names")
  if (!is.null(moment_names)) {
    moment_names <- .gmm_validate_nonempty_chr(moment_names, "moment_names")
  }

  out <- list(
    name = name,
    setup = setup,
    moment_obs = moment_obs,
    jacobian_obs = jacobian_obs,
    start = start,
    coef_names = coef_names,
    moment_names = moment_names,
    p = length(coef_names)
  )
  class(out) <- "rdistnn_gmm_model"
  out
}

#' Test whether an object is a GMM model
#'
#' @param x Object to test.
#' @return `TRUE` when `x` inherits from class `rdistnn_gmm_model`.
#' @export
is_gmm_model <- function(x) inherits(x, "rdistnn_gmm_model")

.gmm_validate_model_object <- function(model, arg = "model") {
  if (!is_gmm_model(model)) {
    .abort_bad_arg(arg, "must be a `rdistnn_gmm_model` object")
  }
  if (!is.list(model)) {
    .abort_bad_arg(arg, "must be a list-based `rdistnn_gmm_model` object")
  }

  name <- .gmm_validate_nonempty_chr(model$name, "model$name")
  if (length(name) != 1L) {
    .abort_bad_arg("model$name", "must be a character scalar")
  }
  if (!is.function(model$setup)) {
    .abort_bad_arg("model$setup", "must be a function")
  }
  if (!is.function(model$moment_obs)) {
    .abort_bad_arg("model$moment_obs", "must be a function")
  }
  if (!is.null(model$jacobian_obs) && !is.function(model$jacobian_obs)) {
    .abort_bad_arg("model$jacobian_obs", "must be NULL or a function")
  }
  if (!is.null(model$start) && !is.function(model$start)) {
    .abort_bad_arg("model$start", "must be NULL or a function")
  }

  coef_names <- .gmm_validate_nonempty_chr(model$coef_names, "model$coef_names")
  moment_names <- NULL
  if (!is.null(model$moment_names)) {
    moment_names <- .gmm_validate_nonempty_chr(model$moment_names, "model$moment_names")
  }
  if (!is.null(model$p)) {
    p <- .validate_scalar_int(model$p, arg = "model$p", min = 1L)
    if (p != length(coef_names)) {
      .abort_bad_arg(
        "model$p",
        "must match length(model$coef_names)",
        expected = sprintf("length(model$coef_names)=%d", length(coef_names)),
        actual = sprintf("model$p=%d", p)
      )
    }
  }

  model$name <- name
  model$coef_names <- coef_names
  model$moment_names <- moment_names
  model$p <- length(coef_names)
  class(model) <- union("rdistnn_gmm_model", class(model))
  model
}

.gmm_dr_iv_validate_sample <- function(sample, n) {
  if (!is.data.frame(sample)) {
    .abort_bad_arg(
      "model_data",
      "must be a data.frame for `gmm_model_dr_iv()`",
      expected = "`data` in the data-frame interface, or `model_data` in the matrix interface"
    )
  }
  if (nrow(sample) != n) {
    .abort_bad_arg(
      "nrow(model_data)",
      "must equal nrow(X)",
      expected = sprintf("nrow(X)=%d", n),
      actual = sprintf("nrow(model_data)=%d", nrow(sample)),
      quoted = FALSE
    )
  }
  invisible(NULL)
}

#' Doubly robust IV GMM model
#'
#' Creates the scalar local IV model used by the GMM IV simulation.  The
#' observation-level moment is
#' \deqn{
#'   (Z_j - \pi_j(X))\{Y - \mu_Y(X) - (D - \mu_D(X))\theta\},
#' }
#' evaluated for each instrument `j`.
#'
#' No nuisance models are fit internally.  All nuisance values must be supplied
#' as columns in `data` or `model_data`.
#'
#' @param y_col,d_col Column names or indices for outcome and endogenous
#'   treatment.
#' @param z_cols Column names or indices for instruments.
#' @param pi_cols Column names or indices for \eqn{E[Z_j \mid X]}.
#' @param muY_col Column name or index for \eqn{E[Y \mid X]}.
#' @param muD_col Column name or index for \eqn{E[D \mid X]}.
#'
#' @return A `rdistnn_gmm_model` object for use with [gmm_dnn()].
#' @examples
#' dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 80, kx = 2, q = 2, seed = 1)
#' model <- gmm_model_dr_iv(
#'   z_cols = c("Z1", "Z2"),
#'   pi_cols = c("piZ1_oracle", "piZ2_oracle"),
#'   muY_col = "muY_oracle",
#'   muD_col = "muD_oracle"
#' )
#' is_gmm_model(model)
#' @export
gmm_model_dr_iv <- function(y_col = "Y",
                            d_col = "D",
                            z_cols,
                            pi_cols,
                            muY_col,
                            muD_col) {
  if (missing(z_cols) || missing(pi_cols) || missing(muY_col) || missing(muD_col)) {
    .abort_bad_arg(
      "z_cols/pi_cols/muY_col/muD_col",
      "must all be provided",
      quoted = FALSE
    )
  }
  force(y_col)
  force(d_col)
  force(pi_cols)
  force(muY_col)
  force(muD_col)
  moment_names <- paste0("dr_iv_", seq_along(z_cols))

  setup <- function(sample, n, ...) {
    .gmm_dr_iv_validate_sample(sample, n)

    y_vec <- .dnn_col_vec(y_col, sample, n, arg = "y_col", ref = "model_data")
    d_vec <- .dnn_col_vec(d_col, sample, n, arg = "d_col", ref = "model_data")
    z_mat <- .dnn_col_mat(z_cols, sample, arg = "z_cols", out_arg = "Z")
    pi_mat <- .dnn_col_mat(pi_cols, sample, arg = "pi_cols", out_arg = "pi")
    muY <- .dnn_col_vec(muY_col, sample, n, arg = "muY_col", ref = "model_data")
    muD <- .dnn_col_vec(muD_col, sample, n, arg = "muD_col", ref = "model_data")

    if (nrow(z_mat) != n || nrow(pi_mat) != n) {
      .abort_bad_arg(
        "z_cols/pi_cols",
        "must resolve to matrices with nrow equal to nrow(X)",
        expected = sprintf("nrow(X)=%d", n),
        actual = sprintf("nrow(Z)=%d, nrow(pi)=%d", nrow(z_mat), nrow(pi_mat)),
        quoted = FALSE
      )
    }
    if (ncol(z_mat) != ncol(pi_mat)) {
      .abort_bad_arg(
        "pi_cols",
        "must contain the same number of columns as `z_cols`",
        expected = sprintf("length(z_cols)=%d", ncol(z_mat)),
        actual = sprintf("length(pi_cols)=%d", ncol(pi_mat))
      )
    }

    .abort_if_na_or_nonfinite(y_vec, arg = "Y")
    .abort_if_na_or_nonfinite(d_vec, arg = "D")
    .abort_if_na_or_nonfinite(z_mat, arg = "Z")
    .abort_if_na_or_nonfinite(pi_mat, arg = "pi")
    .abort_if_na_or_nonfinite(muY, arg = "muY")
    .abort_if_na_or_nonfinite(muD, arg = "muD")

    list(
      y = y_vec,
      d = d_vec,
      z_resid = z_mat - pi_mat,
      muY = muY,
      muD = muD,
      q = ncol(z_mat)
    )
  }

  moment_obs <- function(sample, theta, ctx) {
    theta <- as.numeric(theta)
    if (length(theta) != 1L || !is.finite(theta)) {
      .abort_bad_arg("theta", "must be a finite numeric scalar")
    }
    resid <- ctx$y - ctx$muY - (ctx$d - ctx$muD) * theta
    ctx$z_resid * resid
  }

  jacobian_obs <- function(sample, theta, ctx) {
    theta <- as.numeric(theta)
    if (length(theta) != 1L || !is.finite(theta)) {
      .abort_bad_arg("theta", "must be a finite numeric scalar")
    }
    ctx$z_resid * (-(ctx$d - ctx$muD))
  }

  start <- function(sample, ctx, w, x0) {
    a_obs <- ctx$z_resid * (ctx$y - ctx$muY)
    b_obs <- ctx$z_resid * (ctx$d - ctx$muD)
    a_bar <- as.numeric(crossprod(a_obs, w))
    b_bar <- as.numeric(crossprod(b_obs, w))
    denom <- sum(b_bar^2)
    if (!is.finite(denom) || denom <= .Machine$double.eps) {
      return(0)
    }
    sum(b_bar * a_bar) / denom
  }

  new_gmm_model(
    name = "dr_iv",
    setup = setup,
    moment_obs = moment_obs,
    jacobian_obs = jacobian_obs,
    start = start,
    coef_names = "theta",
    moment_names = moment_names
  )
}
