#' Conditional LATE from externally supplied nuisance predictions
#'
#' An augmented binary-IV estimator with rich adjustment and
#' reduced-dimensional reporting. Nuisance functions are supplied as aligned
#' predictions; the estimator fits no learners.
#'
#' @inheritParams affine_ratio_dnn
#' @inheritParams iv_augmented_signals
#' @param y_col,d_col,a_col,e_col Data-frame columns for outcome, treatment,
#'   instrument and instrument propensity.
#' @param m_y0_col,m_y1_col Data-frame outcome predictions by instrument arm.
#' @param m_d0_col,m_d1_col Data-frame treatment predictions by instrument arm.
#'   The zero-arm prediction may be omitted with one_sided = TRUE.
#'
#' @details
#' The documented sampling framework uses independent observations; cluster
#' and survey-design inference are not supplied. Exact distance ties use
#' lexicographic order of the original covariates, in their supplied column
#' order. Identical covariate rows retain input-order ties; permuting those
#' duplicates can change estimates. Ordering uses no outcomes or nuisances.
#'
#' Nuisance predictions use all pre-instrument adjustment variables, while
#' weights use the reporting coordinates X. The estimate is the ratio of the
#' adjusted local reduced form and first stage, using common weights.
#'
#' The accompanying
#' [paper](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z)
#' gives the IV identifying assumptions and conditions for pointwise inference.
#' A training label or a computational choice of k or weight method does not
#' establish those conditions. Inspect the local first stage alongside the ratio.
#'
#' Covariance keeps paired variation of both frozen signals. The residualized
#' IV-GMM workflow uses a different moment and can have a different target.
#' Reporting uses the common generics; overlap panels describe the instrument
#' and influence panels retain both augmented signals.
#'
#' @return A clate_result inheriting affine_ratio_result. Component roles are
#'   reduced_form and first_stage. Optional nuisance output is the complete
#'   iv_signals construction record, including raw/used propensity.
#' @examples
#' d <- data.frame(
#'   x = 1:6, Y = c(1, 2, 0, 3, 2, 1),
#'   D = c(0, 1, 0, 1, 0, 1), A = c(0, 1, 1, 1, 0, 1),
#'   e = rep(.5, 6), y0 = rep(1, 6), y1 = rep(2, 6), d1 = rep(.6, 6)
#' )
#' fit <- clate_dnn(
#'   data = d, x_cols = "x", y_col = "Y", d_col = "D", a_col = "A",
#'   e_col = "e", m_y0_col = "y0", m_y1_col = "y1", m_d1_col = "d1",
#'   one_sided = TRUE, x0 = 3, s = 2
#' )
#' fit$components
#' fit$theta
#' clate_dnn(
#'   X = matrix(d$x, ncol = 1), Y = d$Y, D = d$D, A = d$A,
#'   e = d$e, m_y0 = d$y0, m_y1 = d$y1, m_d1 = d$d1,
#'   one_sided = TRUE, x0 = 3, s = 2
#' )$theta
#' @export
clate_dnn <- function(
  data = NULL, x_cols = NULL,
  y_col = NULL, d_col = NULL, a_col = NULL, e_col = NULL,
  m_y0_col = NULL, m_y1_col = NULL,
  m_d0_col = NULL, m_d1_col = NULL,
  X = NULL, X0 = NULL, x0 = NULL,
  Y = NULL, D = NULL, A = NULL, e = NULL,
  m_y0 = NULL, m_y1 = NULL, m_d0 = NULL, m_d1 = NULL,
  clip = 1e-8, one_sided = FALSE,
  distance = "l2", method = c("exact", "asymp"), k = 1L,
  s = NULL, tuning = NULL, tuning_mode = c("default", "grid"),
  normalize = TRUE,
  se = TRUE, jackknife = c("fast", "loo", "none"),
  denominator_tol = sqrt(.Machine$double.eps),
  row_ids = NULL, training = "supplied", ctx = NULL,
  return = c("W")
) {
  one_sided <- .validate_scalar_logical(one_sided, "one_sided")
  method <- .affine_weight_method(method)
  tuning_mode <- match.arg(tuning_mode)
  jackknife <- match.arg(jackknife)
  valid <- c(
    "W", "ranks", "ord", "dx", "ess", "signals",
    "deletions", "pseudo_values", "neighbors", "tuning", "nuisance"
  )
  .warn_unknown_return(return, valid = valid)
  return <- intersect(return, valid)
  column_arguments <- list(
    y_col = y_col, d_col = d_col, a_col = a_col, e_col = e_col,
    m_y0_col = m_y0_col, m_y1_col = m_y1_col,
    m_d0_col = m_d0_col, m_d1_col = m_d1_col
  )
  vectors <- list(
    Y = Y, D = D, A = A, e = e, m_y0 = m_y0,
    m_y1 = m_y1, m_d0 = m_d0, m_d1 = m_d1
  )
  resolved_x_cols <- NULL
  if (!is.null(data)) {
    if (!is.null(X) || any(!vapply(vectors, is.null, logical(1L)))) {
      .abort_bad_arg("data/X", "cannot combine data with X or direct sample vectors", quoted = FALSE)
    }
    if (!is.data.frame(data)) .abort_bad_arg("data", "must be a data.frame")
    n <- .affine_sample_size(data, "data")
    resolved_x_cols <- .dnn_resolve_data_cols(x_cols, data, "x_cols")
    x_mat <- .dnn_col_mat(resolved_x_cols, data, "x_cols", out_arg = "X")
    for (nm in names(vectors)) {
      col_arg <- switch(nm,
        Y = "y_col",
        D = "d_col",
        A = "a_col",
        e = "e_col",
        m_y0 = "m_y0_col",
        m_y1 = "m_y1_col",
        m_d0 = "m_d0_col",
        m_d1 = "m_d1_col"
      )
      if (nm == "m_d0" && one_sided && is.null(m_d0_col)) next
      vectors[[nm]] <- .affine_column(column_arguments[[col_arg]], data, n, col_arg)
    }
  } else {
    if (!is.null(x_cols) || any(!vapply(column_arguments, is.null, logical(1L)))) {
      .abort_bad_arg("data", "must be provided when column selectors are used")
    }
    if (is.null(X)) .abort_bad_arg("X", "must be supplied for the matrix interface")
    x_mat <- .validate_matrix_nonempty(X, "X")
    n <- .affine_sample_size(x_mat, "X")
    vectors$Y <- .affine_vector(Y, n, "Y")
  }
  signals <- do.call(iv_augmented_signals, c(
    vectors, list(clip = clip, one_sided = one_sided, row_ids = row_ids, training = training)
  ))
  out <- affine_ratio_dnn(
    X = x_mat, X0 = X0, x0 = x0, numerator = signals$phi_y, denominator = signals$phi_d,
    distance = distance, method = method, k = k, s = s, tuning = tuning,
    tuning_mode = tuning_mode, normalize = normalize,
    se = se, jackknife = jackknife, denominator_tol = denominator_tol,
    row_ids = signals$row_ids, ctx = ctx, return = setdiff(return, "nuisance")
  )
  out$meta$component_roles <- c(numerator = "reduced_form", denominator = "first_stage")
  out$meta$training <- training
  out$meta$one_sided <- one_sided
  out$meta$propensity_clip <- clip
  out$.model$x_cols <- resolved_x_cols
  out$.model$iv <- list(
    A = as.integer(vectors$A), D = as.integer(vectors$D),
    propensity = signals$propensity, one_sided = one_sided, training = training
  )
  if ("nuisance" %in% return) out$nuisance <- signals
  class(out) <- c("clate_result", "affine_ratio_result")
  out
}
