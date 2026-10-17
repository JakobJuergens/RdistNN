#' Construct augmented reduced-form and first-stage signals
#'
#' Constructs two aligned signals from externally supplied nuisance predictions.
#' This experimental interface fits no nuisance models.
#'
#' @param Y Finite numeric outcome vector.
#' @param D Numeric binary treatment vector (0/1).
#' @param A Numeric binary instrument vector (0/1).
#' @param e Instrument propensity predictions in \eqn{[0,1]}, before clipping.
#' @param m_y0,m_y1 Outcome predictions conditional on instrument arm 0 and 1.
#' @param m_d0,m_d1 Treatment predictions conditional on instrument arm 0 and 1,
#'   in \eqn{[0,1]}. m_d0 may be omitted only with one_sided = TRUE.
#' @param clip Propensity clipping threshold in (0,0.5).
#' @param one_sided Explicit structural restriction that treatment is impossible
#'   under instrument 0. Requires D <= A and zero m_d0 if supplied.
#' @param row_ids Optional unique, nonmissing character or numeric row identities.
#'   These record alignment; a single ID vector cannot detect independently
#'   misordered prediction arrays.
#' @param training Descriptive declaration: supplied, oracle,
#'   full_sample, oob or cross_fit. It is not verified from predictions.
#'
#' @details For either response \eqn{R}, the signal is
#' \deqn{\phi_R=m_{R1}-m_{R0}+A(R-m_{R1})/e
#'       -(1-A)(R-m_{R0})/(1-e).}
#' Both responses use the same clipped propensity. Arm subscripts refer to
#' instrument assignment, not treatment received. Treatment predictions are
#' validated, not silently clipped. Outcomes and outcome predictions are not
#' capped. Finite unbounded-outcome computations do not certify the bounded
#' learner conditions of the affine inference theorem.
#'
#' Individual first-stage signals can be negative even under one-sided
#' compliance. They remain unchanged. The CLATE equation averages both signals
#' before dividing; do not average observation-level signal ratios.
#'
#' @return An iv_signals list with phi_y, phi_d, row_ids,
#'   propensity (raw, used, clipped), nuisance
#'   (m_y0, m_y1, m_d0, m_d1) and meta.
#' @examples
#' z <- iv_augmented_signals(
#'   Y = c(1, 2, 0, 3), D = c(0, 1, 0, 1), A = c(0, 1, 1, 1),
#'   e = rep(.5, 4), m_y0 = rep(1, 4), m_y1 = rep(2, 4),
#'   m_d1 = rep(.6, 4), one_sided = TRUE
#' )
#' z$phi_y
#' z$phi_d
#' @export
iv_augmented_signals <- function(
  Y, D, A, e, m_y0, m_y1, m_d0 = NULL, m_d1,
  clip = 1e-8, one_sided = FALSE, row_ids = NULL,
  training = "supplied"
) {
  required <- c("Y", "D", "A", "e", "m_y0", "m_y1", "m_d1")
  absent <- setdiff(required, names(as.list(match.call())[-1L]))
  if (length(absent)) .abort_bad_arg(absent[1L], "must be provided")
  one_sided <- .validate_scalar_logical(one_sided, "one_sided")
  training_values <- c("supplied", "oracle", "full_sample", "oob", "cross_fit")
  if (!is.character(training) || length(training) != 1L ||
      is.na(training) || !training %in% training_values) {
    .abort_bad_arg("training", paste("must be one of", paste(training_values, collapse = ", ")))
  }
  y_vec <- .affine_vector(Y, arg = "Y")
  n <- length(y_vec)
  ids <- .affine_row_ids(row_ids, n)
  d_vec <- .validate_binary01(.affine_vector(D, n, "D"), "D")
  a_vec <- .validate_binary01(.affine_vector(A, n, "A"), "A")
  e_raw <- .affine_vector(e, n, "e")
  e_used <- .cate_validate_and_clip_propensity(e_raw, n, clip = clip, arg = "e")
  my0 <- .affine_vector(m_y0, n, "m_y0")
  my1 <- .affine_vector(m_y1, n, "m_y1")
  md1 <- .affine_vector(m_d1, n, "m_d1")
  md0 <- if (one_sided && is.null(m_d0)) numeric(n) else .affine_vector(m_d0, n, "m_d0")
  for (nm in c("m_d0", "m_d1")) {
    prediction <- if (nm == "m_d0") md0 else md1
    if (any(prediction < 0 | prediction > 1)) {
      .abort_bad_arg(nm, "must lie in [0, 1]")
    }
  }
  if (one_sided) {
    if (any(d_vec > a_vec)) {
      first <- which(d_vec > a_vec)[1L]
      .abort_bad_arg("one_sided", sprintf("requires D <= A; violation at row_id '%s'", ids[first]))
    }
    if (any(md0 != 0)) .abort_bad_arg("m_d0", "must be zero when one_sided = TRUE")
  }
  phi_y <- my1 - my0 + a_vec * (y_vec - my1) / e_used -
    (1 - a_vec) * (y_vec - my0) / (1 - e_used)
  phi_d <- md1 - md0 + a_vec * (d_vec - md1) / e_used -
    (1 - a_vec) * (d_vec - md0) / (1 - e_used)
  if (any(!is.finite(phi_y))) {
    .abort_bad_arg("Y/m_y0/m_y1/e", "produce non-finite augmented outcome signals; rescale the outcome inputs")
  }
  if (any(!is.finite(phi_d))) {
    .abort_bad_arg("D/m_d0/m_d1/e", "produce non-finite augmented first-stage signals")
  }
  out <- list(
    phi_y = as.numeric(phi_y), phi_d = as.numeric(phi_d), row_ids = ids,
    propensity = list(raw = e_raw, used = e_used, clipped = e_raw != e_used),
    nuisance = list(m_y0 = my0, m_y1 = my1, m_d0 = md0, m_d1 = md1),
    meta = list(n = n, clip = as.numeric(clip), one_sided = one_sided, training = training)
  )
  class(out) <- "iv_signals"
  out
}
