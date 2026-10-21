# Retained-state validation, provenance and CATE deletion adapter.

.cate_result_meta <- function(meta) {
  if (!is.list(meta) || !is.character(meta$signal_source) ||
      length(meta$signal_source) != 1L || is.na(meta$signal_source) ||
      !meta$signal_source %in% c("aipw", "supplied") ||
      !is.character(meta$created_with_version) || length(meta$created_with_version) != 1L ||
      is.na(meta$created_with_version) || !nzchar(meta$created_with_version)) {
    .abort_bad_state("This object lacks current CATE provenance metadata. Refit with cate_dnn().")
  }
  meta
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

# Internal: paper-consistent delete-1 jackknife for cate_dnn().
# Implements the delete-1 recomputation in Eq. (3.30)-(3.32) of
# Local_Z_Est_JJuergens.pdf for the scalar (p_alpha = 1) CATE target.
.cate_dnn_jackknife_se_loo <- function(
  x_mat, x0_mat, psi, distance, method, k, s, normalize, ctx,
  jackknife = c("fast", "loo"),
  precomputed_ranks = NULL
) {
  jackknife <- match.arg(jackknife)
  as.numeric(
    .dnn_jackknife_se_loo_matrix(
      x_mat = x_mat,
      x0_mat = x0_mat,
      signal = psi,
      distance = distance,
      method = method,
      k = k,
      s = s,
      normalize = normalize,
      ctx = ctx,
      signal_arg = "pseudo_outcome",
      jackknife = jackknife,
      precomputed_ranks = precomputed_ranks
    )[1L, ]
  )
}
