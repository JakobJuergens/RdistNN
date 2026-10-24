#' @keywords internal
rdistnn_bench_require <- function() {
  if (!requireNamespace("bench", quietly = TRUE)) {
    .dnn_abort("dependency_missing", "Please install the 'bench' package (install.packages('bench')).")
  }
  invisible(TRUE)
}

#' @keywords internal
rdistnn_bench_make_idx <- function(n, n_subsets, subset_size) {
  replicate(n_subsets, sample.int(n, subset_size, replace = FALSE)) |> t()
}

#' @keywords internal
rdistnn_bench_call_eval <- function(score_obj, idx) {
  if (is.function(score_obj)) {
    score_obj(idx)
  } else if (!is.null(score_obj$eval_batch) && is.function(score_obj$eval_batch)) {
    score_obj$eval_batch(idx)
  } else if (!is.null(score_obj$eval) && is.function(score_obj$eval)) {
    apply(idx, 1L, function(r) score_obj$eval(r))
  } else {
    .dnn_abort("argument_invalid",
      "Don't know how to call this score object. Expose $eval_batch() or return a function.",
      argument = "score_obj"
    )
  }
}

#' @keywords internal
rdistnn_bench_maybe_setup <- function(score_obj, idx) {
  setup <- score_obj$setup
  if (!is.null(setup) && is.function(setup)) setup(idx)
  invisible(NULL)
}

#' @keywords internal
rdistnn_bench_parity <- function(a, b, tol, label) {
  if (!isTRUE(all.equal(a, b, tolerance = tol))) {
    .dnn_abort("internal_invariant",
      sprintf("%s: compiled vs R mismatch beyond tolerance (tol=%.2e).", label, tol),
      hint = "Report a reproducible example of this internal inconsistency."
    )
  }
}

#' @keywords internal
rdistnn_bench_print <- function(b) {
  print(b[, c("expression", "min", "median", "itr/sec", "mem_alloc", "gc/sec")])
}
