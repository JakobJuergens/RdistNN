# Shared helpers for checkout benchmarks; never installed as package code.
rdistnn_bench_require <- function() {
  if (!requireNamespace("bench", quietly = TRUE)) {
    stop("Please install the 'bench' package (install.packages('bench')).")
  }
  invisible(TRUE)
}

rdistnn_bench_make_idx <- function(n, n_subsets, subset_size) {
  draws <- vapply(seq_len(n_subsets), function(i)
    sample.int(n, subset_size, replace = FALSE), integer(subset_size))
  t(matrix(draws, nrow = subset_size, ncol = n_subsets))
}

rdistnn_bench_parity <- function(a, b, tol, label) {
  if (!isTRUE(all.equal(a, b, tolerance = tol))) {
    stop(
      sprintf("%s: compiled vs R mismatch beyond tolerance (tol=%.2e).", label, tol)
    )
  }
}

rdistnn_bench_print <- function(b) {
  print(b[, c("expression", "min", "median", "itr/sec", "mem_alloc", "gc/sec")])
}

# Compare the current score APIs and native reduction against explicit row sums.
rdistnn_bench_score <- function(score, data, idx, signal, reps, include_cpp_direct) {
  ctx <- score_setup(score, data)
  expected <- rowSums(matrix(signal[idx], nrow(idx), ncol(idx)))
  batch <- function() score_eval_batch(score, idx, data, ctx)
  scalar <- function() vapply(seq_len(nrow(idx)), function(i)
    score_eval(score, idx[i, ], data, ctx), numeric(1L))
  rdistnn_bench_parity(expected, batch(), 1e-12, "score batch")
  rdistnn_bench_parity(expected, scalar(), 1e-12, "score scalar")
  if (include_cpp_direct) {
    direct <- function() rustat_subset_row_sums_g(idx, signal)
    rdistnn_bench_parity(expected, direct(), 1e-12, "native subset reduction")
    result <- bench::mark(R = rowSums(matrix(signal[idx], nrow(idx), ncol(idx))),
      batch = batch(), scalar = scalar(), CPP = direct(), iterations = reps, check = TRUE)
  } else {
    result <- bench::mark(R = rowSums(matrix(signal[idx], nrow(idx), ncol(idx))),
      batch = batch(), scalar = scalar(), iterations = reps, check = TRUE)
  }
  rdistnn_bench_print(result)
  invisible(result)
}
