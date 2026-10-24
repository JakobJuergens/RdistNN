# inst/bench/bench_cate.R
# Benchmark for the CATE rho *pure-score* with optional C++ helpers.

# --- Minimal helper bootstrapping --------------------------------------------
if (!exists("rdistnn_bench_require", inherits = TRUE)) {
  cand <- c(system.file("bench", "utils-bench.R", package = "RdistNN"),
            file.path("inst", "bench", "utils-bench.R"))
  cand <- cand[file.exists(cand)]
  for (f in cand) {
    try(suppressWarnings(source(f, local = TRUE)), silent = TRUE)
    if (exists("rdistnn_bench_require", inherits = TRUE)) break
  }
  if (!exists("rdistnn_bench_require", inherits = TRUE)) {
    rdistnn_bench_require <- function() {
      if (!requireNamespace("bench", quietly = TRUE))
        stop("Please install the 'bench' package (install.packages('bench')).")
      invisible(TRUE)
    }
  }
  if (!exists("%||%", inherits = TRUE)) `%||%` <- function(a, b) if (is.null(a)) b else a
  if (!exists("rdistnn_bench_make_idx", inherits = TRUE)) {
    rdistnn_bench_make_idx <- function(n, n_subsets, subset_size) {
      t(replicate(n_subsets, sample.int(n, subset_size, replace = FALSE)))
    }
  }
  if (!exists("rdistnn_bench_parity", inherits = TRUE)) {
    rdistnn_bench_parity <- function(a, b, tol, label) {
      if (!isTRUE(all.equal(a, b, tolerance = tol)))
        stop(sprintf("%s: mismatch beyond tol=%g.", label, tol))
    }
  }
  if (!exists("rdistnn_bench_print", inherits = TRUE)) {
    rdistnn_bench_print <- function(b) {
      print(b[, c("expression", "min", "median", "itr/sec", "mem_alloc", "gc/sec")])
    }
  }
}

# --- Benchmark ---------------------------------------------------------------
rdistnn_bench_cate <- function(
    n = 2000L,
    pX = 5L,
    s = 32L,                 # pure-score order
    n_subsets = 512L,
    reps = 5L,
    seed = 1L,
    use_rho_col = FALSE
) {
  rdistnn_bench_require()
  set.seed(seed)

  # Data + pseudo-outcome g = rho
  X <- matrix(rnorm(n * pX), n, pX)
  colnames(X) <- paste0("X", seq_len(pX))
  rho <- rowMeans(X[, 1:min(3L, pX), drop = FALSE]) + rnorm(n, sd = 0.5)

  df <- as.data.frame(X)
  if (use_rho_col) df$rho <- rho

  x0     <- colMeans(X)   # checked only for dimension
  x_cols <- seq_len(pX)

  # Build pure-score (no internal weights)
  sc <- score_cate_rho(
    x0 = x0,
    x_cols = x_cols,
    s = as.integer(s),
    rho_col = if (use_rho_col) "rho" else NULL,
    rho     = if (!use_rho_col) rho else NULL
  )
  ctx <- sc$setup(df)

  if (is.null(ctx$g) || length(ctx$g) != n) stop("ctx$g not set by setup().")
  if (n < s) stop(sprintf("n (%d) must be >= s (%d).", n, s), call. = FALSE)

  idx <- rdistnn_bench_make_idx(n, n_subsets, s)

  # Score + baseline evaluators
  eval_score_batch <- function() as.numeric(sc$eval_batch(idx, df, ctx, NULL))
  eval_score_loop  <- function() apply(idx, 1L, function(r) sc$eval(r, df, ctx, NULL))
  eval_baseline    <- function() {
    gi <- ctx$g
    rowSums(matrix(gi[idx], nrow = nrow(idx)))
  }

  # Benchmark the same native subset reduction used by the production score.
  subset <- get0("rustat_subset_row_sums_g", envir = asNamespace("RdistNN"), inherits = FALSE)
  have_cpp <- is.function(subset)
  eval_cpp_subset <- function() subset(idx, ctx$g)

  # Parity checks
  sb <- eval_score_batch()
  rdistnn_bench_parity(sb, eval_baseline(), 1e-12, "score vs baseline")
  rdistnn_bench_parity(sb, as.numeric(eval_score_loop()), 1e-12, "batch vs loop")
  if (have_cpp) {
    rdistnn_bench_parity(sb, eval_cpp_subset(), 1e-12, "score vs native subset")
  }

  # Benchmark (explicit calls; no do.call tricks)
  if (have_cpp) {
    b <- bench::mark(
      CATE_rho_score_eval_batch = eval_score_batch(),
      CATE_rho_baseline_rowSums = eval_baseline(),
      CATE_rho_loop_eval        = eval_score_loop(),
      CATE_rho_native_subset    = eval_cpp_subset(),
      iterations = reps,
      check = TRUE
    )
  } else {
    b <- bench::mark(
      CATE_rho_score_eval_batch = eval_score_batch(),
      CATE_rho_baseline_rowSums = eval_baseline(),
      CATE_rho_loop_eval        = eval_score_loop(),
      iterations = reps,
      check = TRUE
    )
  }

  rdistnn_bench_print(b)
  invisible(b)
}
