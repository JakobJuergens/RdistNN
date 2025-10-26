# inst/bench/bench_cate.R
# Benchmark for the CATE score evaluator

rdistnn_bench_cate <- function(
    n = 2000L,
    pX = 5L,
    n_subsets = 512L,
    reps = 5L,
    tol = 1e-10,
    seed = 1L,
    # score_cate_rho knobs:
    localizer = c("dnnk","dnn1","kernel"),
    k = 5L,
    h = NULL,
    use_rho_col = FALSE  # TRUE to exercise the rho_col branch
) {
  rdistnn_bench_require()
  set.seed(seed)
  localizer <- match.arg(localizer)

  # --- Toy data with a pseudo-outcome rho ---
  X <- matrix(rnorm(n * pX), n, pX)
  colnames(X) <- paste0("X", seq_len(pX))
  # any numeric vector works for the precomputed pseudo-outcome
  rho <- rowMeans(X[, 1:min(3L, pX), drop = FALSE]) + rnorm(n, sd = 0.5)

  df <- as.data.frame(X)
  if (use_rho_col) df$rho <- rho

  x0 <- colMeans(X)
  x_cols <- seq_len(pX)

  # --- Build score via current API (score_cate_rho) ---
  sc <- score_cate_rho(
    x0 = x0,
    x_cols = x_cols,
    rho_col = if (use_rho_col) "rho" else NULL,
    rho = if (!use_rho_col) rho else NULL,
    localizer = localizer,
    k = as.integer(k),
    h = h
  )

  # setup() returns ctx; grab s from there for the subset order
  ctx_compiled <- sc$setup(df)
  s <- as.integer(ctx_compiled$s %||% 32L)

  # build subset matrix at the kernel order s
  idx <- rdistnn_bench_make_idx(n, n_subsets, s)

  # Force the R fallback by nulling the compiled pointer
  ctx_R <- ctx_compiled
  ctx_R$ptr_ctx <- NULL

  compiled_available <- !is.null(ctx_compiled$ptr_ctx) &&
    exists("rustat_cate_eval_batch", mode = "function")

  # One-shot evaluators
  eval_R <- function() {
    as.numeric(sc$eval_batch(idx, df, ctx_R, NULL))
  }
  eval_CPP <- function() {
    as.numeric(sc$eval_batch(idx, df, ctx_compiled, NULL))
  }

  # Parity (only if compiled path exists)
  if (compiled_available) {
    rdistnn_bench_parity(eval_R(), eval_CPP(), tol, "CATE rho score")
  } else {
    message("Compiled symbols not found; benchmarking R fallback only.")
  }

  # --- Benchmark ---
  if (compiled_available) {
    b <- bench::mark(
      CATE_rho_R     = eval_R(),
      CATE_rho_comp  = eval_CPP(),
      iterations = reps,
      check = TRUE
    )
  } else {
    b <- bench::mark(
      CATE_rho_R     = eval_R(),
      iterations = reps,
      check = TRUE
    )
  }

  rdistnn_bench_print(b)
  invisible(b)
}
