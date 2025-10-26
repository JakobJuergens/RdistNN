# inst/bench/bench_plqr.R
# Benchmark for the PLQR score evaluator

rdistnn_bench_plqr <- function(
    n = 2000L,
    pX = 5L,
    pD = 3L,
    tau = 0.5,
    j = 1L,
    reps = 5L,
    n_subsets = 512L,
    tol = 1e-10,
    seed = 1L,
    beta = NULL,                      # defaults to rep(0, pD) if NULL
    localizer = c("dnnk","dnn1","kernel"),
    k = 5L,
    h = NULL,
    pass_names = FALSE                # set TRUE if you prefer passing names
) {
  rdistnn_bench_require()
  set.seed(seed)
  localizer <- match.arg(localizer)

  # --- Simulate oracle-friendly data (finite & consistent) -------------------
  X <- matrix(rnorm(n * pX), n, pX, dimnames = list(NULL, paste0("X", seq_len(pX))))

  # q_tau(X): smooth scalar of X
  qx <- drop(0.3 * X[, 1L]) + 0.15 * rowSums(X[, 1:min(3L, pX), drop = FALSE])

  # m_D(X): n x pD
  B  <- matrix(0.2, nrow = pD, ncol = pX)
  mD <- X %*% t(B)
  colnames(mD) <- paste0("mD", seq_len(pD))

  # D | X: centered around m_D(X)
  D <- mD + matrix(rnorm(n * pD, sd = 0.5), n, pD, dimnames = list(NULL, paste0("D", seq_len(pD))))

  # beta and Y
  if (is.null(beta)) beta <- rep(0, pD) else beta <- as.numeric(beta)
  e <- rlogis(n) - qlogis(tau)         # ensures P(e <= 0) = tau
  Y <- as.numeric(D %*% beta + qx + e)

  # Assemble data.frame in a *known* column order:
  # [ Y | X(1..pX) | D(1..pD) | q | mD(1..pD) ]
  df <- data.frame(Y = Y, X, D, q = qx, mD, check.names = FALSE)

  # Compute numeric indices deterministically from construction order:
  i_Y   <- 1L
  i_X   <- seq.int(2L, 1L + pX)
  i_D   <- seq.int(2L + pX, 1L + pX + pD)
  i_q   <- 2L + pX + pD
  i_mD  <- seq.int(i_q + 1L, i_q + pD)

  # Optionally derive *names* as well (if you want to pass names)
  x_names  <- names(df)[i_X]
  D_names  <- names(df)[i_D]
  mD_names <- names(df)[i_mD]
  q_name   <- names(df)[i_q]

  # Target point x0
  x0 <- colMeans(X)

  # Choose whether to pass positions or names into the score factory
  x_cols  <- if (pass_names) x_names  else i_X
  D_cols  <- if (pass_names) D_names  else i_D
  mD_cols <- if (pass_names) mD_names else i_mD
  q_col   <- if (pass_names) q_name   else i_q

  # --- Build score via your wrapper -----------------------------------------
  sc <- score_plqr_component(
    x0 = x0, x_cols = x_cols, D_cols = D_cols,
    q_col = q_col, mD_cols = mD_cols,
    tau = tau, j = as.integer(j), beta = beta,
    localizer = localizer, k = as.integer(k), h = h
  )

  # setup() => ctx; use its s for the subset order
  ctx_compiled <- sc$setup(df)
  s <- as.integer(ctx_compiled$s %||% 32L)

  # Subset indices of size s
  idx <- rdistnn_bench_make_idx(n, n_subsets, s)

  # Force R fallback by nulling compiled pointer
  ctx_R <- ctx_compiled
  ctx_R$ptr_ctx <- NULL

  compiled_available <- !is.null(ctx_compiled$ptr_ctx) &&
    exists("rustat_plqr_eval_batch", mode = "function")

  eval_R <- function()  as.numeric(sc$eval_batch(idx, df, ctx_R,        NULL))
  eval_C <- function()  as.numeric(sc$eval_batch(idx, df, ctx_compiled, NULL))

  # Parity if compiled path exists
  if (compiled_available) {
    rdistnn_bench_parity(eval_R(), eval_C(), tol, "PLQR component")
  } else {
    message("Compiled symbols not found; benchmarking R fallback only.")
  }

  # --- Benchmark -------------------------------------------------------------
  if (compiled_available) {
    b <- bench::mark(
      PLQR_R    = eval_R(),
      PLQR_comp = eval_C(),
      iterations = reps,
      check = TRUE
    )
  } else {
    b <- bench::mark(
      PLQR_R    = eval_R(),
      iterations = reps,
      check = TRUE
    )
  }

  rdistnn_bench_print(b)
  invisible(b)
}
