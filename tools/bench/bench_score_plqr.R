# inst/bench/bench_plqr.R
# Benchmark for the PLQR component *pure-score* (score_plqr_component)
# with optional C++ helpers if available.

# --- Minimal helper bootstrap (internal utils or inst/bench) ------------------
if (!exists("rdistnn_bench_require", inherits = TRUE)) {
  cand <- c(system.file("bench", "utils-bench.R", package = "RdistNN"),
            file.path("inst","bench","utils-bench.R"))
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
      print(b[, c("expression","min","median","itr/sec","mem_alloc","gc/sec")])
    }
  }
}

# --- Benchmark function -------------------------------------------------------
rdistnn_bench_plqr <- function(
    n = 2000L,
    pX = 5L,
    pD = 3L,
    s  = 32L,               # kernel order (explicit in pure-score API)
    tau = 0.5,
    j   = 1L,
    beta = NULL,            # defaults to rep(0, pD) if NULL
    n_subsets = 512L,
    reps = 5L,
    seed = 1L,
    pass_names = FALSE      # TRUE: pass names; FALSE: pass numeric indices
) {
  rdistnn_bench_require()
  set.seed(seed)

  # --- Oracle-friendly synthetic data ----------------------------------------
  # X
  X <- matrix(rnorm(n * pX), n, pX, dimnames = list(NULL, paste0("X", seq_len(pX))))

  # q_tau(X): smooth scalar of X
  qx <- drop(0.3 * X[, 1L]) + 0.15 * rowSums(X[, 1:min(3L, pX), drop = FALSE])

  # m_D(X): n x pD (simple dense map of X)
  B  <- matrix(0.2, nrow = pD, ncol = pX)
  mD <- X %*% t(B)
  colnames(mD) <- paste0("mD", seq_len(pD))

  # D | X: centered around m_D(X)
  D <- mD + matrix(rnorm(n * pD, sd = 0.5), n, pD,
                   dimnames = list(NULL, paste0("D", seq_len(pD))))

  # beta and Y so that P(Y <= D^T beta + q) = tau in expectation
  if (is.null(beta)) beta <- rep(0, pD) else beta <- as.numeric(beta)
  e <- rlogis(n) - qlogis(tau)                  # ensures P(e <= 0) = tau
  Y <- as.numeric(D %*% beta + qx + e)

  # Assemble df in a *known* order: [ Y | X | D | q | mD ]
  df <- data.frame(Y = Y, X, D, q = qx, mD, check.names = FALSE)

  # Numeric column positions (robust to name munging)
  i_Y  <- 1L
  i_X  <- seq.int(2L, 1L + pX)
  i_D  <- seq.int(2L + pX, 1L + pX + pD)
  i_q  <- 2L + pX + pD
  i_mD <- seq.int(i_q + 1L, i_q + pD)

  # Or names if requested
  x_cols  <- if (pass_names) names(df)[i_X]  else i_X
  D_cols  <- if (pass_names) names(df)[i_D]  else i_D
  q_col   <- if (pass_names) names(df)[i_q]  else i_q
  mD_cols <- if (pass_names) names(df)[i_mD] else i_mD

  # Target point x0 (dimension check only for the pure score)
  x0 <- colMeans(X)

  # --- Build pure-score -------------------------------------------------------
  sc <- score_plqr_component(
    x0 = x0, x_cols = x_cols, D_cols = D_cols, q_col = q_col, mD_cols = mD_cols,
    beta = beta, j = as.integer(j), tau = tau,
    s = as.integer(s), k = 1L
  )

  # setup() => ctx with g
  ctx <- sc$setup(df)
  if (is.null(ctx$g) || length(ctx$g) != n) stop("ctx$g not set correctly by setup().")
  if (n < s) stop(sprintf("n (%d) must be >= s (%d).", n, s), call. = FALSE)

  # Subset indices
  idx <- rdistnn_bench_make_idx(n, n_subsets, s)

  # --- R evaluators -----------------------------------------------------------
  eval_score_batch <- function() as.numeric(sc$eval_batch(idx, df, ctx, NULL))
  eval_score_loop  <- function() apply(idx, 1L, function(r) sc$eval(r, df, ctx, NULL))
  eval_baseline    <- function() {
    gi <- ctx$g
    rowSums(matrix(gi[idx], nrow = nrow(idx)))
  }

  # --- Optional C++ helpers (robust scalar detection; ctx rebuilt from df) ---
  c_syms <- c("rustat_plqr_ctx_create",
              "rustat_plqr_ctx_set_beta",
              "rustat_plqr_g_all",
              "rustat_plqr_g_indexed")

  in_search <- all(vapply(c_syms, exists, logical(1), mode = "function"))
  in_ns <- FALSE
  ns <- NULL
  if ("RdistNN" %in% loadedNamespaces()) {
    ns <- asNamespace("RdistNN")
    in_ns <- all(vapply(c_syms, exists, logical(1), envir = ns, inherits = FALSE))
  }
  have_cpp <- isTRUE(in_search || in_ns)

  eval_cpp_indexed <- NULL
  eval_cpp_all     <- NULL

  if (have_cpp) {
    # Resolve helpers
    res <- function(name) {
      if (in_ns) get(name, envir = ns, inherits = FALSE) else get(name, inherits = TRUE)
    }
    plqr_ctx_create <- res("rustat_plqr_ctx_create")
    plqr_g_all      <- res("rustat_plqr_g_all")
    plqr_g_indexed  <- res("rustat_plqr_g_indexed")

    # REBUILD inputs in EXACTLY the same order used by the score (handles pass_names=TRUE)
    D_cpp  <- as.matrix(df[, D_cols,  drop = FALSE])
    mD_cpp <- as.matrix(df[, mD_cols, drop = FALSE])
    qx_cpp <- df[[q_col]]  # works for both numeric index and name
    Y_cpp  <- df$Y

    cpp_ctx <- plqr_ctx_create(Y_cpp, D_cpp, qx_cpp, mD_cpp, beta, tau, as.integer(j))

    # Column-major flattening to match gi[idx] semantics
    eval_cpp_indexed <- function() {
      flat <- as.integer(idx)                # column-major
      gi   <- plqr_g_indexed(flat, cpp_ctx)  # length K*s
      rowSums(matrix(gi, nrow = nrow(idx)))  # default fill = column-major
    }
    eval_cpp_all <- function() {
      gi <- plqr_g_all(cpp_ctx)               # length n
      rowSums(matrix(gi[idx], nrow = nrow(idx)))
    }
  }


  # --- Parity checks (fail fast) ---------------------------------------------
  sb <- eval_score_batch()
  rdistnn_bench_parity(sb, eval_baseline(), 1e-12, "score vs baseline")
  rdistnn_bench_parity(sb, as.numeric(eval_score_loop()), 1e-12, "batch vs loop")
  if (have_cpp) {
    rdistnn_bench_parity(sb, eval_cpp_indexed(), 1e-12, "score vs cpp_indexed")
    rdistnn_bench_parity(sb, eval_cpp_all(),     1e-12, "score vs cpp_all")
  }

  # --- Benchmark --------------------------------------------------------------
  if (have_cpp) {
    b <- bench::mark(
      PLQR_score_eval_batch = eval_score_batch(),
      PLQR_baseline_rowSums = eval_baseline(),
      PLQR_loop_eval        = eval_score_loop(),
      PLQR_cpp_indexed      = eval_cpp_indexed(),
      PLQR_cpp_all          = eval_cpp_all(),
      iterations = reps,
      check = TRUE
    )
  } else {
    b <- bench::mark(
      PLQR_score_eval_batch = eval_score_batch(),
      PLQR_baseline_rowSums = eval_baseline(),
      PLQR_loop_eval        = eval_score_loop(),
      iterations = reps,
      check = TRUE
    )
  }

  rdistnn_bench_print(b)
  invisible(b)
}
