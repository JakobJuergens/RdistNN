# tools/bench/bench_jackknife.R
# -------------------------------------------------------------------
# Benchmarks for jackknife_slow_rank() vs jackknife_fast_rank()
# Usage:
#   source("tools/bench/run_bench.R"); rdistnn_bench("jackknife", n = 2000)
#   Rscript tools/bench/run_bench.R jackknife n=2000 s=15 q=4 reps=10 phi_dep=TRUE

rdistnn_bench_jackknife <- function(n = 2000, p = 5, s = 10, q = 1,
                                    reps = 10, seed = 1,
                                    normalize = TRUE,
                                    phi_dep = FALSE,   # if TRUE, use phi(j,t) to force R fallback
                                    csv = NULL,        # optional path to write CSV
                                    verbose = TRUE) {
  if (!requireNamespace("bench", quietly = TRUE)) {
    stop("Install 'bench' to run: install.packages('bench')", call. = FALSE)
  }
  set.seed(seed)

  # --- deterministic top-(s+1) neighbor indices (squared L2), 1-based -----
  .sqdist <- function(A, B) {
    A <- as.matrix(A); B <- as.matrix(B)
    A2 <- rowSums(A^2); B2 <- rowSums(B^2)
    -2 * (A %*% t(B)) + matrix(A2, nrow(A), nrow(B)) +
      matrix(B2, nrow(A), nrow(B), byrow = TRUE)
  }
  .knn_ranks <- function(X, X0 = NULL, k = 5L, exclude_self = TRUE) {
    if (is.null(X0)) X0 <- X
    D <- .sqdist(X0, X)
    Tn <- nrow(D); n <- ncol(D)
    ranks <- matrix(NA_integer_, Tn, k)
    for (t in seq_len(Tn)) {
      d <- D[t, ]
      ord <- order(d, seq_len(n))         # stable tie-break by index
      if (exclude_self && identical(X0, X)) ord <- ord[ord != t]
      ranks[t, ] <- ord[seq_len(k)]
    }
    ranks
  }

  # --- data, ranks, weights, phi ------------------------------------------
  X <- matrix(rnorm(n * p), n, p)
  ranks <- .knn_ranks(X, k = s + 1L, exclude_self = TRUE)

  a <- runif(s); a <- a / mean(a)         # nontrivial L-stat coefficients

  # PHI matrix (n x q)
  PHI <- matrix(rnorm(n * q), n, q)

  if (!phi_dep) {
    phi_obj <- PHI                          # compiled fast path will be used if available
  } else {
    # target-dependent phi(j,t) (forces R fallback in jackknife_fast_rank)
    if (q == 1L) {
      # inject mild t-dependence to keep magnitudes sensible
      phi_obj <- function(j, t) PHI[j, 1] * (1 + 1e-3 * t)
    } else {
      phi_obj <- function(j, t) {
        v <- PHI[j, ]
        v[2] <- v[2] * t
        if (q >= 3L) v[3] <- 1.0
        v
      }
    }
  }

  compiled <- is.matrix(phi_obj) && is.numeric(phi_obj) &&
    is.function(.rdistnn_compiled_fun("rustat_jackknife_fast_rank"))
  if (verbose) {
    message(sprintf("jackknife benchmark — n=%d, p=%d, s=%d, q=%d, reps=%d, phi_dep=%s, compiled=%s",
                    n, p, s, q, reps, phi_dep, compiled))
  }

  # --- correctness warm-up (quick check) -----------------------------------
  ref <- jackknife_slow_rank(ranks, a, phi_obj, normalize = normalize)
  fast <- jackknife_fast_rank(ranks, a, phi_obj, normalize = normalize)
  err <- max(abs(fast$theta0 - ref$theta0), abs(fast$theta_loo - ref$theta_loo))
  if (!is.finite(err) || err > 1e-10) stop(sprintf("fast vs slow mismatch = %.3e", err), call. = FALSE)

  # --- benchmark ------------------------------------------------------------
  bm <- bench::mark(
    slow = jackknife_slow_rank(ranks, a, phi_obj, normalize = normalize),
    fast = jackknife_fast_rank(ranks, a, phi_obj, normalize = normalize),
    iterations = reps,
    check = FALSE
  )
  df <- data.frame(expression = as.character(bm$expression), n, p, s, q, reps,
    phi_dep, compiled_fast = compiled, min_seconds = as.numeric(bm$min),
    median_seconds = as.numeric(bm$median),
    iterations_per_second = as.numeric(bm[["itr/sec"]]),
    allocated_bytes = as.numeric(bm$mem_alloc),
    collections_per_second = as.numeric(bm[["gc/sec"]]))
  if (verbose) print(df)
  if (!is.null(csv)) {
    utils::write.csv(df, file = csv, row.names = FALSE)
    if (verbose) message("Wrote: ", csv)
  }
  invisible(df)
}
