# inst/bench/bench_dist_l2.R
# Benchmarks for dist_l2() and order_by_l2(), with strict parity.
# Runs two separate marks: (1) distances, (2) ordering.

# --- Minimal helpers bootstrap ------------------------------------------------
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
  if (!exists("rdistnn_bench_print", inherits = TRUE)) {
    rdistnn_bench_print <- function(b) {
      print(b[, c("expression","min","median","itr/sec","mem_alloc","gc/sec")])
    }
  }
}

# --- Benchmark function -------------------------------------------------------
rdistnn_bench_dist_l2 <- function(
    n = 50000L,
    p = 8L,
    reps = 5L,
    seed = 1L,
    squared = FALSE,
    bench_order = TRUE,
    include_cpp_direct = TRUE
) {
  rdistnn_bench_require()
  set.seed(seed)

  # Data & query
  X  <- matrix(rnorm(n * p), n, p)
  x0 <- rnorm(p)

  # Pure-R baselines
  dx_R <- function() .dist_l2_R(X, x0, squared = squared)
  ord_R <- function() {
    dx <- .dist_l2_R(X, x0, squared = squared)
    as.integer(order(dx, seq_along(dx)))
  }

  # Public API (uses C++ if available)
  dx_API  <- function() dist_l2(X, x0, squared = squared)
  ord_API <- function() order_by_l2(X, x0, squared = squared)$ord

  # Direct C++ availability
  has_cpp <- exists("rustat_dist_l2", mode = "function") &&
    exists("rustat_order_by_l2", mode = "function")
  dx_CPP  <- if (has_cpp && include_cpp_direct) function() rustat_dist_l2(X, x0, squared) else NULL
  ord_CPP <- if (has_cpp && include_cpp_direct) function() rustat_order_by_l2(X, x0, squared)$ord else NULL

  # ----- Parity checks (fail fast) -------------------------------------------
  # Distances
  dx_r <- dx_R(); dx_a <- dx_API()
  if (!isTRUE(all.equal(dx_r, dx_a, tolerance = 1e-12)))
    stop("dist_l2 parity failure: API vs .dist_l2_R differ beyond 1e-12.")
  if (!is.null(dx_CPP)) {
    dx_c <- dx_CPP()
    if (!isTRUE(all.equal(dx_r, dx_c, tolerance = 1e-12)))
      stop("dist_l2 parity failure: C++ vs .dist_l2_R differ beyond 1e-12.")
  }

  # Ordering
  if (bench_order) {
    or_r <- ord_R(); or_a <- ord_API()
    if (!identical(or_r, or_a))
      stop("order_by_l2 parity failure: API ord differs from R baseline.")
    if (!is.null(ord_CPP)) {
      or_c <- ord_CPP()
      if (!identical(or_r, or_c))
        stop("order_by_l2 parity failure (C++): ord differs from R baseline.")
    }
    # Property: L2 and L2^2 induce identical order
    dx_lin <- .dist_l2_R(X, x0, squared = FALSE)
    dx_sq  <- .dist_l2_R(X, x0, squared = TRUE)
    if (!identical(order(dx_lin, seq_along(dx_lin)),
                   order(dx_sq,  seq_along(dx_sq))))
      stop("Ordering property failed: L2 and L2^2 must induce identical order.")
  }

  # ----- Timings --------------------------------------------------------------
  # Distances (all return the same numeric vector; safe to check=TRUE)
  b_dist <- if (!is.null(dx_CPP)) {
    bench::mark(
      L2_R_dist   = .dist_l2_R(X, x0, squared = squared),
      L2_API_dist = dist_l2(X, x0, squared = squared),
      L2_CPP_dist = rustat_dist_l2(X, x0, squared),
      iterations = reps, check = TRUE
    )
  } else {
    bench::mark(
      L2_R_dist   = .dist_l2_R(X, x0, squared = squared),
      L2_API_dist = dist_l2(X, x0, squared = squared),
      iterations = reps, check = TRUE
    )
  }

  # Ordering (all return the same integer vector; safe to check=TRUE)
  if (bench_order) {
    b_ord <- if (!is.null(ord_CPP)) {
      bench::mark(
        L2_R_order   = { dx <- .dist_l2_R(X, x0, squared); order(dx, seq_along(dx)) },
        L2_API_order = order_by_l2(X, x0, squared = squared)$ord,
        L2_CPP_order = rustat_order_by_l2(X, x0, squared = squared)$ord,
        iterations = reps, check = TRUE
      )
    } else {
      bench::mark(
        L2_R_order   = { dx <- .dist_l2_R(X, x0, squared); order(dx, seq_along(dx)) },
        L2_API_order = order_by_l2(X, x0, squared = squared)$ord,
        iterations = reps, check = TRUE
      )
    }
  }

  # ----- Print ---------------------------------------------------------------
  cat("\nDistances (squared =", squared, ")\n", sep = " ")
  rdistnn_bench_print(b_dist)
  if (bench_order) {
    cat("\nOrdering (squared =", squared, ")\n", sep = " ")
    rdistnn_bench_print(b_ord)
  }

  invisible(if (bench_order) list(dist = b_dist, order = b_ord) else b_dist)
}
