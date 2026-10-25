# Scalar distance benchmarks use the same providers as the public API.
rdistnn_bench_distance <- function(metric = c("l1", "l2"), n = 50000L, p = 8L,
                                    reps = 5L, seed = 1L, include_cpp_direct = TRUE,
                                    squared = FALSE) {
  rdistnn_bench_require()
  metric <- match.arg(metric)
  if (metric == "l1" && squared) stop("squared applies to L2 distances.", call. = FALSE)
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  x0 <- rnorm(p)
  flags <- if (metric == "l2") list(squared = squared) else list()
  builder <- do.call(get(paste0("distance_", metric)), flags)
  spec <- do.call(new_distance, builder)
  portable <- function() do.call(get(paste0(".dist_", metric, "_R")), c(list(X, x0), flags))
  public <- function() distance_eval(spec, X, x0)
  expected <- portable()
  rdistnn_bench_parity(expected, public(), 1e-12, "public distance")
  if (include_cpp_direct) {
    native <- get(paste0("rustat_dist_", metric))
    direct <- function() do.call(native, c(list(X, x0), flags))
    rdistnn_bench_parity(expected, direct(), 1e-12, "native distance")
    result <- bench::mark(R = portable(), API = public(), CPP = direct(),
                           iterations = reps, check = TRUE)
  } else {
    result <- bench::mark(R = portable(), API = public(), iterations = reps, check = TRUE)
  }
  rdistnn_bench_print(result)
  invisible(result)
}
