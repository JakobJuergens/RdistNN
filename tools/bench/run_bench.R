# Run from the repository root after installing RdistNN or devtools::load_all().
# Rscript tools/bench/run_bench.R l1 n=5000L reps=3L
# Rscript tools/bench/run_bench.R cate n=1000L s=8L n_subsets=128L reps=3L
rdistnn_bench <- function(which = c("cate", "plqr", "jackknife", "jk",
                                    "dist_l1", "l1", "dist_l2", "l2"), ...) {
  which <- match.arg(which)
  if (!requireNamespace("RdistNN", quietly = TRUE)) {
    stop("Install RdistNN or run devtools::load_all() before benchmarking.", call. = FALSE)
  }
  env <- new.env(parent = asNamespace("RdistNN"))
  sys.source("tools/bench/utils-bench.R", envir = env)
  if (which %in% c("dist_l1", "l1", "dist_l2", "l2")) {
    sys.source("tools/bench/bench_distances.R", envir = env)
    metric <- if (which %in% c("dist_l1", "l1")) "l1" else "l2"
    return(do.call(env$rdistnn_bench_distance, c(list(metric = metric), list(...))))
  }
  module <- switch(which, cate = "bench_score_cate.R", plqr = "bench_score_plqr.R",
                   jackknife = "bench_jackknife.R", jk = "bench_jackknife.R")
  function_name <- switch(which, cate = "rdistnn_bench_cate", plqr = "rdistnn_bench_plqr",
                          jackknife = "rdistnn_bench_jackknife", jk = "rdistnn_bench_jackknife")
  sys.source(file.path("tools/bench", module), envir = env)
  do.call(get(function_name, envir = env), list(...))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args)) stop("Choose cate, plqr, jackknife, l1 or l2.", call. = FALSE)
  dots <- list()
  for (part in args[-1L]) {
    match <- regexec("^([A-Za-z][A-Za-z0-9_.]*)=(.*)$", part)
    pair <- regmatches(part, match)[[1L]]
    if (!length(pair)) stop("Arguments must use name=value: ", part, call. = FALSE)
    value <- pair[[3L]]
    if (grepl("^[+-]?[0-9]+L$", value)) value <- as.integer(sub("L$", "", value))
    else if (grepl("^[+-]?([0-9]+[.]?[0-9]*|[.][0-9]+)([eE][+-]?[0-9]+)?$", value)) value <- as.numeric(value)
    else if (value %in% c("TRUE", "FALSE")) value <- as.logical(value)
    dots[[pair[[2L]]]] <- value
  }
  do.call(rdistnn_bench, c(list(which = args[[1L]]), dots))
}
