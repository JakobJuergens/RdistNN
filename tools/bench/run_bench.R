# tools/bench/run_bench.R
# Usage from R:
#   devtools::load_all()
#   source("tools/bench/run_bench.R")
#   rdistnn_bench("cate", n = 5000)
#   rdistnn_bench("plqr", reps = 10)
#
# Usage from shell:
#   Rscript tools/bench/run_bench.R cate n=5000 reps=10

rdistnn_bench <- function(
    which = c("cate", "plqr", "jackknife", "jk", "dist_l1", "l1", "dist_l2", "l2"), ...)
  {
  which <- match.arg(which)
  # Ensure utils are available when sourced directly
  if (!exists("rdistnn_bench_require")) {
    source(file.path("tools", "bench", "utils-bench.R"), local = TRUE)
  }

  switch(which,
         cate = {
           source(file.path("tools", "bench", "bench_cate.R"), local = TRUE)
           rdistnn_bench_cate(...)
         },
         plqr = {
           source(file.path("tools", "bench", "bench_plqr.R"), local = TRUE)
           rdistnn_bench_plqr(...)
         },
         jackknife = , jk = {
           source(file.path("tools", "bench", "bench_jackknife.R"), local = TRUE)
           rdistnn_bench_jackknife(...)
         },
         dist_l2 = , l2 = {
           source(file.path("tools", "bench", "bench_dist_l2.R"), local = TRUE)
           rdistnn_bench_dist_l2(...)
         },
         dist_l1 = , l1 = {
           source(file.path("tools", "bench", "bench_dist_l1.R"), local = TRUE)
           rdistnn_bench_dist_l1(...)
         }
  )
}

# If called via Rscript, parse simple key=value args
if (identical(environmentName(topenv()), "R_GlobalEnv")) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args)) {
    which <- args[[1]]
    dots <- list()
    if (length(args) > 1) {
      # parse key=value pairs into a list
      kv <- strsplit(args[-1], "=", fixed = TRUE)
      for (p in kv) {
        key <- p[[1]]
        val <- p[[2]]
        # try numeric, integer, logical, otherwise keep character
        if (grepl("^[0-9]+L$", val)) {
          val <- as.integer(sub("L$", "", val))
        } else if (grepl("^[0-9.]+$", val)) {
          val <- as.numeric(val)
        } else if (val %in% c("TRUE","FALSE")) {
          val <- as.logical(val)
        }
        dots[[key]] <- val
      }
    }
    rdistnn_bench(which, !!!dots)
  }
}
