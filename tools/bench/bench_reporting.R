# Run in a fresh process for each explicitly selected installation:
# Rscript tools/bench/bench_reporting.R LIBRARY OUTPUT_DIR LABEL [REPS]
# Compare numerical.rds exactly before comparing measurements.csv. No timing
# threshold belongs in CI; invocation counts and numeric equality are contracts.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("Supply LIBRARY OUTPUT_DIR LABEL [REPS].")
lib <- normalizePath(args[1L], mustWork = TRUE)
out_dir <- args[2L]
label <- args[3L]
reps <- if (length(args) > 3L) as.integer(args[4L]) else 5L
stopifnot(length(reps) == 1L, is.finite(reps), reps > 0L)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
library(RdistNN)
stopifnot(identical(normalizePath(dirname(find.package("RdistNN"))), lib))
options(warn = 1, mc.cores = 2L)
ns <- asNamespace("RdistNN")
set.seed(8124)

measure <- function(fn, quiet, id) {
  run <- if (quiet) function() suppressWarnings(fn()) else fn
  run() # warm code paths and dependencies before timing
  times <- replicate(reps, system.time(run())[["elapsed"]])
  trace <- tempfile("reporting-allocation-")
  Rprofmem(trace)
  value <- tryCatch(run(), finally = Rprofmem(NULL))
  allocations <- suppressWarnings(as.numeric(sub(" .*", "", readLines(trace, warn = FALSE))))
  unlink(trace)
  list(row = data.frame(label, case = id, quiet, median_seconds = median(times),
    allocated_bytes = sum(allocations, na.rm = TRUE), allocations = sum(!is.na(allocations))),
    value = value)
}

# Instrument only a separate untimed execution, preserving original bindings.
count_calls <- function(fn) {
  names <- c("distance_eval", "distance_eval_batch", ".dnn_distance_groups",
    ".dnn_prepare_distance_groups", ".dnn_compute_group_facts",
    "rustat_dist_l2", "rustat_dist_l2_batch", "rustat_order_by_l2")
  names <- names[vapply(names, exists, logical(1L), envir = ns, inherits = FALSE)]
  originals <- mget(names, ns, inherits = FALSE)
  counts <- setNames(integer(length(names)), names)
  for (name in names) {
    wrapper <- local({
      key <- name
      original <- originals[[name]]
      function(...) { counts[key] <<- counts[key] + 1L; original(...) }
    })
    unlockBinding(name, ns)
    assign(name, wrapper, ns)
    lockBinding(name, ns)
  }
  on.exit(for (name in names) {
    unlockBinding(name, ns); assign(name, originals[[name]], ns); lockBinding(name, ns)
  })
  suppressWarnings(fn())
  counts
}

rows <- counters <- numerical <- list()
for (native in c(FALSE, TRUE)) {
  options(RdistNN.use_compiled = native)
  if (native && !is.loaded("_RdistNN_rustat_dist_l2_batch", PACKAGE = "RdistNN")) {
    stop("The selected installed package lacks its required compiled distance engine.")
  }
  for (size in list(c(64L, 1L), c(2000L, 8L), c(32000L, 16L))) {
    n <- size[1L]; m <- size[2L]
    for (tied in c(FALSE, TRUE)) {
      X <- matrix(if (tied) rep(seq_len(n %/% 4L), each = 4L) else seq_len(n), ncol = 1L)
      X0 <- matrix(seq(.125, .375, length.out = m), ncol = 1L)
      Y <- sin(seq_len(n) / 10)
      fit_call <- function() np_dnn(X = X, X0 = X0, Y = Y, s = 8L, se = FALSE,
        return = c("W", "ranks", "dx"))
      id <- sprintf("fit_n%d_m%d_tied%s_native%s", n, m, tied, native)
      counts <- count_calls(fit_call)
      counters[[id]] <- data.frame(case = id, operation = names(counts), calls = unname(counts))
      for (quiet in c(FALSE, TRUE)) {
        res <- measure(fit_call, quiet, id)
        rows[[length(rows) + 1L]] <- res$row
      }
      f <- res$value
      numerical[[id]] <- list(mu = f$mu, W = f$W, ranks = f$.model$ranks, groups = f$.model$groups)
    }
  }
  X <- matrix(seq_len(2000L)); X0 <- matrix(seq(.125, .375, length.out = 8L))
  Y <- sin(seq_len(nrow(X)) / 10)
  grid_spec <- new_tuning("reporting_benchmark", eval = function(...) list(s = 8L),
    eval_grid = function(...) data.frame(s = 2:12))
  grid_call <- function() np_dnn(X = X, X0 = X0, Y = Y, tuning = grid_spec,
    tuning_mode = "grid", se = FALSE, return = c("W", "ranks", "dx"))
  grid <- suppressWarnings(grid_call())
  cases <- list(grid = grid_call)
  if (requireNamespace("ggplot2", quietly = TRUE)) cases$dashboard <- function() plot_dashboard(grid)
  for (name in names(cases)) {
    id <- paste0(name, "_native", native)
    counts <- count_calls(cases[[name]])
    counters[[id]] <- data.frame(case = id, operation = names(counts), calls = unname(counts))
    for (quiet in c(FALSE, TRUE)) {
      res <- measure(cases[[name]], quiet, id)
      rows[[length(rows) + 1L]] <- res$row
    }
    if (name == "grid") numerical[[id]] <- res$value$tuning_grid[c("candidates", "estimates", "diagnostics")]
    else numerical[[id]] <- lapply(if (is.list(res$value) && !inherits(res$value, "ggplot")) res$value else list(res$value),
      function(panel) if (inherits(panel, "ggplot")) panel$data else NULL)
  }
}
write.csv(do.call(rbind, rows), file.path(out_dir, "measurements.csv"), row.names = FALSE)
write.csv(do.call(rbind, counters), file.path(out_dir, "invocations.csv"), row.names = FALSE)
saveRDS(numerical, file.path(out_dir, "numerical.rds"))
writeLines(c(capture.output(sessionInfo()), paste("Library:", lib), paste("Label:", label)),
  file.path(out_dir, "environment.txt"))
