# inst/paper-replication/cate/run.R
#
# Full paper replication script for the CATE simulation study.
# Reference: "Robust Inference for Conditional Z-Estimation via
#            Distributional Nearest Neighbors" (Juergens)
#            Paper/Parts/05_Applications/05b_CATE_Simulation.tex
#
# Assumes RdistNN is installed.  Run from any working directory:
#
#   Rscript system.file("paper-replication/cate/run.R", package = "RdistNN")
#
# or interactively:
#
#   source(system.file("paper-replication/cate/run.R", package = "RdistNN"))
#
# Results are written to ./results/sim_cate_raw.rds and
# ./results/sim_cate_table.csv in the current working directory.

library(RdistNN)

# ── Paper simulation grid ─────────────────────────────────────────────────────
# Mirrors the design grid in 05b_CATE_Simulation.tex.
# Extend kx_vals to c(1L, 2L) once p=1 DGPs are finalised.

PAPER_GRID <- list(
  designs     = c("dgp1_linear", "dgp2_nonlinear", "dgp3_sparse_overlap"),
  n_vals      = c(500L, 1000L, 2000L),
  kx_vals     = 2L,
  n_reps      = 500L,
  master_seed = 42L,
  k           = 5L,
  s_fn        = function(n, kx) floor(n^((kx + 1L) / (kx + 2L))),
  alpha       = 0.05,
  sigma       = 0.5,
  clip        = 1e-3
)

# ── Run DNN estimator ─────────────────────────────────────────────────────────

cat("=== DNN Z-estimator (oracle nuisances) ===\n")

raw_dnn <- sim_cate_run(
  designs     = PAPER_GRID$designs,
  n_vals      = PAPER_GRID$n_vals,
  kx_vals     = PAPER_GRID$kx_vals,
  n_reps      = PAPER_GRID$n_reps,
  master_seed = PAPER_GRID$master_seed,
  k           = PAPER_GRID$k,
  s_fn        = PAPER_GRID$s_fn,
  alpha       = PAPER_GRID$alpha,
  sigma       = PAPER_GRID$sigma,
  clip        = PAPER_GRID$clip,
  verbose     = TRUE
)

table_dnn <- summarise_sim_results(raw_dnn)

cat("\n=== DNN results ===\n")
print(table_dnn, digits = 3)

# ── Run benchmark estimators (optional) ──────────────────────────────────────
# Uncomment after implementing benchmarks.R.
#
# bench_script <- system.file(
#   "paper-replication/cate/benchmarks.R", package = "RdistNN"
# )
# source(bench_script)
#
# raw_bench  <- run_benchmarks(PAPER_GRID)
# table_comp <- merge(table_dnn, raw_bench, by = c("design", "n", "kx"),
#                     suffixes = c("_dnn", "_bench"))

# ── Save results ──────────────────────────────────────────────────────────────

out_dir <- file.path(getwd(), "results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

saveRDS(raw_dnn,    file.path(out_dir, "sim_cate_raw.rds"))
write.csv(table_dnn, file.path(out_dir, "sim_cate_table.csv"),
          row.names = FALSE)

cat("\nResults saved to", out_dir, "\n")
