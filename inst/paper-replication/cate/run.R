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
# This package-side script mirrors the manuscript-facing main production grid
# only. The chunked manifest layer, the nuisance-efficiency run, and the
# fixed-lambda sweep wrappers live in Cond_Z_Estimation/Codes/Simulations/.

library(RdistNN)

# ── Paper simulation grid ─────────────────────────────────────────────────────
# Mirrors the current manuscript-facing main grid in 05b_CATE_Simulation.tex:
# four retained designs, four sample sizes, and 10,000 replications per cell.
# The focused nuisance-efficiency design dgp6_nuisance_efficiency is implemented
# in the package but is intended for a separate full-sample-versus-cross-fit run.

PAPER_GRID <- list(
  designs     = c("dgp1_linear", "dgp2_nonlinear",
                  "dgp4_clustered", "dgp5_correlated"),
  n_vals      = c(100L, 250L, 500L, 1000L),
  kx_vals     = 2L,
  n_reps      = 10000L,
  master_seed = 42L,
  k           = 1L,  # one-nearest-neighbor DNN kernel used in the paper theory
  s_fn        = function(n, kx) floor(n^((kx + 1L) / (kx + 2L))),
  # Main manuscript comparison: same-sample DNN, cross-fitted DNN, fixed-
  # bandwidth local constant, and interaction OLS. Diagnostic oracle /
  # Mahalanobis variants remain available in the package but are not part of the
  # production table-and-figure run.
  methods     = c("dnn_fullsample", "dnn_crossfit",
                  "local_constant", "ols_interact"),
  alpha       = 0.05,
  sigma       = 0.5,
  clip        = 1e-3,
  n_folds     = 5L
)

# ── Run simulation ────────────────────────────────────────────────────────────

cat("Available methods:", paste(sim_cate_methods(), collapse = ", "), "\n\n")
cat("=== Running simulation ===\n")
cat("Methods:", paste(PAPER_GRID$methods, collapse = ", "), "\n\n")

# Parallelise across reps (optional):
#   future::plan(future::multisession)
raw <- sim_cate_run(
  designs     = PAPER_GRID$designs,
  n_vals      = PAPER_GRID$n_vals,
  kx_vals     = PAPER_GRID$kx_vals,
  n_reps      = PAPER_GRID$n_reps,
  master_seed = PAPER_GRID$master_seed,
  k           = PAPER_GRID$k,
  s_fn        = PAPER_GRID$s_fn,
  methods     = PAPER_GRID$methods,
  alpha       = PAPER_GRID$alpha,
  sigma       = PAPER_GRID$sigma,
  clip        = PAPER_GRID$clip,
  n_folds     = PAPER_GRID$n_folds,
  verbose     = TRUE
)

# ── Summarise ─────────────────────────────────────────────────────────────────

table_sim <- summarise_sim_results(raw)

cat("\n=== Simulation results ===\n")
print(table_sim, digits = 3)

# ── Save results ──────────────────────────────────────────────────────────────

out_dir <- file.path(getwd(), "results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

saveRDS(raw,       file.path(out_dir, "sim_cate_raw.rds"))
write.csv(table_sim, file.path(out_dir, "sim_cate_table.csv"),
          row.names = FALSE)

cat("\nResults saved to", out_dir, "\n")
