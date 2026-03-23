# sim-cate-run.R
# Single-replicate estimator and Monte Carlo harness for the CATE simulation.
#
# Implements the proposed DNN Z-estimator with oracle nuisances.
# Competitor estimators live in inst/paper-replication/cate/benchmarks.R.
#
# Public API:
#   sim_cate_one_rep()  -- one replicate, returns a named numeric vector
#   sim_cate_run()      -- full MC grid, returns a data.frame of raw results
#
# Parallelism:
#   sim_cate_run() uses future.apply::future_lapply() when the future.apply
#   package is installed, falling back silently to lapply otherwise.
#   Users control the execution backend by calling future::plan() before
#   sim_cate_run():
#
#     # Personal computer (all local cores):
#     future::plan(future::multisession)
#     raw <- sim_cate_run(...)
#
#     # SLURM single-node (cores allocated via --cpus-per-task):
#     ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 1))
#     future::plan(future::multisession, workers = ncores)
#     raw <- sim_cate_run(...)
#
#   For multi-node SLURM array parallelism, see
#   inst/paper-replication/cate/hpc/task.R and submit.sh, which dispatch
#   one sim_cate_run() call per cell via SLURM_ARRAY_TASK_ID.

# ── Single replicate ──────────────────────────────────────────────────────────

#' Run one replicate of the DNN CATE estimator
#'
#' Draws one dataset from the named DGP (via [sim_cate_data()]), fits the DNN
#' Z-estimator with oracle nuisances, and returns per-replicate performance
#' metrics aligned with the paper's evaluation criteria.
#'
#' @param n Integer sample size.
#' @param kx Integer covariate dimension (\eqn{\ge 2}).
#' @param k Integer DNN neighbour count.
#' @param s Integer kernel order (undersmoothing parameter).
#' @param design Character DGP name; passed to [sim_cate_data()].
#' @param alpha Nominal coverage level (default 0.05).
#' @param seed Optional integer seed for this replicate.
#' @param sigma Noise SD forwarded to the DGP.
#' @param clip Propensity clipping bound forwarded to the DGP.
#'
#' @return Named numeric vector with elements `bias`, `sq_err`, `se`,
#'   `cover` (0/1), `ci_len`, `neff`.
#'
#' @examples
#' rep1 <- sim_cate_one_rep(
#'   n = 200, kx = 2, k = 5, s = 12,
#'   design = "dgp2_nonlinear", seed = 1
#' )
#' rep1
#'
#' @export
sim_cate_one_rep <- function(n, kx, k, s, design,
                              alpha = 0.05, seed = NULL,
                              sigma = 0.5, clip = 1e-3) {
  dgp    <- sim_cate_data(design, n = n, kx = kx,
                          sigma = sigma, clip = clip, seed = seed)
  tau0   <- dgp$tau_x0
  dat    <- dgp$data
  x0     <- dgp$x0
  x_cols <- paste0("X", seq_len(kx))

  fit <- cate_dnn(
    data      = dat,
    x_cols    = x_cols,
    y_col     = "Y",
    t_col     = "D",
    e_col     = "e_oracle",
    m0_col    = "m0_oracle",
    m1_col    = "m1_oracle",
    x0        = x0,
    k         = k,
    s         = s,
    se        = TRUE,
    jackknife = "fast",
    return    = "ess"
  )

  tau_hat <- fit$tau
  se_hat  <- fit$se
  neff    <- fit$ess
  z       <- stats::qnorm(1 - alpha / 2)
  covered <- (tau_hat - z * se_hat <= tau0) & (tau0 <= tau_hat + z * se_hat)

  c(bias   = tau_hat - tau0,
    sq_err = (tau_hat - tau0)^2,
    se     = se_hat,
    cover  = as.numeric(covered),
    ci_len = 2 * z * se_hat,
    neff   = neff)
}

# ── Monte Carlo grid ──────────────────────────────────────────────────────────

#' Run the CATE simulation grid
#'
#' Expands a grid of DGP designs, sample sizes, and covariate dimensions, then
#' runs `n_reps` independent Monte Carlo replications for each cell using
#' [sim_cate_one_rep()].  Returns a tidy data.frame of raw per-replicate
#' results.  Use [summarise_sim_results()] to aggregate to one row per cell.
#'
#' @section Parallelism:
#' When the `future.apply` package is installed, the rep loop within each cell
#' is parallelised via `future.apply::future_lapply()`.  Control the backend
#' by calling `future::plan()` before `sim_cate_run()`:
#'
#' ```r
#' # All local cores (PC or single SLURM node):
#' future::plan(future::multisession)
#' raw <- sim_cate_run(...)
#' ```
#'
#' Without `future.apply` the function falls back silently to sequential
#' `lapply`.  For multi-node HPC, see
#' `system.file("paper-replication/cate/hpc", package = "RdistNN")`.
#'
#' @section Fault tolerance:
#' Each replicate is wrapped in `tryCatch`.  A failed rep records `NA` for all
#' metrics and emits a warning; the rest of the run continues unaffected.
#' [summarise_sim_results()] reports `n_failed` per cell.
#'
#' @param designs Character vector of DGP names.
#' @param n_vals Integer vector of sample sizes.
#' @param kx_vals Integer vector of covariate dimensions.
#' @param n_reps Number of Monte Carlo replications per cell.
#' @param master_seed Integer seed used to draw the per-replicate seeds.
#'   Every cell uses the same seed vector so noise realisations are comparable
#'   across DGPs and sample sizes.
#' @param k Integer DNN neighbour count.
#' @param s_fn Function `f(n, kx)` returning the integer kernel order.
#'   Default satisfies the undersmoothing condition
#'   \eqn{s \asymp n^{(k_x+1)/(k_x+2)}}.
#' @param alpha Nominal coverage level.
#' @param sigma Noise SD forwarded to the DGP.
#' @param clip Propensity clipping bound forwarded to the DGP.
#' @param verbose Print one progress line per cell (default `TRUE`).
#'
#' @return A `data.frame` with one row per replicate and columns:
#'   `design`, `n`, `kx`, `s`, `seed`, `bias`, `sq_err`, `se`,
#'   `cover`, `ci_len`, `neff`.  Failed reps have `NA` for all metric columns.
#'
#' @seealso [sim_cate_one_rep()], [summarise_sim_results()]
#'
#' @examples
#' \dontrun{
#' # Sequential (no extra packages needed):
#' raw <- sim_cate_run(designs = "dgp2_nonlinear", n_vals = 500L, n_reps = 10L)
#'
#' # Parallel on a PC:
#' future::plan(future::multisession)
#' raw <- sim_cate_run(designs = "dgp2_nonlinear", n_vals = 500L, n_reps = 10L)
#' future::plan(future::sequential)   # reset
#'
#' summarise_sim_results(raw)
#' }
#'
#' @export
sim_cate_run <- function(
    designs     = c("dgp1_linear", "dgp2_nonlinear", "dgp3_sparse_overlap"),
    n_vals      = c(500L, 1000L, 2000L),
    kx_vals     = 2L,
    n_reps      = 500L,
    master_seed = 42L,
    k           = 5L,
    s_fn        = function(n, kx) floor(n^((kx + 1L) / (kx + 2L))),
    alpha       = 0.05,
    sigma       = 0.5,
    clip        = 1e-3,
    verbose     = TRUE
) {
  # Rep-level parallel map: future_lapply when available, lapply otherwise.
  # We set future.seed = FALSE because seeds are managed explicitly inside
  # sim_cate_one_rep() via sim_cate_data(..., seed = sd).
  .lapply <- if (requireNamespace("future.apply", quietly = TRUE)) {
    function(X, FUN) future.apply::future_lapply(X, FUN, future.seed = FALSE)
  } else {
    lapply
  }

  # All cells use the same seed vector so noise realisations are aligned
  # across DGPs and sample sizes (fair comparison).
  set.seed(master_seed)
  seeds <- sample.int(1e6L, n_reps)

  cells <- expand.grid(
    design = designs,
    n      = n_vals,
    kx     = kx_vals,
    stringsAsFactors = FALSE
  )

  all_rows <- vector("list", nrow(cells))

  for (i in seq_len(nrow(cells))) {
    d  <- cells$design[i]
    n  <- cells$n[i]
    kx <- cells$kx[i]
    s  <- s_fn(n, kx)

    if (verbose)
      message(sprintf("  design=%-22s  n=%4d  kx=%d  s=%d", d, n, kx, s))

    reps <- .lapply(as.list(seeds), function(sd) {
      result <- tryCatch(
        sim_cate_one_rep(n = n, kx = kx, k = k, s = s,
                         design = d, alpha = alpha,
                         sigma = sigma, clip = clip, seed = sd),
        error = function(e) {
          warning(sprintf(
            "sim_cate_run: rep seed %d failed [design=%s n=%d kx=%d]: %s",
            sd, d, n, kx, conditionMessage(e)
          ), call. = FALSE)
          NULL
        }
      )

      if (is.null(result)) {
        data.frame(design = d, n = n, kx = kx, s = s, seed = sd,
                   bias = NA_real_, sq_err = NA_real_, se = NA_real_,
                   cover = NA_real_, ci_len = NA_real_, neff = NA_real_,
                   stringsAsFactors = FALSE)
      } else {
        data.frame(design = d, n = n, kx = kx, s = s, seed = sd,
                   bias   = result[["bias"]],
                   sq_err = result[["sq_err"]],
                   se     = result[["se"]],
                   cover  = result[["cover"]],
                   ci_len = result[["ci_len"]],
                   neff   = result[["neff"]],
                   stringsAsFactors = FALSE)
      }
    })

    all_rows[[i]] <- do.call(rbind, reps)
  }

  do.call(rbind, all_rows)
}
