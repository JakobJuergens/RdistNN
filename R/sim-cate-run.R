# sim-cate-run.R
# Single-replicate estimator and Monte Carlo harness for the CATE simulation.
#
# Public API:
#   sim_cate_one_rep()  -- one replicate: draws data, runs all requested methods,
#                          returns a data.frame with one row per method
#   sim_cate_run()      -- full MC grid, returns a data.frame of raw results
#
# Available methods are listed by sim_cate_methods() and implemented in
# sim-cate-benchmarks.R.
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

#' Run one replicate of CATE estimators
#'
#' Draws one dataset from the named DGP (via [sim_cate_data()]), then runs each
#' requested method on the **same** dataset and returns per-replicate performance
#' metrics aligned with the paper's evaluation criteria.
#'
#' @param n Integer sample size.
#' @param kx Integer covariate dimension (\eqn{\ge 2}).
#' @param k Integer DNN neighbour count.
#' @param s Integer kernel order (undersmoothing parameter).
#' @param design Character DGP name; passed to [sim_cate_data()].
#' @param methods Character vector of method names; see [sim_cate_methods()].
#'   Default: `"dnn_oracle"`.
#' @param alpha Nominal coverage level (default 0.05).
#' @param seed Optional integer seed for this replicate.
#' @param sigma Noise SD forwarded to the DGP.
#' @param clip Propensity clipping bound forwarded to the DGP and cross-fitting.
#' @param n_folds Number of cross-fitting folds for the `dnn_crossfit` method.
#'
#' @return A `data.frame` with one row per method and columns
#'   `method`, `bias`, `sq_err`, `se`, `cover` (0/1), `ci_len`, `neff`.
#'   Failed methods have `NA` for all metric columns.
#'
#' @examples
#' rep1 <- sim_cate_one_rep(
#'   n = 200, kx = 2, k = 5, s = 12,
#'   design = "dgp2_nonlinear",
#'   methods = c("dnn_oracle", "ols_interact"),
#'   seed = 1
#' )
#' rep1
#'
#' @export
sim_cate_one_rep <- function(n, kx, k, s, design,
                              methods   = "dnn_oracle",
                              alpha     = 0.05,
                              seed      = NULL,
                              sigma     = 0.5,
                              clip      = 1e-3,
                              n_folds   = 5L) {
  dgp    <- sim_cate_data(design, n = n, kx = kx,
                          sigma = sigma, clip = clip, seed = seed)
  x_cols <- paste0("X", seq_len(kx))

  na_row <- c(tau_hat = NA_real_, tau0 = NA_real_,
              bias = NA_real_, sq_err = NA_real_, se = NA_real_,
              cover = NA_real_, ci_len = NA_real_, neff = NA_real_)

  rows <- lapply(methods, function(meth) {
    fn  <- .get_sim_method(meth)
    res <- tryCatch(
      fn(dgp = dgp, x_cols = x_cols, k = k, s = s, alpha = alpha,
         n_folds = n_folds, clip = clip),
      error = function(e) {
        warning(sprintf("sim_cate_one_rep: method '%s' failed: %s",
                        meth, conditionMessage(e)), call. = FALSE)
        NULL
      }
    )
    if (is.null(res)) {
      na_row
    } else {
      # Ensure neff is present even if method omits it
      if (!"neff" %in% names(res)) res <- c(res, neff = NA_real_)
      res[names(na_row)]
    }
  })

  do.call(rbind, lapply(seq_along(methods), function(i) {
    r <- rows[[i]]
    data.frame(
      method  = methods[i],
      tau_hat = r[["tau_hat"]],
      tau0    = r[["tau0"]],
      bias    = r[["bias"]],
      sq_err  = r[["sq_err"]],
      se      = r[["se"]],
      cover   = r[["cover"]],
      ci_len  = r[["ci_len"]],
      neff    = r[["neff"]],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))
}

# ── Monte Carlo grid ──────────────────────────────────────────────────────────

#' Run the CATE simulation grid
#'
#' Expands a grid of DGP designs, sample sizes, and covariate dimensions, then
#' runs `n_reps` independent Monte Carlo replications for each cell using
#' [sim_cate_one_rep()].  Each replicate draws one dataset and evaluates **all**
#' requested `methods` on the same draw, enabling like-for-like comparisons.
#' Returns a tidy data.frame of raw per-replicate, per-method results.  Use
#' [summarise_sim_results()] to aggregate to one row per `(method, design, n, kx)`.
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
#' Each replicate is wrapped in `tryCatch`.  A total rep failure (e.g. DGP
#' draw error) records `NA` for all metrics and all methods and emits a warning;
#' per-method failures are handled inside [sim_cate_one_rep()].
#' [summarise_sim_results()] reports `n_failed` per `(method, design, n, kx)`.
#'
#' @param designs Character vector of DGP names.
#' @param n_vals Integer vector of sample sizes.
#' @param kx_vals Integer vector of covariate dimensions.
#' @param n_reps Number of Monte Carlo replications per cell.
#' @param master_seed Integer seed used to draw the per-replicate seeds.
#'   Every cell uses the same seed vector so noise realisations are comparable
#'   across DGPs, sample sizes, and methods (fair comparison).
#' @param k Integer DNN neighbour count.
#' @param s_fn Function `f(n, kx)` returning the integer kernel order.
#'   Default satisfies the undersmoothing condition
#'   \eqn{s \asymp n^{(k_x+1)/(k_x+2)}}.
#' @param tuning Optional kernel-order tuning specification. When supplied,
#'   `s` is resolved per simulation cell from `tuning_eval(..., mode = "default")`
#'   and `s_fn` must be omitted.
#' @param methods Character vector of method names; see [sim_cate_methods()].
#'   Default: `"dnn_oracle"`.
#' @param alpha Nominal coverage level.
#' @param sigma Noise SD forwarded to the DGP.
#' @param clip Propensity clipping bound forwarded to the DGP and cross-fitting.
#' @param n_folds Number of cross-fitting folds for the `dnn_crossfit` method.
#' @param verbose Print one progress line per cell (default `TRUE`).
#'
#' @return A `data.frame` with one row per replicate-method combination and
#'   columns: `method`, `design`, `n`, `kx`, `s`, `seed`,
#'   `bias`, `sq_err`, `se`, `cover`, `ci_len`, `neff`.
#'   Failed entries have `NA` for all metric columns.
#'
#' @seealso [sim_cate_one_rep()], [summarise_sim_results()], [sim_cate_methods()]
#'
#' @examples
#' \dontrun{
#' # Single method (default):
#' raw <- sim_cate_run(designs = "dgp2_nonlinear", n_vals = 500L, n_reps = 10L)
#'
#' # Multiple methods on the same draws:
#' future::plan(future::multisession)
#' raw <- sim_cate_run(
#'   designs = "dgp2_nonlinear", n_vals = 500L, n_reps = 10L,
#'   methods = c("dnn_oracle", "dnn_crossfit", "ols_interact", "causal_forest")
#' )
#' future::plan(future::sequential)
#'
#' summarise_sim_results(raw)
#' }
#'
#' @export
sim_cate_run <- function(
    designs     = c("dgp1_linear", "dgp2_nonlinear", "dgp3_sparse_overlap",
                    "dgp4_clustered", "dgp5_correlated"),
    n_vals      = c(500L, 1000L, 2000L),
    kx_vals     = 2L,
    n_reps      = 500L,
    master_seed = 42L,
    k           = 5L,
    s_fn        = function(n, kx) floor(n^((kx + 1L) / (kx + 2L))),
    tuning      = NULL,
    methods     = "dnn_oracle",
    alpha       = 0.05,
    sigma       = 0.5,
    clip        = 1e-3,
    n_folds     = 5L,
    verbose     = TRUE
) {
  s_fn_arg <- if (missing(s_fn)) NULL else s_fn
  if (!is.null(tuning) && !is.null(s_fn_arg)) {
    .abort_bad_arg(
      "tuning/s_fn",
      "cannot provide both `tuning` and `s_fn`; omit `s_fn` when `tuning` is used",
      quoted = FALSE
    )
  }
  if (is.null(tuning) && is.null(s_fn_arg)) {
    s_fn <- function(n, kx) floor(n^((kx + 1L) / (kx + 2L)))
  } else if (is.null(tuning)) {
    s_fn <- s_fn_arg
  }

  # Rep-level parallel map: future_lapply when available, lapply otherwise.
  # Use future.seed = TRUE for parallel-safe reproducible RNG streams in the
  # workers. Replicate-specific randomness remains anchored by explicit unique
  # per-rep seeds passed as `seed = sd` into sim_cate_one_rep().
  .lapply <- if (requireNamespace("future.apply", quietly = TRUE)) {
    function(X, FUN) future.apply::future_lapply(X, FUN, future.seed = TRUE)
  } else {
    lapply
  }

  # All cells use the same seed vector so noise realisations are aligned
  # across DGPs, sample sizes, and methods (fair comparison).
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
    s  <- if (!is.null(tuning)) {
      tuning_eval(
        tuning = tuning,
        n = n,
        kx = kx,
        k = k,
        se = TRUE,
        jackknife = "fast",
        mode = "default"
      )$s
    } else {
      s_fn(n, kx)
    }

    if (verbose)
      message(sprintf("  design=%-22s  n=%4d  kx=%d  s=%d  methods=[%s]",
                      d, n, kx, s, paste(methods, collapse = ", ")))

    reps <- .lapply(as.list(seeds), function(sd) {
      result <- tryCatch(
        sim_cate_one_rep(n = n, kx = kx, k = k, s = s,
                         design = d, methods = methods,
                         alpha = alpha, sigma = sigma,
                         clip = clip, seed = sd, n_folds = n_folds),
        error = function(e) {
          warning(sprintf(
            "sim_cate_run: rep seed %d failed [design=%s n=%d kx=%d]: %s",
            sd, d, n, kx, conditionMessage(e)
          ), call. = FALSE)
          NULL
        }
      )

      if (is.null(result)) {
        # Total rep failure: NA rows for every method
        do.call(rbind, lapply(methods, function(meth) {
          data.frame(design = d, n = n, kx = kx, s = s, seed = sd,
                     method = meth, tau_hat = NA_real_, tau0 = NA_real_,
                     bias = NA_real_, sq_err = NA_real_, se = NA_real_,
                     cover = NA_real_, ci_len = NA_real_, neff = NA_real_,
                     stringsAsFactors = FALSE, row.names = NULL)
        }))
      } else {
        # Prepend cell identity columns (scalar assignment recycles safely)
        result$design <- d
        result$n      <- n
        result$kx     <- kx
        result$s      <- s
        result$seed   <- sd
        result[, c("design", "n", "kx", "s", "seed", "method",
                   "tau_hat", "tau0", "bias", "sq_err", "se",
                   "cover", "ci_len", "neff"),
               drop = FALSE]
      }
    })

    all_rows[[i]] <- do.call(rbind, reps)
  }

  do.call(rbind, all_rows)
}
