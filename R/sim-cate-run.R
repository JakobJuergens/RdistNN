# sim-cate-run.R
# Single-replicate estimator and Monte Carlo harness for the CATE simulation.
#
# Public API:
#   sim_cate_one_rep()  -- one replicate: draws data, runs all requested methods,
#                          returns a data.frame with one or more rows per method
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

# ── Shared row helpers ────────────────────────────────────────────────────────

.sim_cate_raw_schema <- function() {
  c(
    "method", "design", "n", "kx", "k", "tuning_source",
    "candidate_id", "lambda", "s", "seed",
    .sim_cate_metric_names(),
    "status"
  )
}

.sim_cate_is_dnn_method <- function(method) {
  grepl("^dnn_", method)
}

.sim_cate_tuning_default <- function(tuning, n, kx, k) {
  tuning_eval(
    tuning = tuning,
    n = n,
    kx = kx,
    k = k,
    se = TRUE,
    jackknife = "fast",
    mode = "default"
  )
}

.sim_cate_grid_candidates <- function(tuning, n, kx, k) {
  resolved <- .dnn_resolve_s(
    s = NULL,
    tuning = tuning,
    n = n,
    p = kx,
    k = k,
    se = TRUE,
    jackknife = "fast"
  )
  resolved$tuning$grid$candidates
}

.sim_cate_failure_rows <- function(method, n, kx, k, s, tuning, tuning_mode,
                                   alpha, by_candidate = FALSE) {
  if (isTRUE(by_candidate) && identical(tuning_mode, "grid") &&
      .sim_cate_is_dnn_method(method)) {
    candidates <- .sim_cate_grid_candidates(tuning = tuning, n = n, kx = kx, k = k)
    return(.sim_metric(
      tau_hat = rep(NA_real_, nrow(candidates)),
      se_hat = rep(NA_real_, nrow(candidates)),
      tau0 = NA_real_,
      alpha = alpha,
      neff = rep(NA_real_, nrow(candidates)),
      candidate_id = candidates$candidate_id,
      lambda = candidates$lambda,
      s = candidates$s,
      status = "failed"
    ))
  }

  .sim_metric(
    tau_hat = NA_real_,
    se_hat = NA_real_,
    tau0 = NA_real_,
    alpha = alpha,
    s = if (.sim_cate_is_dnn_method(method)) s else NA_integer_,
    status = "failed"
  )
}

.sim_cate_normalize_method_rows <- function(res) {
  if (!is.data.frame(res)) {
    res <- as.list(res)
    res <- as.data.frame(res, check.names = FALSE, stringsAsFactors = FALSE)
  }
  missing <- setdiff(.sim_cate_row_names(), names(res))
  if (length(missing) > 0L) {
    for (col in missing) {
      res[[col]] <- if (identical(col, "status")) "ok" else NA_real_
    }
  }
  res[, .sim_cate_row_names(), drop = FALSE]
}

.sim_cate_cell_config <- function(tuning, tuning_mode, s_fn, n, kx, k) {
  if (identical(tuning_mode, "grid")) {
    candidates <- .sim_cate_grid_candidates(tuning = tuning, n = n, kx = kx, k = k)
    return(list(
      s = NA_integer_,
      tuning_source = "tuning",
      candidates = candidates
    ))
  }

  if (!is.null(tuning)) {
    default <- .sim_cate_tuning_default(tuning = tuning, n = n, kx = kx, k = k)
    return(list(
      s = as.integer(default$s),
      tuning_source = "tuning",
      candidates = NULL
    ))
  }

  list(
    s = as.integer(s_fn(n, kx)),
    tuning_source = "s_fn",
    candidates = NULL
  )
}

.sim_cate_add_cell_columns <- function(result, method, d, n, kx, k,
                                       tuning_source, seed) {
  result$method <- method
  result$design <- d
  result$n <- n
  result$kx <- kx
  result$k <- k
  result$tuning_source <- tuning_source
  result$seed <- seed
  result[, .sim_cate_raw_schema(), drop = FALSE]
}

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
#' @param s Integer kernel order (undersmoothing parameter). May be omitted
#'   when `tuning` is supplied.
#' @param design Character DGP name; passed to [sim_cate_data()].
#' @param methods Character vector of method names; see [sim_cate_methods()].
#'   Default: `"dnn_oracle"`.
#' @param tuning Optional kernel-order tuning specification.
#' @param tuning_mode One of `"default"` or `"grid"`. `"default"` evaluates
#'   one scalar kernel order. `"grid"` requires `tuning` and emits one row per
#'   DNN candidate while non-DNN benchmarks remain one row per replicate.
#' @param alpha Nominal coverage level (default 0.05).
#' @param seed Optional integer seed for this replicate.
#' @param sigma Noise SD forwarded to the DGP.
#' @param clip Propensity clipping bound forwarded to the DGP and cross-fitting.
#' @param n_folds Number of cross-fitting folds for the `dnn_crossfit` method.
#'
#' @return A `data.frame` with one or more rows per method and columns
#'   `method`, `candidate_id`, `lambda`, `s`, `tau_hat`, `tau0`, `bias`,
#'   `sq_err`, `se`, `cover` (0/1), `ci_len`, `neff`, and `status`. Methods
#'   with explicit localization weights also report the realized localized
#'   truth, branch drift, estimator remainder, studentized error, physical
#'   90-percent-mass radius, treatment-arm effective sample sizes, and the
#'   top-five score-contribution share.
#'   Failed methods have `NA` for all metric columns and `status = "failed"`.
#'
#' @examples
#' rep1 <- sim_cate_one_rep(
#'   n = 200, kx = 2, k = 1, s = 12,
#'   design = "dgp2_nonlinear",
#'   methods = c("dnn_fullsample", "ols_interact"),
#'   seed = 1
#' )
#' rep1
#'
#' @export
sim_cate_one_rep <- function(n, kx, k, s, design,
                              methods   = "dnn_oracle",
                              tuning    = NULL,
                              tuning_mode = c("default", "grid"),
                              alpha     = 0.05,
                              seed      = NULL,
                              sigma     = 0.5,
                              clip      = 1e-3,
                              n_folds   = 5L) {
  s_missing <- missing(s)
  if (s_missing) {
    s <- NULL
  }
  tuning_mode <- match.arg(tuning_mode)
  if (identical(tuning_mode, "grid") && is.null(tuning)) {
    .abort_bad_arg(
      "tuning_mode",
      "requires `tuning` when `tuning_mode = \"grid\"`",
      quoted = FALSE
    )
  }
  if (!identical(tuning_mode, "grid") && is.null(s) && is.null(tuning)) {
    .abort_bad_arg(
      "s",
      "must be supplied unless `tuning` is supplied",
      quoted = FALSE
    )
  }
  if (!identical(tuning_mode, "grid") && is.null(s) && !is.null(tuning)) {
    s <- .sim_cate_tuning_default(tuning = tuning, n = n, kx = kx, k = k)$s
  }

  dgp    <- sim_cate_data(design, n = n, kx = kx,
                          sigma = sigma, clip = clip, seed = seed)
  x_cols <- paste0("X", seq_len(kx))
  fns <- lapply(methods, .get_sim_method)

  rows <- lapply(seq_along(methods), function(i) {
    meth <- methods[[i]]
    fn <- fns[[i]]
    res <- tryCatch(
      fn(dgp = dgp, x_cols = x_cols, k = k, s = s,
         tuning = tuning, tuning_mode = tuning_mode, alpha = alpha,
         n_folds = n_folds, clip = clip),
      error = function(e) {
        warning(sprintf("sim_cate_one_rep: method '%s' failed: %s",
                        meth, conditionMessage(e)), call. = FALSE)
        NULL
      }
    )
    if (is.null(res)) {
      res <- .sim_cate_failure_rows(
        method = meth,
        n = n,
        kx = kx,
        k = k,
        s = s,
        tuning = tuning,
        tuning_mode = tuning_mode,
        alpha = alpha,
        by_candidate = identical(tuning_mode, "grid")
      )
    } else {
      res <- .sim_cate_normalize_method_rows(res)
    }
    res$method <- meth
    res[, c("method", .sim_cate_row_names()), drop = FALSE]
  })

  do.call(rbind, rows)
}

# ── Monte Carlo grid ──────────────────────────────────────────────────────────

#' Run the CATE simulation grid
#'
#' Expands a grid of DGP designs, sample sizes, and covariate dimensions, then
#' runs `n_reps` independent Monte Carlo replications for each cell using
#' [sim_cate_one_rep()].  Each replicate draws one dataset and evaluates **all**
#' requested `methods` on the same draw, enabling like-for-like comparisons.
#' Returns a tidy data.frame of raw per-replicate, per-method results.  Use
#' [summarise_sim_results()] to aggregate to one row per
#' `(method, design, n, kx, candidate_id, lambda, s)`.
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
#' [summarise_sim_results()] reports `n_failed` per candidate-aware group.
#'
#' @param designs Character vector of DGP names.
#' @param n_vals Integer vector of sample sizes.
#' @param kx_vals Integer vector of covariate dimensions.
#' @param n_reps Number of Monte Carlo replications per cell.
#' @param master_seed Integer seed used to draw the per-replicate seeds.
#'   Every cell uses the same seed vector so noise realisations are comparable
#'   across DGPs, sample sizes, and methods (fair comparison).
#' @param k Integer DNN neighbour count. The paper-facing DNN default is
#'   `k = 1L`, meaning one nearest neighbor is selected inside each size-`s`
#'   subsample; larger values use the package's generalized top-`k` weights.
#'   The fixed-`k` extension is essentially direct, but the current paper
#'   derivations are written for `k = 1`.
#' @param s_fn Function `f(n, kx)` returning the integer kernel order.
#'   Default satisfies the undersmoothing condition
#'   \eqn{s \asymp n^{(k_x+1)/(k_x+2)}}.
#' @param tuning Optional kernel-order tuning specification. When supplied,
#'   `s` is resolved per simulation cell from `tuning_eval(..., mode = "default")`
#'   and `s_fn` must be omitted.
#' @param tuning_mode One of `"default"` or `"grid"`. `"default"` preserves
#'   scalar simulation behavior. `"grid"` requires `tuning`; DNN methods emit
#'   one row per tuning candidate from the integrated [cate_dnn()] grid payload,
#'   while non-DNN benchmarks emit one row with `candidate_id = NA` and `s = NA`.
#' @param methods Character vector of method names; see [sim_cate_methods()].
#'   Default: `"dnn_oracle"`.
#' @param alpha Nominal coverage level.
#' @param sigma Noise SD forwarded to the DGP.
#' @param clip Propensity clipping bound forwarded to the DGP and cross-fitting.
#' @param n_folds Number of cross-fitting folds for the `dnn_crossfit` method.
#' @param verbose Print one progress line per cell (default `TRUE`).
#'
#' @return A `data.frame` with raw rows and columns `method`, `design`, `n`,
#'   `kx`, `k`, `tuning_source`, `candidate_id`, `lambda`, `s`, `seed`,
#'   `tau_hat`, `tau0`, `bias`, `sq_err`, `se`, `cover`, `ci_len`, `neff`,
#'   localization and studentization diagnostics described in
#'   [sim_cate_one_rep()], and `status`. Failed entries have `NA` for all
#'   metric columns.
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
#'   methods = c("dnn_fullsample", "dnn_crossfit", "ols_interact", "causal_forest")
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
    k           = 1L,
    s_fn        = function(n, kx) floor(n^((kx + 1L) / (kx + 2L))),
    tuning      = NULL,
    tuning_mode = c("default", "grid"),
    methods     = "dnn_oracle",
    alpha       = 0.05,
    sigma       = 0.5,
    clip        = 1e-3,
    n_folds     = 5L,
    verbose     = TRUE
) {
  tuning_mode <- match.arg(tuning_mode)
  s_fn_arg <- if (missing(s_fn)) NULL else s_fn
  if (!is.null(tuning) && !is.null(s_fn_arg)) {
    .abort_bad_arg(
      "tuning/s_fn",
      "cannot provide both `tuning` and `s_fn`; omit `s_fn` when `tuning` is used",
      quoted = FALSE
    )
  }
  if (identical(tuning_mode, "grid") && is.null(tuning)) {
    .abort_bad_arg(
      "tuning_mode",
      "requires `tuning` when `tuning_mode = \"grid\"`",
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
    cell <- .sim_cate_cell_config(
      tuning = tuning,
      tuning_mode = tuning_mode,
      s_fn = s_fn,
      n = n,
      kx = kx,
      k = k
    )

    if (verbose) {
      if (identical(tuning_mode, "grid")) {
        skipped <- sum(!cell$candidates$feasible)
        message(sprintf(
          "  design=%-22s  n=%4d  kx=%d  s in {%s}  skipped=%d  methods=[%s]",
          d, n, kx, paste(sort(unique(cell$candidates$s)), collapse = ", "),
          skipped, paste(methods, collapse = ", ")
        ))
      } else {
        message(sprintf("  design=%-22s  n=%4d  kx=%d  s=%d  methods=[%s]",
                        d, n, kx, cell$s, paste(methods, collapse = ", ")))
      }
    }

    reps <- .lapply(as.list(seeds), function(sd) {
      result <- tryCatch(
        sim_cate_one_rep(n = n, kx = kx, k = k, s = cell$s,
                         design = d, methods = methods,
                         tuning = tuning, tuning_mode = tuning_mode,
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
        # Total rep failure: shape rows by the rows each method would have
        # emitted in this cell.
        do.call(rbind, lapply(methods, function(meth) {
          fail <- .sim_cate_failure_rows(
            method = meth,
            n = n,
            kx = kx,
            k = k,
            s = cell$s,
            tuning = tuning,
            tuning_mode = tuning_mode,
            alpha = alpha,
            by_candidate = TRUE
          )
          .sim_cate_add_cell_columns(
            fail,
            method = meth,
            d = d,
            n = n,
            kx = kx,
            k = k,
            tuning_source = cell$tuning_source,
            seed = sd
          )
        }))
      } else {
        do.call(rbind, lapply(split(result, seq_len(nrow(result))), function(row) {
          .sim_cate_add_cell_columns(
            row,
            method = row$method[[1L]],
            d = d,
            n = n,
            kx = kx,
            k = k,
            tuning_source = cell$tuning_source,
            seed = sd
          )
        }))
      }
    })

    all_rows[[i]] <- do.call(rbind, reps)
  }

  do.call(rbind, all_rows)
}
