# sim-iv-run.R
# Single-replicate estimator and Monte Carlo harness for the IV-GMM simulation.

.sim_iv_metric_names <- function() {
  c(
    "theta_hat", "theta0", "bias", "sq_err", "se", "cover", "ci_len",
    "neff", "J", "df", "pval_J", "reject_J"
  )
}

.sim_iv_na_row <- function(theta0 = NA_real_, alpha = 0.05) {
  .sim_metric_iv(
    theta_hat = NA_real_,
    se_hat = NA_real_,
    theta0 = theta0,
    alpha = alpha,
    neff = NA_real_,
    J = NA_real_,
    df = NA_integer_
  )
}

#' Run one replicate of IV-GMM estimators
#'
#' Draws one dataset from a named IV DGP via [sim_iv_data()], then evaluates all
#' requested methods on that same draw.
#'
#' @param n Integer sample size.
#' @param kx Integer covariate dimension.
#' @param q Integer number of instruments requested from the DGP.
#' @param k,s Integer DNN tuning parameters.
#' @param design IV DGP name; passed to [sim_iv_data()].
#' @param methods Character vector of method names; see [sim_iv_methods()].
#' @param alpha Nominal test size for confidence intervals and J-tests.
#' @param seed Optional integer seed.
#' @param sigma Noise SD forwarded to [sim_iv_data()].
#' @param clip Clipping bound forwarded to DGPs and simulation nuisance fits.
#' @param n_folds Number of cross-fitting folds for cross-fit methods.
#'
#' @return A `data.frame` with one row per method and columns
#'   `method`, `theta_hat`, `theta0`, `bias`, `sq_err`, `se`, `cover`,
#'   `ci_len`, `neff`, `J`, `df`, `pval_J`, and `reject_J`.
#' @examples
#' sim_iv_one_rep(
#'   n = 120, kx = 2, q = 2, k = 1, s = 10,
#'   design = "dgp2_iv_homog_overid",
#'   methods = c("dnn_iv_oracle", "iv_oracle"),
#'   seed = 1
#' )
#' @export
sim_iv_one_rep <- function(n, kx, q, k, s, design,
                           methods = "dnn_iv_oracle",
                           alpha = 0.05,
                           seed = NULL,
                           sigma = 0.5,
                           clip = 1e-3,
                           n_folds = 5L) {
  dgp <- sim_iv_data(
    design = design,
    n = n,
    kx = kx,
    q = q,
    sigma = sigma,
    clip = clip,
    seed = seed
  )
  x_cols <- paste0("X", seq_len(kx))
  metric_names <- .sim_iv_metric_names()
  na_row <- .sim_iv_na_row(theta0 = dgp$theta_x0, alpha = alpha)

  rows <- lapply(methods, function(meth) {
    fn <- .get_sim_iv_method(meth)
    res <- tryCatch(
      fn(
        dgp = dgp,
        x_cols = x_cols,
        k = k,
        s = s,
        alpha = alpha,
        n_folds = n_folds,
        clip = clip
      ),
      error = function(e) {
        warning(sprintf(
          "sim_iv_one_rep: method '%s' failed: %s",
          meth, conditionMessage(e)
        ), call. = FALSE)
        NULL
      }
    )
    if (is.null(res)) {
      na_row
    } else {
      missing <- setdiff(metric_names, names(res))
      if (length(missing) > 0L) {
        res <- c(res, stats::setNames(rep(NA_real_, length(missing)), missing))
      }
      res[metric_names]
    }
  })

  do.call(rbind, lapply(seq_along(methods), function(i) {
    r <- rows[[i]]
    data.frame(
      method = methods[i],
      theta_hat = r[["theta_hat"]],
      theta0 = r[["theta0"]],
      bias = r[["bias"]],
      sq_err = r[["sq_err"]],
      se = r[["se"]],
      cover = r[["cover"]],
      ci_len = r[["ci_len"]],
      neff = r[["neff"]],
      J = r[["J"]],
      df = r[["df"]],
      pval_J = r[["pval_J"]],
      reject_J = r[["reject_J"]],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))
}

.sim_iv_cells <- function(designs, n_vals, kx_vals, q_vals) {
  rows <- list()
  idx <- 1L
  for (design in designs) {
    q_for_design <- if (identical(design, "dgp1_iv_justid_linear")) 1L else q_vals
    for (n in n_vals) {
      for (kx in kx_vals) {
        for (q in q_for_design) {
          rows[[idx]] <- data.frame(
            design = design,
            n = n,
            kx = kx,
            q = q,
            stringsAsFactors = FALSE
          )
          idx <- idx + 1L
        }
      }
    }
  }
  do.call(rbind, rows)
}

#' Run the IV-GMM simulation grid
#'
#' Expands a grid over DGP designs, sample sizes, covariate dimensions, and
#' instrument counts, then runs [sim_iv_one_rep()] for each replicate.  Each
#' cell uses the same replicate seed vector so methods are compared on aligned
#' draws.
#'
#' @param designs Character vector of IV DGP names.
#' @param n_vals Integer vector of sample sizes.
#' @param kx_vals Integer vector of covariate dimensions.
#' @param q_vals Integer vector of requested instrument counts.  The
#'   just-identified design is run once with realized `q = 1`.
#' @param n_reps Number of Monte Carlo replications per cell.
#' @param master_seed Seed used to generate replicate seeds.
#' @param k Integer DNN neighbor count.
#' @param s_fn Function `f(n, kx)` returning kernel order. The finite-sample
#'   default is \eqn{\lfloor 0.4 n^{(k_x+1)/(k_x+2)} \rfloor}.
#' @param tuning Optional tuning specification.  When supplied, omit `s_fn`.
#' @param methods Character vector of methods; see [sim_iv_methods()].
#' @param alpha Nominal test size.
#' @param sigma Noise SD forwarded to [sim_iv_data()].
#' @param clip Clipping bound forwarded to DGPs and nuisance fits.
#' @param n_folds Number of cross-fitting folds.
#' @param verbose Print one progress line per cell.
#'
#' @return Raw replicate-level `data.frame`.
#' @seealso [sim_iv_one_rep()], [summarise_iv_results()], [sim_iv_methods()]
#' @examples
#' \dontrun{
#' raw <- sim_iv_run(
#'   designs = "dgp2_iv_homog_overid",
#'   n_vals = 200,
#'   q_vals = 2,
#'   n_reps = 10,
#'   methods = c("dnn_iv_oracle", "iv_oracle")
#' )
#' summarise_iv_results(raw)
#' }
#' @export
sim_iv_run <- function(
    designs = c("dgp1_iv_justid_linear", "dgp2_iv_homog_overid",
                "dgp3_iv_heterog_overid", "dgp4_iv_weak_local",
                "dgp5_iv_correlated_instruments",
                "dgp6_iv_nonlinear_firststage"),
    n_vals = c(500L, 1000L, 2000L),
    kx_vals = 2L,
    q_vals = c(1L, 2L, 3L),
    n_reps = 500L,
    master_seed = 42L,
    k = 1L,
    s_fn = function(n, kx) floor(0.4 * n^((kx + 1L) / (kx + 2L))),
    tuning = NULL,
    methods = "dnn_iv_oracle",
    alpha = 0.05,
    sigma = 0.5,
    clip = 1e-3,
    n_folds = 5L,
    verbose = TRUE
) {
  n_reps <- .validate_scalar_int(n_reps, arg = "n_reps", min = 1L)
  k <- .validate_scalar_int(k, arg = "k", min = 1L)
  q_vals <- vapply(q_vals, .validate_scalar_int, integer(1L), arg = "q_vals", min = 1L)
  if (any(q_vals > 3L)) {
    .abort_bad_arg("q_vals", "must be <= 3 for the Phase-A IV simulation DGPs")
  }

  s_fn_arg <- if (missing(s_fn)) NULL else s_fn
  if (!is.null(tuning) && !is.null(s_fn_arg)) {
    .abort_bad_arg(
      "tuning/s_fn",
      "cannot provide both `tuning` and `s_fn`; omit `s_fn` when `tuning` is used",
      quoted = FALSE
    )
  }
  if (is.null(tuning) && is.null(s_fn_arg)) {
    s_fn <- function(n, kx) floor(0.4 * n^((kx + 1L) / (kx + 2L)))
  } else if (is.null(tuning)) {
    s_fn <- s_fn_arg
  }

  .lapply <- if (requireNamespace("future.apply", quietly = TRUE)) {
    function(X, FUN) future.apply::future_lapply(X, FUN, future.seed = TRUE)
  } else {
    lapply
  }

  set.seed(master_seed)
  seeds <- sample.int(1e6L, n_reps)
  cells <- .sim_iv_cells(
    designs = designs,
    n_vals = n_vals,
    kx_vals = kx_vals,
    q_vals = unique(q_vals)
  )
  metric_names <- .sim_iv_metric_names()
  all_rows <- vector("list", nrow(cells))

  for (i in seq_len(nrow(cells))) {
    d <- cells$design[i]
    n <- cells$n[i]
    kx <- cells$kx[i]
    q <- cells$q[i]
    s <- if (!is.null(tuning)) {
      tuning_eval(
        tuning = tuning,
        n = n,
        kx = kx,
        k = k,
        se = FALSE,
        jackknife = "none",
        mode = "default"
      )$s
    } else {
      s_fn(n, kx)
    }

    if (verbose) {
      message(sprintf(
        "  design=%-32s n=%4d kx=%d q=%d s=%d methods=[%s]",
        d, n, kx, q, s, paste(methods, collapse = ", ")
      ))
    }

    reps <- .lapply(as.list(seeds), function(sd) {
      result <- tryCatch(
        sim_iv_one_rep(
          n = n,
          kx = kx,
          q = q,
          k = k,
          s = s,
          design = d,
          methods = methods,
          alpha = alpha,
          sigma = sigma,
          clip = clip,
          seed = sd,
          n_folds = n_folds
        ),
        error = function(e) {
          warning(sprintf(
            "sim_iv_run: rep seed %d failed [design=%s n=%d kx=%d q=%d]: %s",
            sd, d, n, kx, q, conditionMessage(e)
          ), call. = FALSE)
          NULL
        }
      )

      if (is.null(result)) {
        do.call(rbind, lapply(methods, function(meth) {
          row <- as.list(.sim_iv_na_row(alpha = alpha)[metric_names])
          data.frame(
            design = d, n = n, kx = kx, q = q, s = s, seed = sd,
            method = meth,
            as.data.frame(row, check.names = FALSE),
            stringsAsFactors = FALSE,
            row.names = NULL
          )
        }))
      } else {
        result$design <- d
        result$n <- n
        result$kx <- kx
        result$q <- q
        result$s <- s
        result$seed <- sd
        result[, c("design", "n", "kx", "q", "s", "seed", "method", metric_names),
               drop = FALSE]
      }
    })

    all_rows[[i]] <- do.call(rbind, reps)
  }

  do.call(rbind, all_rows)
}
