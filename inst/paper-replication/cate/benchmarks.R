# inst/paper-replication/cate/benchmarks.R
#
# Competing estimators for the CATE simulation study.
# Reference: 05b_CATE_Simulation.tex §Estimators
#
# Each estimator follows the same interface as sim_cate_one_rep():
#
#   one_rep_XXX(n, kx, design, alpha, seed, sigma, clip)
#     -> named numeric vector: bias, sq_err, se, cover, ci_len, with GRF
#        nuisance-route metadata attached for forest estimators
#        (neff is not meaningful for non-DNN methods)
#
# All competitor packages are listed under Suggests in DESCRIPTION and loaded
# with requireNamespace(..., quietly = TRUE) to avoid hard dependencies.

library(RdistNN)

# ── Helper: check and attach a suggested package ─────────────────────────────

.bench_require <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(pkg, " is required for this benchmark. Install with: ",
         'install.packages("', pkg, '")')
  invisible(TRUE)
}

# ── Shared: draw data and extract standard fields ────────────────────────────

.bench_draw <- function(n, kx, design, sigma, clip, seed) {
  dgp    <- sim_cate_data(design, n = n, kx = kx,
                          sigma = sigma, clip = clip, seed = seed)
  x_cols <- paste0("X", seq_len(kx))
  list(
    dat    = dgp$data,
    x0     = dgp$x0,
    tau0   = dgp$tau_x0,
    X      = as.matrix(dgp$data[, x_cols]),
    D      = dgp$data$D,
    Y      = dgp$data$Y,
    e      = dgp$data$e_oracle,
    m0     = dgp$data$m0_oracle,
    m1     = dgp$data$m1_oracle
  )
}

.bench_result <- function(tau_hat, se_hat, tau0, alpha) {
  z       <- stats::qnorm(1 - alpha / 2)
  covered <- (tau_hat - z * se_hat <= tau0) & (tau0 <= tau_hat + z * se_hat)
  c(bias   = tau_hat - tau0,
    sq_err = (tau_hat - tau0)^2,
    se     = se_hat,
    cover  = as.numeric(covered),
    ci_len = 2 * z * se_hat)
}

# ── Modern GRF causal forest with internally estimated nuisances ─────────────

#' @param ... Forwarded to grf::causal_forest
one_rep_causal_forest <- function(n, kx, design, alpha = 0.05, seed = NULL,
                                  sigma = 0.5, clip = 1e-3, ...) {
  .bench_require("grf")
  d <- .bench_draw(n, kx, design, sigma, clip, seed)
  if (!is.null(seed)) set.seed(seed)

  cf  <- grf::causal_forest(X = d$X, Y = d$Y, W = d$D, ...)
  pred <- stats::predict(cf, newdata = matrix(d$x0, nrow = 1L),
                         estimate.variance = TRUE)
  tau_hat <- pred$predictions[1L]
  se_hat  <- sqrt(pred$variance.estimates[1L])

  out <- .bench_result(tau_hat, se_hat, d$tau0, alpha)
  attr(out, "method_variant") <- "grf_internal_nuisance"
  attr(out, "forest_nuisance") <- "internal"
  attr(out, "y_hat_supplied") <- FALSE
  attr(out, "w_hat_supplied") <- FALSE
  out
}

# ── Modern: local linear forest (Friedberg et al.) ───────────────────────────

one_rep_llf <- function(n, kx, design, alpha = 0.05, seed = NULL,
                        sigma = 0.5, clip = 1e-3, ...) {
  .bench_require("grf")
  d <- .bench_draw(n, kx, design, sigma, clip, seed)
  if (!is.null(seed)) set.seed(seed)

  cf  <- grf::causal_forest(X = d$X, Y = d$Y, W = d$D,
                             ll.lambda = 0.1, ll.weight.penalty = FALSE, ...)
  pred <- stats::predict(cf, newdata = matrix(d$x0, nrow = 1L),
                         ll.lambda = 0.1, estimate.variance = TRUE,
                         linear.correction.variables = seq_len(kx))
  tau_hat <- pred$predictions[1L]
  se_hat  <- sqrt(pred$variance.estimates[1L])

  out <- .bench_result(tau_hat, se_hat, d$tau0, alpha)
  attr(out, "method_variant") <- "grf_local_linear_internal_nuisance"
  attr(out, "forest_nuisance") <- "internal"
  attr(out, "y_hat_supplied") <- FALSE
  attr(out, "w_hat_supplied") <- FALSE
  out
}

# ── Classical: naive OLS (Y ~ T + X) ─────────────────────────────────────────

one_rep_ols_naive <- function(n, kx, design, alpha = 0.05, seed = NULL,
                               sigma = 0.5, clip = 1e-3) {
  d   <- .bench_draw(n, kx, design, sigma, clip, seed)
  dat <- d$dat

  fmla    <- stats::as.formula(
    paste("Y ~ D +", paste(paste0("X", seq_len(kx)), collapse = " + "))
  )
  fit     <- stats::lm(fmla, data = dat)
  coefs   <- stats::coef(summary(fit))
  tau_hat <- coefs["D", "Estimate"]
  se_hat  <- coefs["D", "Std. Error"]

  .bench_result(tau_hat, se_hat, d$tau0, alpha)
}

# ── Classical: interaction OLS (Y ~ T + X + T*X) ─────────────────────────────

one_rep_ols_interact <- function(n, kx, design, alpha = 0.05, seed = NULL,
                                  sigma = 0.5, clip = 1e-3) {
  d   <- .bench_draw(n, kx, design, sigma, clip, seed)
  dat <- d$dat
  x0  <- d$x0

  x_terms <- paste0("X", seq_len(kx))
  fmla    <- stats::as.formula(paste(
    "Y ~ D * (", paste(x_terms, collapse = " + "), ")"
  ))
  fit <- stats::lm(fmla, data = dat)

  # tau(x0) = b_D + b_{D:X1}*x0[1] + ... + b_{D:Xkx}*x0[kx]
  b     <- stats::coef(fit)
  tau_hat <- b["D"] + sum(b[paste0("D:X", seq_len(kx))] * x0, na.rm = TRUE)

  # Delta method SE (approximate)
  V     <- stats::vcov(fit)
  idx   <- c("D", paste0("D:X", seq_len(kx)))
  idx   <- intersect(idx, rownames(V))
  grad  <- c(1, x0[seq_along(idx) - 1L])
  se_hat <- sqrt(as.numeric(t(grad) %*% V[idx, idx] %*% grad))

  .bench_result(tau_hat, se_hat, d$tau0, alpha)
}

# ── TODO: additional benchmarks ───────────────────────────────────────────────
#
# one_rep_dml_kernel()   -- DML + local-linear kernel-smoothed AIPW
# one_rep_dr_learner()   -- DR-learner (Kennedy); pointwise eval at x0
# one_rep_knn_matching() -- k-NN matching (Abadie-Imbens), k in {1, 5}
# one_rep_ps_matching()  -- propensity-score nearest-neighbour matching
# one_rep_ipw()          -- Hajek-normalised IPW (global, not localised)

# ── Grid runner ───────────────────────────────────────────────────────────────

#' Run all implemented benchmark estimators over the paper grid
#'
#' @param grid A list with the same fields as PAPER_GRID in run.R.
#' @return A data.frame with columns method, design, n, kx, bias, rmse,
#'   se_mean, cover, ci_len (one row per estimator x cell).
run_benchmarks <- function(grid) {
  set.seed(grid$master_seed)
  seeds <- sample.int(1e6L, grid$n_reps)

  estimators <- list(
    causal_forest = one_rep_causal_forest,
    llf           = one_rep_llf,
    ols_naive     = one_rep_ols_naive,
    ols_interact  = one_rep_ols_interact
  )

  cells <- expand.grid(
    design = grid$designs,
    n      = grid$n_vals,
    kx     = grid$kx_vals,
    stringsAsFactors = FALSE
  )

  out <- vector("list", nrow(cells) * length(estimators))
  idx <- 1L

  for (i in seq_len(nrow(cells))) {
    d  <- cells$design[i]
    n  <- cells$n[i]
    kx <- cells$kx[i]

    for (nm in names(estimators)) {
      fn <- estimators[[nm]]
      message(sprintf("  %-20s  design=%-22s  n=%4d  kx=%d", nm, d, n, kx))

      reps <- lapply(seeds, function(sd) {
        tryCatch(
          fn(n = n, kx = kx, design = d, alpha = grid$alpha,
             sigma = grid$sigma, clip = grid$clip, seed = sd),
          error = function(e) rep(NA_real_, 5L)
        )
      })
      reps <- do.call(rbind, reps)

      out[[idx]] <- data.frame(
        method  = nm,
        method_label = if (nm == "causal_forest") {
          "Causal forest (GRF, internal nuisances)"
        } else if (nm == "llf") {
          "Local linear forest (GRF, internal nuisances)"
        } else {
          nm
        },
        method_variant = if (nm == "causal_forest") {
          "grf_internal_nuisance"
        } else if (nm == "llf") {
          "grf_local_linear_internal_nuisance"
        } else {
          NA_character_
        },
        forest_nuisance = if (nm %in% c("causal_forest", "llf")) "internal" else NA_character_,
        y_hat_supplied = if (nm %in% c("causal_forest", "llf")) FALSE else NA,
        w_hat_supplied = if (nm %in% c("causal_forest", "llf")) FALSE else NA,
        design  = d,
        n       = n,
        kx      = kx,
        bias    = mean(reps[, "bias"],   na.rm = TRUE),
        rmse    = sqrt(mean(reps[, "sq_err"], na.rm = TRUE)),
        se_mean = mean(reps[, "se"],     na.rm = TRUE),
        cover   = mean(reps[, "cover"],  na.rm = TRUE),
        ci_len  = mean(reps[, "ci_len"], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  do.call(rbind, out)
}
