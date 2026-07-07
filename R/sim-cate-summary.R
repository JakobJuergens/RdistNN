# sim-cate-summary.R
# Aggregation utility for raw per-replicate CATE simulation results.
#
# Public API: summarise_sim_results()

.sim_cate_mean <- function(x) {
  if (length(x) == 0L || all(is.na(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

.sim_cate_rmse <- function(x) {
  if (length(x) == 0L || all(is.na(x))) {
    return(NA_real_)
  }
  sqrt(mean(x, na.rm = TRUE))
}

.sim_cate_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) {
    return(NA_real_)
  }
  stats::sd(x)
}

.sim_cate_se_to_sd <- function(se_mean, emp_sd) {
  if (!is.finite(se_mean) || !is.finite(emp_sd) || se_mean <= 0 || emp_sd <= 0) {
    return(NA_real_)
  }
  se_mean / emp_sd
}

#' Summarise raw per-replicate CATE simulation results
#'
#' Aggregates a data.frame of raw replicate-level results (as returned by
#' [sim_cate_run()]) to one summary row per
#' `(method, design, n, kx, candidate_id, kappa, s)` cell. Reports bias, RMSE,
#' mean SE, empirical coverage, mean CI length, empirical sampling dispersion of
#' `tau_hat`, the ratio of mean SE to empirical dispersion, and (optionally)
#' mean effective sample size. Failed or skipped rows are excluded from all
#' averages and dispersion calculations; their count is reported in `n_failed`.
#'
#' @param raw A `data.frame` with at least the columns `method`, `design`, `n`,
#'   `kx`, `candidate_id`, `kappa`, `s`, `bias`, `sq_err`, `se`, `cover`,
#'   `ci_len`, and `status`. Optional columns `neff`, `tau_hat`, and `tau0`
#'   (returned by [sim_cate_run()]) are handled gracefully: `neff` is averaged
#'   into the summary; `tau_hat` is summarized through its empirical standard
#'   deviation; `tau0` remains a per-rep quantity and is not aggregated.
#'
#' @return A `data.frame` with one row per unique group and columns
#'   `method`, `design`, `n`, `kx`, `candidate_id`, `kappa`, `s`, `n_reps`,
#'   `n_failed`, `bias`, `rmse`, `se_mean`, `cover`, `ci_len`, `emp_sd`,
#'   `se_to_sd` (and `neff` when present in `raw`).
#'
#' @seealso [sim_cate_run()]
#'
#' @examples
#' \dontrun{
#' raw <- sim_cate_run(designs = "dgp1_linear", n_vals = 500L, n_reps = 10L,
#'                     methods = c("dnn_oracle", "ols_interact"))
#' summarise_sim_results(raw)
#' }
#'
#' @export
summarise_sim_results <- function(raw) {
  required <- c(
    "method", "design", "n", "kx", "candidate_id", "kappa", "s",
    "bias", "sq_err", "se", "cover", "ci_len", "status"
  )
  missing  <- setdiff(required, names(raw))
  if (length(missing) > 0L)
    stop("raw is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)

  has_neff    <- "neff" %in% names(raw)
  has_tau_hat <- "tau_hat" %in% names(raw)
  group_cols <- c("method", "design", "n", "kx", "candidate_id", "kappa", "s")

  cells <- unique(raw[, group_cols, drop = FALSE])
  cells <- cells[do.call(order, unname(as.list(cells[, group_cols, drop = FALSE]))),
                 , drop = FALSE]
  rownames(cells) <- NULL

  out <- vector("list", nrow(cells))

  for (i in seq_len(nrow(cells))) {
    # Row selector across all grouping columns
    idx <- rep(TRUE, nrow(raw))
    for (col in group_cols) {
      value <- cells[[col]][i]
      idx <- idx & if (is.na(value)) is.na(raw[[col]]) else raw[[col]] == value
    }
    sub <- raw[idx, , drop = FALSE]

    n_failed <- sum(is.na(sub$status) | sub$status != "ok" | is.na(sub$bias))
    ok <- !is.na(sub$status) & sub$status == "ok" & !is.na(sub$bias)

    se_mean <- .sim_cate_mean(sub$se[ok])
    emp_sd <- if (has_tau_hat) .sim_cate_sd(sub$tau_hat[ok]) else NA_real_

    row <- data.frame(
      design   = cells$design[i],
      n        = cells$n[i],
      kx       = cells$kx[i],
      candidate_id = cells$candidate_id[i],
      kappa    = cells$kappa[i],
      s        = cells$s[i],
      n_reps   = nrow(sub),
      n_failed = n_failed,
      bias     = .sim_cate_mean(sub$bias[ok]),
      rmse     = .sim_cate_rmse(sub$sq_err[ok]),
      se_mean  = se_mean,
      cover    = .sim_cate_mean(sub$cover[ok]),
      ci_len   = .sim_cate_mean(sub$ci_len[ok]),
      emp_sd   = emp_sd,
      se_to_sd = .sim_cate_se_to_sd(se_mean, emp_sd),
      stringsAsFactors = FALSE
    )
    row$method <- cells$method[i]
    if (has_neff)   row$neff   <- .sim_cate_mean(sub$neff[ok])
    out[[i]] <- row
  }

  result <- do.call(rbind, out)
  rownames(result) <- NULL

  # Put method first; neff last.
  id_cols     <- c("method", "design", "n", "kx", "candidate_id", "kappa", "s")
  other_cols  <- setdiff(names(result), id_cols)
  result      <- result[, c(id_cols, other_cols), drop = FALSE]
  result
}
