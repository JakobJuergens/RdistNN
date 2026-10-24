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

.sim_cate_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  as.numeric(stats::quantile(x, probs = probability, names = FALSE, type = 8))
}

.sim_cate_skewness <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3L) return(NA_real_)
  centered <- x - mean(x)
  scale <- sqrt(mean(centered^2))
  if (!is.finite(scale) || scale <= 0) NA_real_ else mean((centered / scale)^3)
}

.sim_cate_excess_kurtosis <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4L) return(NA_real_)
  centered <- x - mean(x)
  scale <- sqrt(mean(centered^2))
  if (!is.finite(scale) || scale <= 0) NA_real_ else mean((centered / scale)^4) - 3
}

#' Summarise raw per-replicate CATE simulation results
#'
#' Aggregates a data.frame of raw replicate-level results (as returned by
#' [sim_cate_run()]) to one summary row per
#' `(method, design, n, kx, candidate_id, lambda, s)` cell and, when present,
#' per GRF nuisance configuration. Reports bias, RMSE,
#' mean SE, empirical coverage, mean CI length, empirical sampling dispersion of
#' `tau_hat`, the ratio of mean SE to empirical dispersion, and (optionally)
#' mean effective sample size. When the raw rows carry localization diagnostics,
#' the summary also reports mean realized localization drift, the mean estimator
#' remainder, physical neighborhood radius, arm-specific effective sample
#' sizes, and shape diagnostics for the studentized statistic. Failed or
#' skipped rows are excluded from all averages and dispersion calculations;
#' their count is reported in `n_failed`. Here `localized_truth` is the realized
#' DNN-weighted truth, not the deterministic population branch. Under
#' complete-DNN exchangeability its expectation equals \eqn{\alpha_s(x_0)}, so
#' mean realized localization drift estimates the population branch
#' displacement. The estimator remainder contains centered score variation and
#' fitted-state effects.
#'
#' @param raw A `data.frame` with at least the columns `method`, `design`, `n`,
#'   `kx`, `candidate_id`, `lambda`, `s`, `bias`, `sq_err`, `se`, `cover`,
#'   `ci_len`, and `status`. Optional columns `neff`, `tau_hat`, and `tau0`
#'   (returned by [sim_cate_run()]) are handled gracefully: `neff` is averaged
#'   into the summary; `tau_hat` is summarized through its empirical standard
#'   deviation; `tau0` remains a per-rep quantity and is not aggregated.
#'   `method_variant`, `forest_nuisance`, `y_hat_supplied`, and
#'   `w_hat_supplied` extend the grouping key. These route columns are required;
#'   methods without GRF nuisances record typed missing values.
#'
#' @return A `data.frame` with one row per unique group and columns
#'   `method`, `design`, `n`, `kx`, `candidate_id`, `lambda`, `s`, `n_reps`,
#'   `n_failed`, `bias`, `rmse`, `se_mean`, `cover`, `ci_len`, `emp_sd`,
#'   `se_to_sd`, `bias_to_se`, localization, support, and studentized-shape
#'   diagnostics (and `neff` when present in `raw`). The output always includes
#'   the four GRF route fields, with typed missing values when they do not apply.
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
    "method", "design", "n", "kx", "candidate_id", "lambda", "s",
    "bias", "sq_err", "se", "cover", "ci_len", "status",
    "method_variant", "forest_nuisance", "y_hat_supplied", "w_hat_supplied"
  )
  missing  <- setdiff(required, names(raw))
  if (length(missing) > 0L)
    .dnn_abort("argument_invalid",
      paste0("raw is missing required columns: ", paste(missing, collapse = ", ")),
      argument = "raw"
    )

  route_cols <- c(
    "method_variant", "forest_nuisance", "y_hat_supplied", "w_hat_supplied"
  )
  raw$method_variant <- as.character(raw$method_variant)
  raw$forest_nuisance <- as.character(raw$forest_nuisance)
  raw$y_hat_supplied <- as.logical(raw$y_hat_supplied)
  raw$w_hat_supplied <- as.logical(raw$w_hat_supplied)

  has_neff    <- "neff" %in% names(raw)
  has_tau_hat <- "tau_hat" %in% names(raw)
  group_cols <- c(
    "method", "design", "n", "kx", "candidate_id", "lambda", "s",
    route_cols
  )

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
      lambda    = cells$lambda[i],
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
    diagnostic_mean <- function(name) {
      if (!name %in% names(sub)) NA_real_ else .sim_cate_mean(sub[[name]][ok])
    }
    diagnostic_quantile <- function(name, probability) {
      if (!name %in% names(sub)) NA_real_ else .sim_cate_quantile(sub[[name]][ok], probability)
    }
    studentized <- if ("studentized_error" %in% names(sub)) {
      sub$studentized_error[ok]
    } else {
      numeric()
    }
    row$bias_to_se <- if (is.finite(row$bias) && is.finite(se_mean) && se_mean > 0) {
      row$bias / se_mean
    } else {
      NA_real_
    }
    row$localized_truth_mean <- diagnostic_mean("localized_truth")
    row$localization_drift <- diagnostic_mean("localization_drift")
    row$estimator_remainder <- diagnostic_mean("estimator_remainder")
    row$score_sd0 <- diagnostic_mean("score_sd0")
    row$standardized_localization_drift <- diagnostic_mean("standardized_localization_drift")
    row$physical_radius_90_mean <- diagnostic_mean("physical_radius_90")
    row$physical_radius_90_p25 <- diagnostic_quantile("physical_radius_90", 0.25)
    row$physical_radius_90_p50 <- diagnostic_quantile("physical_radius_90", 0.50)
    row$physical_radius_90_p75 <- diagnostic_quantile("physical_radius_90", 0.75)
    row$treated_ess <- diagnostic_mean("treated_ess")
    row$control_ess <- diagnostic_mean("control_ess")
    row$top5_score_share <- diagnostic_mean("top5_score_share")
    row$studentized_mean <- .sim_cate_mean(studentized)
    row$studentized_sd <- .sim_cate_sd(studentized)
    row$studentized_skewness <- .sim_cate_skewness(studentized)
    row$studentized_excess_kurtosis <- .sim_cate_excess_kurtosis(studentized)
    row$studentized_q025 <- .sim_cate_quantile(studentized, 0.025)
    row$studentized_q05 <- .sim_cate_quantile(studentized, 0.05)
    row$studentized_q50 <- .sim_cate_quantile(studentized, 0.50)
    row$studentized_q95 <- .sim_cate_quantile(studentized, 0.95)
    row$studentized_q975 <- .sim_cate_quantile(studentized, 0.975)
    row$studentized_abs_q95 <- .sim_cate_quantile(abs(studentized), 0.95)
    row$method <- cells$method[i]
    for (col in route_cols) row[[col]] <- cells[[col]][i]
    if (has_neff)   row$neff   <- .sim_cate_mean(sub$neff[ok])
    out[[i]] <- row
  }

  result <- do.call(rbind, out)
  rownames(result) <- NULL

  # Put method first; neff last.
  id_cols <- c(
    "method", "design", "n", "kx", "candidate_id", "lambda", "s",
    route_cols
  )
  other_cols  <- setdiff(names(result), id_cols)
  result      <- result[, c(id_cols, other_cols), drop = FALSE]
  result
}
