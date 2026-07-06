# sim-cate-summary.R
# Aggregation utility for raw per-replicate CATE simulation results.
#
# Public API: summarise_sim_results()

#' Summarise raw per-replicate CATE simulation results
#'
#' Aggregates a data.frame of raw replicate-level results (as returned by
#' [sim_cate_run()]) to one summary row per
#' `(method, design, n, kx, candidate_id, kappa, s)` cell. Reports bias, RMSE,
#' mean SE, empirical coverage, mean CI length, and (optionally) mean effective
#' sample size. Failed or skipped rows are excluded from all averages; their
#' count is reported in `n_failed`.
#'
#' @param raw A `data.frame` with at least the columns `method`, `design`, `n`,
#'   `kx`, `candidate_id`, `kappa`, `s`, `bias`, `sq_err`, `se`, `cover`,
#'   `ci_len`, and `status`. Optional columns `neff` and `tau_hat`/`tau0`
#'   (returned by [sim_cate_run()]) are handled gracefully: `neff` is averaged
#'   into the summary; `tau_hat` and `tau0` are per-rep quantities and are not
#'   aggregated (they live only in the raw data).
#'
#' @return A `data.frame` with one row per unique group and columns
#'   `method`, `design`, `n`, `kx`, `candidate_id`, `kappa`, `s`, `n_reps`,
#'   `n_failed`, `bias`, `rmse`, `se_mean`, `cover`, `ci_len` (and `neff`
#'   when present in `raw`).
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

  has_neff   <- "neff"   %in% names(raw)
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

    row <- data.frame(
      design   = cells$design[i],
      n        = cells$n[i],
      kx       = cells$kx[i],
      candidate_id = cells$candidate_id[i],
      kappa    = cells$kappa[i],
      s        = cells$s[i],
      n_reps   = nrow(sub),
      n_failed = n_failed,
      bias     = mean(sub$bias,    na.rm = TRUE),
      rmse     = sqrt(mean(sub$sq_err, na.rm = TRUE)),
      se_mean  = mean(sub$se,      na.rm = TRUE),
      cover    = mean(sub$cover,   na.rm = TRUE),
      ci_len   = mean(sub$ci_len,  na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    row$method <- cells$method[i]
    if (has_neff)   row$neff   <- mean(sub$neff, na.rm = TRUE)
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
