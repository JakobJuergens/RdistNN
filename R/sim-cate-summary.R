# sim-cate-summary.R
# Aggregation utility for raw per-replicate CATE simulation results.
#
# Public API: summarise_sim_results()

#' Summarise raw per-replicate CATE simulation results
#'
#' Aggregates a data.frame of raw replicate-level results (as returned by
#' [sim_cate_run()]) to one summary row per `(design, n, kx)` cell.
#' Reports bias, RMSE, mean SE, empirical coverage, mean CI length, and
#' (optionally) mean effective sample size.  Failed replicates (rows where
#' metric columns are `NA`) are excluded from all averages; their count is
#' reported in `n_failed`.
#'
#' @param raw A `data.frame` with at least the columns `design`, `n`, `kx`,
#'   `bias`, `sq_err`, `se`, `cover`, `ci_len`.  An optional `neff` column
#'   is included in the output when present.
#'
#' @return A `data.frame` with one row per unique `(design, n, kx)` cell and
#'   columns `design`, `n`, `kx`, `n_reps`, `n_failed`, `bias`, `rmse`,
#'   `se_mean`, `cover`, `ci_len` (and `neff` when present in `raw`).
#'
#' @seealso [sim_cate_run()]
#'
#' @examples
#' \dontrun{
#' raw <- sim_cate_run(designs = "dgp1_linear", n_vals = 500L, n_reps = 10L)
#' summarise_sim_results(raw)
#' }
#'
#' @export
summarise_sim_results <- function(raw) {
  required <- c("design", "n", "kx", "bias", "sq_err", "se", "cover", "ci_len")
  missing  <- setdiff(required, names(raw))
  if (length(missing) > 0L)
    stop("raw is missing required columns: ", paste(missing, collapse = ", "))

  has_neff <- "neff" %in% names(raw)

  cells <- unique(raw[, c("design", "n", "kx"), drop = FALSE])
  cells <- cells[order(cells$design, cells$n, cells$kx), , drop = FALSE]
  rownames(cells) <- NULL

  out <- vector("list", nrow(cells))

  for (i in seq_len(nrow(cells))) {
    idx <- raw$design == cells$design[i] &
           raw$n      == cells$n[i]      &
           raw$kx     == cells$kx[i]
    sub <- raw[idx, , drop = FALSE]

    n_failed <- sum(is.na(sub$bias))

    row <- data.frame(
      design   = cells$design[i],
      n        = cells$n[i],
      kx       = cells$kx[i],
      n_reps   = nrow(sub),
      n_failed = n_failed,
      bias     = mean(sub$bias,    na.rm = TRUE),
      rmse     = sqrt(mean(sub$sq_err, na.rm = TRUE)),
      se_mean  = mean(sub$se,      na.rm = TRUE),
      cover    = mean(sub$cover,   na.rm = TRUE),
      ci_len   = mean(sub$ci_len,  na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    if (has_neff) row$neff <- mean(sub$neff, na.rm = TRUE)
    out[[i]] <- row
  }

  do.call(rbind, out)
}
