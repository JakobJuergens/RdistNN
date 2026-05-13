# sim-iv-summary.R
# Aggregation utility for raw IV-GMM simulation results.

.sim_iv_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

.sim_iv_median <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  stats::median(x, na.rm = TRUE)
}

#' Summarise raw IV-GMM simulation results
#'
#' Aggregates the raw per-replicate output from [sim_iv_run()] to one row per
#' `(method, design, n, kx, q)` cell.
#'
#' @param raw A data frame with columns returned by [sim_iv_run()].
#'
#' @return A data frame with bias, RMSE, mean standard error, coverage, mean CI
#'   length, effective sample size, and overidentification diagnostics.
#' @seealso [sim_iv_run()]
#' @examples
#' raw <- sim_iv_run(
#'   designs = "dgp1_iv_justid_linear",
#'   n_vals = 100,
#'   q_vals = c(1, 2),
#'   n_reps = 2,
#'   methods = "iv_oracle",
#'   verbose = FALSE
#' )
#' summarise_iv_results(raw)
#' @export
summarise_iv_results <- function(raw) {
  required <- c(
    "method", "design", "n", "kx", "q", "bias", "sq_err", "se",
    "cover", "ci_len", "neff", "J", "pval_J", "reject_J"
  )
  missing <- setdiff(required, names(raw))
  if (length(missing) > 0L) {
    stop("raw is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  group_cols <- c("method", "design", "n", "kx", "q")
  cells <- unique(raw[, group_cols, drop = FALSE])
  cells <- cells[do.call(order, unname(as.list(cells[, group_cols, drop = FALSE]))),
                 , drop = FALSE]
  rownames(cells) <- NULL

  out <- vector("list", nrow(cells))
  for (i in seq_len(nrow(cells))) {
    idx <- rep(TRUE, nrow(raw))
    for (col in group_cols) {
      idx <- idx & (raw[[col]] == cells[[col]][i])
    }
    sub <- raw[idx, , drop = FALSE]
    n_failed <- sum(is.na(sub$bias))

    row <- data.frame(
      method = cells$method[i],
      design = cells$design[i],
      n = cells$n[i],
      kx = cells$kx[i],
      q = cells$q[i],
      n_reps = nrow(sub),
      n_failed = n_failed,
      bias = .sim_iv_mean(sub$bias),
      rmse = if (all(is.na(sub$sq_err))) NA_real_ else sqrt(mean(sub$sq_err, na.rm = TRUE)),
      se_mean = .sim_iv_mean(sub$se),
      cover = .sim_iv_mean(sub$cover),
      ci_len = .sim_iv_mean(sub$ci_len),
      neff = .sim_iv_mean(sub$neff),
      J_median = .sim_iv_median(sub$J),
      pval_J_median = .sim_iv_median(sub$pval_J),
      reject_J = .sim_iv_mean(sub$reject_J),
      stringsAsFactors = FALSE
    )
    out[[i]] <- row
  }

  result <- do.call(rbind, out)
  rownames(result) <- NULL
  result
}
