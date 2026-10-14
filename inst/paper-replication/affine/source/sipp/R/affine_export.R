# Portable affine evidence and displays, exclusively from saved fits.
export_sipp_401k_affine <- function(fit, diagnostics, output_dir, profile_ids = NULL,
                                   allow_validation = FALSE) {
  if (!identical(fit$schema_version, 2L) ||
      !setequal(unique(fit$estimates$method), sipp_401k_estimation_methods("affine")))
    stop("Affine export requires schema2 and target-compatible DNN methods.", call. = FALSE)
  validation <- isTRUE(fit$run_info$smoke) || (!is.null(fit$run_info$sample_n) &&
    !is.na(fit$run_info$sample_n)) || fit$source$n_fit != fit$source$n_analysis
  if (validation && !allow_validation)
    stop("Smoke/subsample validation fits cannot be exported as accepted empirical evidence.", call. = FALSE)
  fingerprint <- sipp_401k_fit_fingerprint(fit)
  if (!identical(fingerprint, diagnostics$fit_fingerprint))
    stop("Diagnostic fingerprint does not match this fitted state.", call. = FALSE)
  frozen <- analyze_sipp_401k_frozen(fit)
  if (is.null(profile_ids)) {
    interior <- c("age40_inc_p25", "age40_inc_p50", "age40_inc_p75")
    profile_ids <- if (all(interior %in% fit$profiles$profile_id)) interior else fit$profiles$profile_id
  }
  if (!length(profile_ids) || anyDuplicated(profile_ids) ||
      any(!profile_ids %in% fit$profiles$profile_id))
    stop("Selected profile identifiers must exist and be distinct.", call. = FALSE)
  if (dir.exists(output_dir) || file.exists(output_dir))
    stop("Export directory already exists; choose a new directory.", call. = FALSE)
  if (!dir.create(output_dir, recursive = TRUE))
    stop("Cannot reserve export directory.", call. = FALSE)
  tables <- list(estimates_all = fit$estimates,
    estimates_selected = fit$estimates[fit$estimates$profile_id %in% profile_ids, ],
    support = diagnostics$support, information = frozen$information,
    covariance_routes = frozen$covariance_routes, paired_training = frozen$paired,
    component_means = frozen$component_means, component_covariance = frozen$component_covariance,
    tuning_all = frozen$tuning)
  for (name in names(tables)) write.csv(tables[[name]],
    file.path(output_dir, paste0(name, ".csv")), row.names = FALSE)
  for (scale in unique(fit$estimates$outcome_scale))
    render_sipp_401k_affine_displays(fit, frozen, output_dir, profile_ids, scale, validation)
  files <- list.files(output_dir, full.names = TRUE)
  names(files) <- basename(files)
  selection <- list(schema_version = 2L, purpose = if (validation) "validation" else "empirical",
    run_id = fit$run_info$run_id, fit_fingerprint = fingerprint,
    selected_profiles = profile_ids, all_profiles = fit$profiles$profile_id,
    method_metadata = fit$method_metadata, specification = fit$specification,
    inference = frozen$scope, generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    files = lapply(files, sipp_401k_file_hash))
  jsonlite::write_json(selection, file.path(output_dir, "selection.json"), pretty = TRUE,
                       auto_unbox = TRUE, null = "null", na = "null", digits = NA)
  invisible(selection)
}

render_sipp_401k_affine_displays <- function(fit, frozen, output_dir, profiles, scale, validation) {
  methods <- sipp_401k_estimation_methods()
  colors <- c("#a65d28", "#17688c", "#7f6188")
  labels <- c("All-tree", "Single fit/OOB", paste0(fit$config$n_folds, "-fold"))
  p <- fit$profiles[match(profiles, fit$profiles$profile_id), ]
  incomes <- paste0("$", format(round(p$inc), big.mark = ",", scientific = FALSE, trim = TRUE))
  rows <- fit$estimates[fit$estimates$outcome_scale == scale & fit$estimates$profile_id %in% profiles, ]
  multiplier <- switch(scale, raw_dollars = .001, capped_dollars = .001,
                       positive_assets = 100, asinh_dollars = 1)
  units <- switch(scale, raw_dollars = "$1,000", capped_dollars = "capped $1,000",
                  positive_assets = "percentage points", asinh_dollars = "asinh dollars")
  outcome_title <- switch(scale, raw_dollars = "Net financial assets",
    capped_dollars = "Capped net financial assets", positive_assets = "Positive net financial assets",
    asinh_dollars = "Asinh net financial assets")
  prefix <- if (validation) "VALIDATION: " else ""
  draw_intervals <- function() {
    par(mar = c(6.5, 5, 3.5, 1), las = 1)
    limits <- range(c(rows$ci_lo, rows$ci_hi, 0)) * multiplier
    if (diff(limits) == 0) limits <- c(-1, 1)
    plot(NA, xlim = limits, ylim = c(.4, length(profiles) + 1),
      yaxt = "n", ylab = "", xlab = paste0("Eligibility contrast (", units, ")"),
      main = paste0(prefix, outcome_title, ": eligibility contrast"), bty = "l")
    y <- rev(seq_along(profiles)); axis(2, at = y, labels = incomes)
    abline(v = 0, col = "grey60")
    for (k in seq_along(methods)) {
      r <- rows[rows$method == methods[k], ]; r <- r[match(profiles, r$profile_id), ]
      yy <- y + c(.14, 0, -.14)[k]
      segments(r$ci_lo * multiplier, yy, r$ci_hi * multiplier, yy, col = colors[k])
      points(r$tau * multiplier, yy, col = colors[k], pch = c(15,16,17)[k])
    }
    legend("topleft", labels, col = colors, pch = c(15,16,17), bty = "n", cex = .8)
    mtext("Pointwise normal-reference intervals; same final local weights across training schemes.",
          side = 1, line = 4.8, cex = .75)
  }
  draw_concentration <- function() {
    par(mar = c(6.5, 5, 3.5, 1), las = 1)
    plot(NA, xlim = c(0, 100), ylim = c(.4, length(profiles) + 1), yaxt = "n",
      ylab = "", xlab = "Squared centered pseudo-value share (%)",
      main = paste0(prefix, outcome_title, ": frozen covariance"), bty = "l")
    y <- rev(seq_along(profiles)); axis(2, at = y, labels = incomes)
    for (k in seq_along(methods)) {
      r <- frozen$information[frozen$information$outcome_scale == scale &
        frozen$information$method == methods[k], ]; r <- r[match(profiles, r$profile_id), ]
      yy <- y + c(.14, 0, -.14)[k]
      segments(100 * r$jk_top1_sq_share, yy, 100 * r$jk_top5_sq_share, yy, col = colors[k])
      points(100 * r$jk_top1_sq_share, yy, pch = 16, col = colors[k])
      points(100 * r$jk_top5_sq_share, yy, pch = 1, col = colors[k])
    }
    legend("topleft", labels, col = colors, pch = 16, bty = "n", cex = .8)
    mtext("Filled: largest one; open: largest five. Fitted covariance attribution by deletion index.",
          side = 1, line = 4.8, cex = .75)
  }
  for (kind in c("intervals", "concentration")) {
    pdf(file.path(output_dir, paste0(scale, "_", kind, ".pdf")), width = 8, height = 4.5)
    tryCatch(if (kind == "intervals") draw_intervals() else draw_concentration(),
             finally = dev.off())
  }
}
