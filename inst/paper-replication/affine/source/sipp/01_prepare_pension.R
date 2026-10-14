# Prepare the SIPP 401(k) sample for RdistNN CATE analysis.
#
# This script intentionally keeps all objects in memory. It does not write
# derived data, tables, figures, or caches.

.sipp_prepare_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.sipp_prepare_file <- if (length(.sipp_prepare_files)) tail(.sipp_prepare_files, 1L)[[1L]] else NULL
if (is.null(.sipp_prepare_file)) {
  .sipp_prepare_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
}
source(file.path(dirname(normalizePath(.sipp_prepare_file)), "R", "analysis_specification.R"),
       local = TRUE)

required_packages <- c("hdm")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install required package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

load_pension_data <- function() {
  pension_env <- new.env(parent = emptyenv())
  utils::data("pension", package = "hdm", envir = pension_env)
  as.data.frame(pension_env$pension)
}

check_required_columns <- function(data, cols) {
  missing_cols <- setdiff(cols, names(data))
  if (length(missing_cols) > 0L) {
    stop(
      "The pension data are missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

as_binary01 <- function(x, name) {
  bad <- !is.na(x) & !(x %in% c(0, 1))

  if (any(bad)) {
    stop(sprintf("Column `%s` must be coded as 0/1.", name), call. = FALSE)
  }

  as.integer(x)
}

standardize_with_stats <- function(x, center = mean(x), scale = stats::sd(x)) {
  if (!is.finite(center) || !is.finite(scale) || scale <= 0) {
    stop("Cannot standardize a variable with non-finite or zero scale.", call. = FALSE)
  }

  (x - center) / scale
}

standardize_sample <- function(x) {
  center <- mean(x)
  scale <- stats::sd(x)
  out <- standardize_with_stats(x, center = center, scale = scale)
  attr(out, "scaled:center") <- center
  attr(out, "scaled:scale") <- scale
  out
}

make_profile_grid <- function(data, age_value = 40, spec = make_sipp_401k_analysis_spec()) {
  income_probs <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)
  income_values <- if (!is.null(spec$profile_incomes)) spec$profile_incomes else as.numeric(
    stats::quantile(data$inc, probs = income_probs, names = FALSE, type = 7)
  )

  profile_grid <- data.frame(
    profile_id = if (is.null(spec$profile_incomes))
      sprintf("%s_p%02d", if (spec$reporting == "income") "income" else paste0("age", age_value, "_inc"),
              round(100 * income_probs)) else
      paste0(if (spec$reporting == "income") "income_" else paste0("age", age_value, "_inc_"),
             format(income_values, scientific = FALSE, trim = TRUE)),
    age = rep(if (spec$reporting == "income") NA_real_ else age_value, length(income_values)),
    inc = income_values,
    income_quantile = if (is.null(spec$profile_incomes)) income_probs else rep(NA_real_, length(income_values)),
    row.names = NULL
  )

  profile_grid$asinh_inc <- asinh(profile_grid$inc)
  profile_grid$age_z <- standardize_with_stats(
    profile_grid$age,
    center = attr(data$age_z, "scaled:center"),
    scale = attr(data$age_z, "scaled:scale")
  )
  profile_grid$asinh_inc_z <- standardize_with_stats(
    profile_grid$asinh_inc,
    center = attr(data$asinh_inc_z, "scaled:center"),
    scale = attr(data$asinh_inc_z, "scaled:scale")
  )

  profile_grid
}

prepare_sipp_401k <- function(spec = make_sipp_401k_analysis_spec(), raw = NULL) {
  spec <- do.call(make_sipp_401k_analysis_spec, spec)
  if (is.null(raw)) raw <- load_pension_data()
  n_raw <- nrow(raw)
  original_id <- seq_len(n_raw)

  required_cols <- c(
    "net_tfa", "e401", "p401", "inc", "age", "educ", "fsize", "marr",
    "male", "twoearn", "nohs", "hs", "smcol", "col",
    if (spec$adjustment != "demographic") c("db", "hown"),
    if (spec$adjustment == "ira") "pira"
  )
  check_required_columns(raw, required_cols)

  keep <- stats::complete.cases(raw[, required_cols])
  dropped <- sum(!keep)
  age_keep <- if (spec$reporting == "age40") raw$age == 40 else rep(TRUE, n_raw)
  keep <- keep & !is.na(age_keep) & age_keep
  raw <- raw[keep, , drop = FALSE]
  if (nrow(raw) < 2L) stop("Too few complete households for this specification.", call. = FALSE)
  optional_binary <- function(name) if (name %in% names(raw)) as_binary01(raw[[name]], name) else rep(NA_integer_, nrow(raw))

  data <- data.frame(
    row_id = as.integer(original_id[keep]),
    Y = as.numeric(raw$net_tfa),
    Y_asinh = asinh(as.numeric(raw$net_tfa)),
    D = as_binary01(raw$e401, "e401"),
    P = as_binary01(raw$p401, "p401"),
    net_tfa = as.numeric(raw$net_tfa),
    e401 = as_binary01(raw$e401, "e401"),
    p401 = as_binary01(raw$p401, "p401"),
    age = as.numeric(raw$age),
    inc = as.numeric(raw$inc),
    educ = as.numeric(raw$educ),
    fsize = as.numeric(raw$fsize),
    marr = as_binary01(raw$marr, "marr"),
    male = as_binary01(raw$male, "male"),
    twoearn = as_binary01(raw$twoearn, "twoearn"),
    db = optional_binary("db"),
    hown = optional_binary("hown"),
    pira = optional_binary("pira"),
    nohs = as_binary01(raw$nohs, "nohs"),
    hs = as_binary01(raw$hs, "hs"),
    smcol = as_binary01(raw$smcol, "smcol"),
    col = as_binary01(raw$col, "col"),
    row.names = NULL
  )

  data$asinh_inc <- asinh(data$inc)
  data$age2 <- data$age^2
  data$asinh_inc2 <- data$asinh_inc^2
  data$Y_positive <- as.numeric(data$Y > 0)
  data$Y_capped <- pmax(-spec$asset_cap, pmin(spec$asset_cap, data$Y))

  age_scale <- if (spec$reporting == "age_income") {
    if (is.null(spec$metric_scales)) stats::sd(data$age) else spec$metric_scales[1L]
  } else 1
  income_scale <- if (is.null(spec$metric_scales)) stats::sd(data$asinh_inc) else tail(spec$metric_scales, 1L)
  data$age_z <- standardize_with_stats(data$age, scale = age_scale)
  data$asinh_inc_z <- standardize_with_stats(data$asinh_inc, scale = income_scale)
  attr(data$age_z, "scaled:center") <- mean(data$age)
  attr(data$age_z, "scaled:scale") <- age_scale
  attr(data$asinh_inc_z, "scaled:center") <- mean(data$asinh_inc)
  attr(data$asinh_inc_z, "scaled:scale") <- income_scale

  x_cols <- if (spec$reporting == "age_income") c("age_z", "asinh_inc_z") else "asinh_inc_z"
  baseline_w_cols <- sipp_401k_adjustment_columns(spec$adjustment)
  sensitivity_w_cols <- unique(c(baseline_w_cols, "pira"))

  profiles <- make_profile_grid(data, spec = spec)
  X0 <- as.matrix(profiles[, x_cols, drop = FALSE])
  rownames(X0) <- profiles$profile_id

  list(
    data = data,
    x_cols = x_cols,
    y_col = "Y",
    t_col = "D",
    outcome_cols = c(raw_dollars = "Y", asinh_dollars = "Y_asinh",
                     positive_assets = "Y_positive", capped_dollars = "Y_capped"),
    participation_col = "P",
    baseline_w_cols = baseline_w_cols,
    sensitivity_w_cols = sensitivity_w_cols,
    profiles = profiles,
    X0 = X0,
    specification = c(spec, list(
      profile_rule = if (is.null(spec$profile_incomes)) "sample_quantiles" else "literal",
      target = switch(spec$reporting, age_income = "Eligibility CATE conditional on age and income",
                      income = "Eligibility CATE conditional on income, averaged over age",
                      age40 = "Eligibility CATE conditional on income within age40 households"),
      adjustment_columns = baseline_w_cols,
      metric = list(columns = x_cols, centers = if (length(x_cols) == 2L)
                      c(mean(data$age), mean(data$asinh_inc)) else mean(data$asinh_inc),
                    scales = if (length(x_cols) == 2L) c(age_scale, income_scale) else income_scale,
                    scale_rule = if (is.null(spec$metric_scales)) "sample_sd" else "fixed",
                    tie_rule = "input row order"),
      outcomes = sipp_401k_outcome_definitions(spec$asset_cap))),
    source = list(
      package = "hdm",
      dataset = "pension",
      hdm_version = as.character(utils::packageVersion("hdm")),
      n_raw = n_raw,
      n_analysis = nrow(data),
      n_dropped_incomplete = dropped,
      n_excluded_age = sum(!is.na(age_keep) & !age_keep),
      original_row_ids = original_id
    )
  )
}

print_sipp_401k_setup <- function(x) {
  data <- x$data

  cat("SIPP 401(k) setup for RdistNN\n")
  cat("--------------------------------\n")
  cat("Dataset: hdm::pension\n")
  cat("hdm version:", x$source$hdm_version, "\n")
  cat("Rows:", x$source$n_analysis, "\n")
  cat("Dropped incomplete rows:", x$source$n_dropped_incomplete, "\n")
  cat("Outcome column:", x$y_col, "(net_tfa dollars)\n")
  cat("Treatment column:", x$t_col, "(e401 eligibility)\n")
  cat("RdistNN x_cols:", paste(x$x_cols, collapse = ", "), "\n\n")

  summary_table <- data.frame(
    statistic = c(
      "treatment_share",
      "participation_share",
      "age_min",
      "age_median",
      "age_max",
      "income_p05",
      "income_p10",
      "income_p25",
      "income_median",
      "income_p75",
      "income_p90",
      "income_p95",
      "net_tfa_median"
    ),
    value = c(
      mean(data$D),
      mean(data$P),
      min(data$age),
      stats::median(data$age),
      max(data$age),
      stats::quantile(data$inc, 0.05, names = FALSE),
      stats::quantile(data$inc, 0.10, names = FALSE),
      stats::quantile(data$inc, 0.25, names = FALSE),
      stats::median(data$inc),
      stats::quantile(data$inc, 0.75, names = FALSE),
      stats::quantile(data$inc, 0.90, names = FALSE),
      stats::quantile(data$inc, 0.95, names = FALSE),
      stats::median(data$Y)
    ),
    row.names = NULL
  )

  print(summary_table, row.names = FALSE)
  cat("\nProfiles prepared for future X0 use:\n")
  print(x$profiles[, c("profile_id", "age", "inc", "age_z", "asinh_inc_z")], row.names = FALSE)

  invisible(x)
}

run_sipp_401k_prepare_script <- function() {
  sipp_401k <- prepare_sipp_401k()
  print_sipp_401k_setup(sipp_401k)
  invisible(sipp_401k)
}

.sipp_401k_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]
if (!is.na(.sipp_401k_file_arg) &&
    basename(sub("^--file=", "", .sipp_401k_file_arg)) == "01_prepare_pension.R") {
  sipp_401k <- run_sipp_401k_prepare_script()
}
rm(.sipp_401k_file_arg)
