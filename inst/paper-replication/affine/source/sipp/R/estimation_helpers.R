# Shared estimation helpers for the SIPP 401(k) CATE illustration.

.sipp_helper_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.sipp_helper_dir <- dirname(normalizePath(tail(.sipp_helper_files, 1L)[[1L]]))
source(file.path(.sipp_helper_dir, "analysis_specification.R"), local = TRUE)

sipp_401k_estimation_methods <- function(family = "affine") {
  sipp_401k_method_metadata(family)$method
}

sipp_401k_method_metadata <- function(family = "affine") {
  if (family == "affine") return(data.frame(
    method = c("dnn_fullsample", "dnn_oob", "dnn_crossfit"),
    role = c("full-fit comparison", "historical OOB benchmark", "sensitivity"),
    family = rep("DNN AIPW conditional CATE", 3),
    nuisance = c("full-sample all-tree grf predictions", "single-fit own-arm OOB grf predictions",
                 "cross-fitted held-out grf predictions"),
    localization = rep("common DNN rank weights in declared reporting variables", 3),
    inference = rep("frozen-score DNN delete-one jackknife; pointwise normal reference", 3)))
  if (family != "legacy") stop("Unknown method family.", call. = FALSE)
  data.frame(
    method = c(
      "dnn_fullsample",
      "dnn_crossfit",
      "causal_forest"
    ),
    role = c(
      "primary",
      "sensitivity",
      "competitor"
    ),
    family = c(
      "DNN conditional-Z CATE",
      "DNN conditional-Z CATE",
      "causal forest"
    ),
    nuisance = c(
      "full-sample OOB grf nuisances",
      "cross-fitted grf nuisances",
      "cross-fitted grf Y.hat/W.hat"
    ),
    localization = c(
      "DNN rank weights in profile variables",
      "DNN rank weights in profile variables",
      "forest prediction at profile variables"
    ),
    inference = c(
      "DNN delete-one jackknife",
      "DNN delete-one jackknife",
      "grf variance estimate"
    ),
    stringsAsFactors = FALSE
  )
}

make_sipp_401k_fit_config <- function(
  seed = 20260709L,
  n_folds = 5L,
  propensity_clip = 0.01,
  alpha = 0.05,
  rf_num_trees = 2000L,
  cf_num_trees = 4000L,
  dnn_k = 1L,
  dnn_lambda = 0.4,
  dnn_lambda_grid = c(0.3, 0.4, 0.5, 0.75),
  method_family = "affine",
  save_models = TRUE,
  num_threads = 4L
) {
  list(
    seed = as.integer(seed),
    n_folds = as.integer(n_folds),
    propensity_clip = propensity_clip,
    alpha = alpha,
    rf_num_trees = as.integer(rf_num_trees),
    cf_num_trees = as.integer(cf_num_trees),
    dnn_k = as.integer(dnn_k),
    dnn_lambda = as.numeric(dnn_lambda),
    dnn_lambda_grid = as.numeric(dnn_lambda_grid),
    method_family = method_family,
    save_models = save_models,
    num_threads = as.integer(num_threads)
  )
}

require_sipp_401k_fit_packages <- function(packages = c("hdm", "grf", "RdistNN")) {
  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_packages) > 0L) {
    stop(
      "Install required package(s): ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

detect_sipp_401k_rdistnn_source <- function(script_dir = getwd(), local_repo = NULL) {
  if (is.null(local_repo)) {
    local_repo <- file.path(script_dir, "..", "..", "..", "..", "RdistNN")
  }

  installed_version <- if (requireNamespace("RdistNN", quietly = TRUE)) {
    as.character(utils::packageVersion("RdistNN"))
  } else {
    NA_character_
  }
  installed_path <- tryCatch(
    normalizePath(find.package("RdistNN"), winslash = "/", mustWork = TRUE),
    error = function(e) NA_character_
  )

  local_present <- dir.exists(local_repo) &&
    file.exists(file.path(local_repo, "DESCRIPTION"))
  local_repo_path <- if (local_present) {
    normalizePath(local_repo, winslash = "/", mustWork = TRUE)
  } else {
    NA_character_
  }
  local_desc <- if (local_present) {
    read.dcf(file.path(local_repo_path, "DESCRIPTION"))
  } else {
    matrix(character(), nrow = 0L, ncol = 0L)
  }
  local_version <- if (local_present && "Version" %in% colnames(local_desc)) {
    unname(local_desc[1L, "Version"])
  } else {
    NA_character_
  }

  list(
    installed_version = installed_version,
    installed_path = installed_path,
    local_present = local_present,
    local_repo_path = local_repo_path,
    local_version = local_version,
    local_git_sha = sipp_401k_git_value(local_repo_path, c("rev-parse", "--short", "HEAD")),
    local_git_dirty = sipp_401k_git_dirty(local_repo_path)
  )
}

sipp_401k_git_value <- function(repo_path, args) {
  if (is.na(repo_path) || !dir.exists(repo_path)) {
    return(NA_character_)
  }
  out <- tryCatch(
    suppressWarnings(system2("git", c("-C", repo_path, args), stdout = TRUE, stderr = TRUE)),
    error = function(e) character()
  )
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0L) {
    return(NA_character_)
  }
  if (length(out) == 0L) {
    return(NA_character_)
  }
  out[1L]
}

sipp_401k_git_dirty <- function(repo_path) {
  if (is.na(repo_path) || !dir.exists(repo_path)) {
    return(NA)
  }
  out <- tryCatch(
    suppressWarnings(system2("git", c("-C", repo_path, "status", "--porcelain"),
      stdout = TRUE, stderr = TRUE
    )),
    error = function(e) character()
  )
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0L) {
    return(NA)
  }
  length(out) > 0L
}

parse_sipp_401k_fit_args <- function(args = character()) {
  defaults <- make_sipp_401k_fit_config()
  opts <- list(
    smoke = FALSE,
    sample_n = NA_integer_,
    outcomes = names(c(raw_dollars = "Y", asinh_dollars = "Y_asinh")),
    cache_dir = NULL,
    run_id = NULL,
    spec = make_sipp_401k_analysis_spec(),
    write_cache = TRUE,
    config = defaults
  )

  for (arg in args) {
    key <- sub("=.*$", "", arg)
    value <- sub("^[^=]*=", "", arg)
    if (key %in% c("--adjustment", "--reporting", "--profile-incomes",
                   "--metric-scales", "--asset-cap")) {
      field <- switch(key, "--adjustment" = "adjustment", "--reporting" = "reporting",
                      "--profile-incomes" = "profile_incomes",
                      "--metric-scales" = "metric_scales", "--asset-cap" = "asset_cap")
      opts$spec[[field]] <- if (field %in% c("adjustment", "reporting")) value else
        suppressWarnings(as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]]))
      next
    }
    if (key == "--method-family") { opts$config$method_family <- value; next }
    if (key == "--run-id") { opts$run_id <- value; next }
    if (key == "--num-threads") {
      opts$config$num_threads <- suppressWarnings(as.numeric(value)); next
    }
    if (key %in% c("--dnn-lambda", "--dnn-lambda-grid")) {
      field <- if (key == "--dnn-lambda") "dnn_lambda" else "dnn_lambda_grid"
      opts$config[[field]] <- suppressWarnings(as.numeric(strsplit(value,",",fixed=TRUE)[[1L]]))
      next
    }
    if (arg == "--no-models") { opts$config$save_models <- FALSE; next }
    if (identical(arg, "--smoke")) {
      opts$smoke <- TRUE
      next
    }
    if (identical(arg, "--no-cache")) {
      opts$write_cache <- FALSE
      next
    }
    if (startsWith(arg, "--sample-n=")) {
      opts$sample_n <- as.integer(sub("^--sample-n=", "", arg))
      next
    }
    if (startsWith(arg, "--cache-dir=")) {
      opts$cache_dir <- sub("^--cache-dir=", "", arg)
      next
    }
    if (startsWith(arg, "--rf-num-trees=")) {
      opts$config$rf_num_trees <- as.integer(sub("^--rf-num-trees=", "", arg))
      next
    }
    if (startsWith(arg, "--cf-num-trees=")) {
      opts$config$cf_num_trees <- as.integer(sub("^--cf-num-trees=", "", arg))
      next
    }
    if (startsWith(arg, "--n-folds=")) {
      opts$config$n_folds <- suppressWarnings(as.numeric(sub("^--n-folds=", "", arg)))
      next
    }
    if (startsWith(arg, "--seed=")) {
      opts$config$seed <- as.integer(sub("^--seed=", "", arg))
      next
    }
    if (startsWith(arg, "--outcomes=")) {
      raw <- strsplit(sub("^--outcomes=", "", arg), ",", fixed = TRUE)[[1L]]
      raw <- trimws(raw)
      map <- c(raw = "raw_dollars", asinh = "asinh_dollars",
               positive = "positive_assets", capped = "capped_dollars")
      unknown <- setdiff(raw, names(map))
      if (length(unknown) > 0L) {
        stop(
          "Unknown outcome option(s): ",
          paste(unknown, collapse = ", "),
          ". Use raw, asinh, positive or capped.",
          call. = FALSE
        )
      }
      opts$outcomes <- unname(map[raw])
      next
    }

    stop("Unknown argument for 02_fit_cate_profiles.R: ", arg, call. = FALSE)
  }

  if (isTRUE(opts$smoke)) {
    if (is.na(opts$sample_n)) {
      opts$sample_n <- 300L
    }
    if (!any(startsWith(args, "--outcomes="))) opts$outcomes <- "raw_dollars"
    opts$config$n_folds <- min(opts$config$n_folds, 2L)
    opts$config$rf_num_trees <- min(opts$config$rf_num_trees, 50L)
    opts$config$cf_num_trees <- min(opts$config$cf_num_trees, 50L)
  }

  validate_sipp_401k_fit_options(opts)
}

validate_sipp_401k_fit_options <- function(opts) {
  opts$spec <- do.call(make_sipp_401k_analysis_spec, opts$spec)
  if (!opts$config$method_family %in% c("affine", "legacy"))
    stop("Method family must be affine or legacy.", call. = FALSE)
  if (!is.logical(opts$smoke) || length(opts$smoke) != 1L || is.na(opts$smoke)) {
    stop("`smoke` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.na(opts$sample_n) && (!is.finite(opts$sample_n) || opts$sample_n < 50L)) {
    stop("`sample_n` must be at least 50 when supplied.", call. = FALSE)
  }
  if (length(opts$outcomes) < 1L) {
    stop("At least one outcome scale must be selected.", call. = FALSE)
  }
  if (anyDuplicated(opts$outcomes) ||
      any(!opts$outcomes %in% c("raw_dollars", "asinh_dollars", "positive_assets", "capped_dollars"))) {
    stop("Outcome scales must be `raw_dollars` or `asinh_dollars`.", call. = FALSE)
  }
  if (!is.finite(opts$config$n_folds) || opts$config$n_folds < 2L ||
      opts$config$n_folds != floor(opts$config$n_folds)) {
    stop("`n_folds` must be at least 2.", call. = FALSE)
  }
  if (any(!is.finite(c(opts$config$rf_num_trees, opts$config$cf_num_trees))) ||
      opts$config$rf_num_trees < 1L || opts$config$cf_num_trees < 1L) {
    stop("Tree counts must be positive.", call. = FALSE)
  }
  if (length(opts$config$dnn_lambda) != 1L ||
      !is.finite(opts$config$dnn_lambda) || opts$config$dnn_lambda <= 0) {
    stop("`dnn_lambda` must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.finite(opts$config$num_threads) || opts$config$num_threads < 1L ||
      opts$config$num_threads != floor(opts$config$num_threads))
    stop("Thread count must be a positive finite integer.", call. = FALSE)
  if (length(opts$config$dnn_lambda_grid) < 1L ||
      anyDuplicated(opts$config$dnn_lambda_grid) ||
      any(!is.finite(opts$config$dnn_lambda_grid)) ||
      any(opts$config$dnn_lambda_grid <= 0) ||
      !any(abs(opts$config$dnn_lambda_grid - opts$config$dnn_lambda) < 1e-12)) {
    stop(
      "`dnn_lambda_grid` must be positive, finite, and contain `dnn_lambda`.",
      call. = FALSE
    )
  }

  opts
}

apply_sipp_401k_fit_options <- function(sipp, opts) {
  missing_outcomes <- setdiff(opts$outcomes, names(sipp$outcome_cols))
  if (length(missing_outcomes) > 0L) {
    stop(
      "Selected outcome scale(s) not available: ",
      paste(missing_outcomes, collapse = ", "),
      call. = FALSE
    )
  }
  sipp$outcome_cols <- sipp$outcome_cols[opts$outcomes]

  if (!("row_id" %in% names(sipp$data))) {
    sipp$data$row_id <- seq_len(nrow(sipp$data))
  }
  fit_sample_rows <- sipp$data$row_id

  if (!is.na(opts$sample_n) && opts$sample_n < nrow(sipp$data)) {
    set.seed(opts$config$seed)
    idx <- sample.int(nrow(sipp$data), opts$sample_n)
    sipp$data <- sipp$data[idx, , drop = FALSE]
    fit_sample_rows <- sipp$data$row_id
  }

  sipp$source$n_fit <- nrow(sipp$data)
  sipp$source$fit_sample_rows <- as.integer(fit_sample_rows)
  sipp$source$outcome_scales <- names(sipp$outcome_cols)
  sipp
}

make_analysis_data_payload <- function(sipp) {
  cols <- unique(c(
    "row_id",
    "Y",
    sipp$t_col,
    sipp$participation_col,
    "age",
    "inc",
    "asinh_inc",
    sipp$x_cols,
    unname(sipp$outcome_cols)
  ))
  cols <- cols[!is.na(cols) & nzchar(cols)]
  cols <- intersect(cols, names(sipp$data))
  payload <- sipp$data[, cols, drop = FALSE]
  rownames(payload) <- NULL
  payload
}

sipp_401k_kish_effective_size <- function(weights) {
  weights <- as.numeric(weights)
  if (length(weights) < 1L) {
    return(NA_real_)
  }
  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop("Kish effective-size weights must be finite and nonnegative.", call. = FALSE)
  }
  weight_mass <- sum(weights)
  squared_mass <- sum(weights^2)
  if (weight_mass <= 0 || squared_mass <= 0) {
    return(NA_real_)
  }
  weight_mass^2 / squared_mass
}

sipp_401k_top_absolute_contribution_share <- function(contributions, top_n = 5L) {
  contributions <- abs(as.numeric(contributions))
  if (length(contributions) < 1L || any(!is.finite(contributions))) {
    stop("Score contributions must be finite.", call. = FALSE)
  }
  top_n <- as.integer(top_n)
  if (length(top_n) != 1L || is.na(top_n) || top_n < 1L) {
    stop("top_n must be a positive integer.", call. = FALSE)
  }
  total <- sum(contributions)
  if (total <= 0) {
    return(NA_real_)
  }
  sum(sort(contributions, decreasing = TRUE)[seq_len(min(top_n, length(contributions)))]) /
    total
}

profile_weight_diagnostics <- function(weights, d, e_hat, profiles, method,
                                       outcome_scale, x_mat = NULL,
                                       x0_mat = NULL, raw_data = NULL,
                                       pseudo_outcome = NULL, tau_hat = NULL,
                                       nonzero_tol = 1e-12, top_n = 5L,
                                       radius_mass = 0.90) {
  w_mat <- as.matrix(weights)
  if (nrow(w_mat) != length(d) && ncol(w_mat) == length(d)) {
    w_mat <- t(w_mat)
  }
  if (nrow(w_mat) != length(d)) {
    stop("Weight diagnostics require one weight row per observation.", call. = FALSE)
  }
  if (ncol(w_mat) != nrow(profiles)) {
    stop("Weight diagnostics require one weight column per profile.", call. = FALSE)
  }
  if (length(e_hat) != length(d)) {
    stop("Propensity diagnostics must align with the treatment vector.", call. = FALSE)
  }
  if (!is.finite(radius_mass) || length(radius_mass) != 1L ||
      radius_mass <= 0 || radius_mass > 1) {
    stop("`radius_mass` must be in (0, 1].", call. = FALSE)
  }

  has_geometry <- !is.null(x_mat) || !is.null(x0_mat)
  if (has_geometry) {
    x_mat <- as.matrix(x_mat)
    x0_mat <- as.matrix(x0_mat)
    if (nrow(x_mat) != length(d) || ncol(x_mat) != ncol(x0_mat) ||
        nrow(x0_mat) != nrow(profiles)) {
      stop("Geometry diagnostics require aligned observation and profile matrices.", call. = FALSE)
    }
  }
  if (!is.null(raw_data) && nrow(raw_data) != length(d)) {
    stop("Raw-data diagnostics must align with the treatment vector.", call. = FALSE)
  }
  if (!is.null(pseudo_outcome) && length(pseudo_outcome) != length(d)) {
    stop("Pseudo-outcome diagnostics must align with the treatment vector.", call. = FALSE)
  }
  if (!is.null(tau_hat) && length(tau_hat) != nrow(profiles)) {
    stop("Profile estimates must align with the profile rows.", call. = FALSE)
  }

  d <- as.integer(d)
  if (any(!(d %in% c(0L, 1L)))) {
    stop("Treatment diagnostics require a 0/1 treatment vector.", call. = FALSE)
  }

  out <- data.frame(
    outcome_scale = outcome_scale,
    profile_id = profiles$profile_id,
    method = method,
    ess = NA_real_,
    treated_ess = NA_real_,
    control_ess = NA_real_,
    treated_mass = NA_real_,
    untreated_mass = NA_real_,
    max_weight = NA_real_,
    top5_weight_share = NA_real_,
    n_nonzero = NA_integer_,
    propensity_min = NA_real_,
    propensity_max = NA_real_,
    radius_90_z = NA_real_,
    mean_distance_z = NA_real_,
    age_min_90 = NA_real_,
    age_max_90 = NA_real_,
    income_min_90 = NA_real_,
    income_max_90 = NA_real_,
    max_abs_score_share = NA_real_,
    top5_abs_score_share = NA_real_,
    row.names = NULL
  )

  for (j in seq_len(ncol(w_mat))) {
    w <- as.numeric(w_mat[, j])
    w[!is.finite(w)] <- 0
    if (any(w < -nonzero_tol)) {
      stop("Weight diagnostics require nonnegative localization weights.", call. = FALSE)
    }
    w <- pmax(w, 0)
    active <- w > nonzero_tol
    total_w <- sum(w)

    out$n_nonzero[j] <- sum(active)
    if (total_w <= 0 || sum(w^2) <= 0) {
      next
    }

    w_norm <- w / total_w
    top_idx <- seq_len(min(top_n, length(w_norm)))
    sorted_w <- sort(w_norm, decreasing = TRUE)
    prop_range <- range(e_hat[active])

    out$ess[j] <- sipp_401k_kish_effective_size(w)
    out$treated_mass[j] <- sum(w_norm[d == 1L])
    out$untreated_mass[j] <- sum(w_norm[d == 0L])
    treated_w <- w_norm[d == 1L]
    control_w <- w_norm[d == 0L]
    # Descriptive within-arm Kish effective sizes:
    #   n_eff,t = (sum_{i:T_i=t} w_i)^2 / sum_{i:T_i=t} w_i^2.
    # They need not sum to the overall DNN effective size and do not include
    # inverse-propensity factors or conditional residual variances from the
    # AIPW variance.
    out$treated_ess[j] <- sipp_401k_kish_effective_size(treated_w)
    out$control_ess[j] <- sipp_401k_kish_effective_size(control_w)
    out$max_weight[j] <- max(w_norm)
    out$top5_weight_share[j] <- sum(sorted_w[top_idx])
    out$propensity_min[j] <- prop_range[1L]
    out$propensity_max[j] <- prop_range[2L]

    if (has_geometry) {
      centered <- sweep(x_mat, 2L, x0_mat[j, ], FUN = "-")
      distance <- sqrt(rowSums(centered^2))
      distance_order <- order(distance)
      cumulative_mass <- cumsum(w_norm[distance_order])
      radius_index <- which(cumulative_mass >= radius_mass)[1L]
      radius <- distance[distance_order[radius_index]]
      inside_weight_ball <- distance <= radius + sqrt(.Machine$double.eps)
      out$radius_90_z[j] <- radius
      out$mean_distance_z[j] <- sum(w_norm * distance)

      if (!is.null(raw_data) && all(c("age", "inc") %in% names(raw_data))) {
        # These are separate coordinatewise summaries among observations in
        # the weight ball, not the bounds of a rectangle contained in the ball.
        age_range <- range(raw_data$age[inside_weight_ball], na.rm = TRUE)
        income_range <- range(raw_data$inc[inside_weight_ball], na.rm = TRUE)
        out$age_min_90[j] <- age_range[1L]
        out$age_max_90[j] <- age_range[2L]
        out$income_min_90[j] <- income_range[1L]
        out$income_max_90[j] <- income_range[2L]
      }
    }

    if (!is.null(pseudo_outcome) && !is.null(tau_hat)) {
      # This is concentration of the centered weighted score terms, not an
      # influence-function or delete-one influence diagnostic.
      score_contribution <- w_norm * (as.numeric(pseudo_outcome) - tau_hat[j])
      total_abs_score <- sum(abs(score_contribution))
      if (is.finite(total_abs_score) && total_abs_score > 0) {
        out$max_abs_score_share[j] <- max(abs(score_contribution)) / total_abs_score
        out$top5_abs_score_share[j] <-
          sipp_401k_top_absolute_contribution_share(score_contribution, top_n)
      }
    }
  }

  out
}

make_sipp_401k_diagnostic_config <- function(
  min_ess = 10,
  min_arm_mass = 0.05,
  max_top5_weight_share = 0.80,
  propensity_lower = 0.01,
  propensity_upper = 0.99,
  max_se_ratio = 3,
  max_tau_gap_se_units = 2
) {
  config <- list(
    min_ess = min_ess,
    min_arm_mass = min_arm_mass,
    max_top5_weight_share = max_top5_weight_share,
    propensity_lower = propensity_lower,
    propensity_upper = propensity_upper,
    max_se_ratio = max_se_ratio,
    max_tau_gap_se_units = max_tau_gap_se_units
  )
  validate_sipp_401k_diagnostic_config(config)
}

validate_sipp_401k_diagnostic_config <- function(config) {
  positive <- c("min_ess", "max_se_ratio", "max_tau_gap_se_units")
  for (nm in positive) {
    if (!is.finite(config[[nm]]) || config[[nm]] <= 0) {
      stop(sprintf("Diagnostic threshold `%s` must be positive.", nm), call. = FALSE)
    }
  }

  unit_interval <- c(
    "min_arm_mass", "max_top5_weight_share",
    "propensity_lower", "propensity_upper"
  )
  for (nm in unit_interval) {
    if (!is.finite(config[[nm]]) || config[[nm]] < 0 || config[[nm]] > 1) {
      stop(sprintf("Diagnostic threshold `%s` must be in [0, 1].", nm), call. = FALSE)
    }
  }
  if (config$propensity_lower >= config$propensity_upper) {
    stop("`propensity_lower` must be below `propensity_upper`.", call. = FALSE)
  }

  config
}

parse_sipp_401k_diagnostic_args <- function(args = character(), script_dir = getwd()) {
  cache_dir <- file.path(script_dir, "cache")
  opts <- list(
    fit_cache = file.path(cache_dir, "02_fit_cate_profiles.rds"),
    output = file.path(cache_dir, "03_profile_diagnostics.rds"),
    write_cache = TRUE,
    config = make_sipp_401k_diagnostic_config()
  )

  for (arg in args) {
    if (identical(arg, "--no-cache")) {
      opts$write_cache <- FALSE
      next
    }
    if (startsWith(arg, "--cache-dir=")) {
      cache_dir <- sub("^--cache-dir=", "", arg)
      opts$fit_cache <- file.path(cache_dir, "02_fit_cate_profiles.rds")
      opts$output <- file.path(cache_dir, "03_profile_diagnostics.rds")
      next
    }
    if (startsWith(arg, "--fit-cache=")) {
      opts$fit_cache <- sub("^--fit-cache=", "", arg)
      next
    }
    if (startsWith(arg, "--output=")) {
      opts$output <- sub("^--output=", "", arg)
      next
    }
    if (startsWith(arg, "--min-ess=")) {
      opts$config$min_ess <- as.numeric(sub("^--min-ess=", "", arg))
      next
    }
    if (startsWith(arg, "--min-arm-mass=")) {
      opts$config$min_arm_mass <- as.numeric(sub("^--min-arm-mass=", "", arg))
      next
    }
    if (startsWith(arg, "--max-top5-weight-share=")) {
      opts$config$max_top5_weight_share <- as.numeric(
        sub("^--max-top5-weight-share=", "", arg)
      )
      next
    }
    if (startsWith(arg, "--propensity-lower=")) {
      opts$config$propensity_lower <- as.numeric(sub("^--propensity-lower=", "", arg))
      next
    }
    if (startsWith(arg, "--propensity-upper=")) {
      opts$config$propensity_upper <- as.numeric(sub("^--propensity-upper=", "", arg))
      next
    }
    if (startsWith(arg, "--max-se-ratio=")) {
      opts$config$max_se_ratio <- as.numeric(sub("^--max-se-ratio=", "", arg))
      next
    }
    if (startsWith(arg, "--max-tau-gap-se-units=")) {
      opts$config$max_tau_gap_se_units <- as.numeric(
        sub("^--max-tau-gap-se-units=", "", arg)
      )
      next
    }

    stop("Unknown argument for 03_profile_diagnostics.R: ", arg, call. = FALSE)
  }

  if (!any(startsWith(args, "--output=")))
    opts$output <- file.path(dirname(opts$fit_cache), "03_profile_diagnostics.rds")
  opts$config <- validate_sipp_401k_diagnostic_config(opts$config)
  opts
}

make_sipp_401k_profile_diagnostics <- function(
  fit,
  config = make_sipp_401k_diagnostic_config()
) {
  config <- validate_sipp_401k_diagnostic_config(config)
  required <- c("estimates", "profile_diagnostics")
  missing <- setdiff(required, names(fit))
  if (length(missing) > 0L) {
    stop(
      "Fit object is missing required component(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  support <- flag_sipp_401k_support_diagnostics(fit$profile_diagnostics, config)
  competitor_comparison <- compare_sipp_401k_dnn_competitor(fit$estimates, config)
  tuning_summary <- summarize_sipp_401k_tuning_grids(fit$tuning_grids)
  outcome_scale_coverage <- summarize_sipp_401k_outcome_scale_coverage(fit$estimates)
  profile_screen <- summarize_sipp_401k_profile_screen(support)

  result <- list(
    support = support,
    competitor_comparison = competitor_comparison,
    tuning_summary = tuning_summary,
    outcome_scale_coverage = outcome_scale_coverage,
    profile_screen = profile_screen,
    config = config,
    source = fit$source
  )
  result$profile_screen$passes_mechanical_checks <- result$profile_screen$usable_profile
  result$fit_fingerprint <- sipp_401k_fit_fingerprint(fit)
  if (!is.null(fit$dnn_weights)) result$frozen <- analyze_sipp_401k_frozen(fit)
  result
}

flag_sipp_401k_support_diagnostics <- function(profile_diagnostics, config) {
  support <- profile_diagnostics
  if (all(c("raw_propensity_min", "raw_propensity_max") %in% names(support))) {
    propensity_min <- support$raw_propensity_min; propensity_max <- support$raw_propensity_max
    support$propensity_screen_source <- "raw predictions"
  } else {
    propensity_min <- support$propensity_min; propensity_max <- support$propensity_max
    support$propensity_screen_source <- "clipped predictions; raw unavailable"
  }
  support$low_ess <- is.na(support$ess) | support$ess < config$min_ess
  support$weak_treated_mass <- is.na(support$treated_mass) |
    support$treated_mass < config$min_arm_mass
  support$weak_untreated_mass <- is.na(support$untreated_mass) |
    support$untreated_mass < config$min_arm_mass
  support$concentrated_weights <- is.na(support$top5_weight_share) |
    support$top5_weight_share > config$max_top5_weight_share
  support$propensity_boundary <- is.na(propensity_min) | is.na(propensity_max) |
    propensity_min < config$propensity_lower | propensity_max > config$propensity_upper

  flag_cols <- c(
    "low_ess", "weak_treated_mass", "weak_untreated_mass",
    "concentrated_weights", "propensity_boundary"
  )
  support$flag_count <- rowSums(support[, flag_cols, drop = FALSE])
  support$usable_profile_method <- support$flag_count == 0L
  support
}

compare_sipp_401k_dnn_competitor <- function(
  estimates,
  config,
  dnn_methods = c("dnn_fullsample", "dnn_crossfit"),
  competitor_method = "causal_forest"
) {
  competitor <- estimates[
    estimates$method == competitor_method,
    c("outcome_scale", "profile_id", "tau", "se"),
    drop = FALSE
  ]
  names(competitor)[names(competitor) == "tau"] <- "competitor_tau"
  names(competitor)[names(competitor) == "se"] <- "competitor_se"

  out <- list()
  for (dnn_method in dnn_methods) {
    dnn <- estimates[
      estimates$method == dnn_method,
      c("outcome_scale", "profile_id", "tau", "se"),
      drop = FALSE
    ]
    names(dnn)[names(dnn) == "tau"] <- "dnn_tau"
    names(dnn)[names(dnn) == "se"] <- "dnn_se"
    rows <- merge(dnn, competitor, by = c("outcome_scale", "profile_id"), all = FALSE)
    if (nrow(rows) == 0L) {
      next
    }
    rows$dnn_method <- dnn_method
    rows$competitor_method <- competitor_method
    rows$se_ratio <- rows$dnn_se / rows$competitor_se
    rows$se_disagreement <- !is.finite(rows$se_ratio) |
      rows$se_ratio > config$max_se_ratio |
      rows$se_ratio < 1 / config$max_se_ratio
    rows$tau_gap <- rows$dnn_tau - rows$competitor_tau
    rows$tau_gap_se_units <- abs(rows$tau_gap) /
      sqrt(rows$dnn_se^2 + rows$competitor_se^2)
    rows$large_tau_gap <- !is.finite(rows$tau_gap_se_units) |
      rows$tau_gap_se_units > config$max_tau_gap_se_units
    out[[dnn_method]] <- rows[
      , c(
        "outcome_scale", "profile_id", "dnn_method", "competitor_method",
        "dnn_tau", "competitor_tau", "tau_gap", "dnn_se", "competitor_se",
        "se_ratio", "tau_gap_se_units", "se_disagreement", "large_tau_gap"
      ),
      drop = FALSE
    ]
  }

  if (length(out) == 0L) {
    return(data.frame())
  }
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

summarize_sipp_401k_tuning_grids <- function(tuning_grids) {
  cols <- c(
    "outcome_scale", "candidate_id", "lambda", "s", "is_default",
    "feasible", "ess", "max_weight", "top5_weight_share", "status", "reason"
  )
  if (is.null(tuning_grids) || length(tuning_grids) == 0L) {
    return(empty_sipp_401k_frame(cols))
  }

  out <- list()
  for (scale_name in names(tuning_grids)) {
    grid <- tuning_grids[[scale_name]]
    candidates <- grid$candidates
    if (is.null(candidates) || nrow(candidates) == 0L) {
      next
    }
    diagnostics <- grid$diagnostics
    if (!is.null(diagnostics) && nrow(diagnostics) > 0L) {
      drop_diag <- intersect(c("status", "reason"), names(diagnostics))
      diagnostics <- diagnostics[, setdiff(names(diagnostics), drop_diag), drop = FALSE]
      candidates <- merge(candidates, diagnostics, by = "candidate_id", all.x = TRUE)
    }
    candidates$outcome_scale <- scale_name
    for (nm in setdiff(cols, names(candidates))) {
      candidates[[nm]] <- NA
    }
    out[[scale_name]] <- candidates[, cols, drop = FALSE]
  }

  if (length(out) == 0L) {
    return(empty_sipp_401k_frame(cols))
  }
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

summarize_sipp_401k_outcome_scale_coverage <- function(estimates) {
  keys <- unique(estimates[, c("profile_id", "method"), drop = FALSE])
  out <- data.frame(
    profile_id = keys$profile_id,
    method = keys$method,
    n_outcome_scales = NA_integer_,
    outcome_scales = NA_character_,
    has_raw = NA,
    has_asinh = NA,
    has_both = NA,
    row.names = NULL
  )

  for (i in seq_len(nrow(keys))) {
    scales <- unique(estimates$outcome_scale[
      estimates$profile_id == keys$profile_id[i] &
        estimates$method == keys$method[i]
    ])
    out$n_outcome_scales[i] <- length(scales)
    out$outcome_scales[i] <- paste(sort(scales), collapse = ",")
    out$has_raw[i] <- "raw_dollars" %in% scales
    out$has_asinh[i] <- "asinh_dollars" %in% scales
    out$has_both[i] <- out$has_raw[i] && out$has_asinh[i]
  }

  out
}

summarize_sipp_401k_profile_screen <- function(support) {
  keys <- unique(support[, c("outcome_scale", "profile_id"), drop = FALSE])
  flag_cols <- c(
    "low_ess", "weak_treated_mass", "weak_untreated_mass",
    "concentrated_weights", "propensity_boundary"
  )
  out <- data.frame(
    outcome_scale = keys$outcome_scale,
    profile_id = keys$profile_id,
    n_methods = NA_integer_,
    flagged_methods = NA_character_,
    any_flag = NA,
    usable_profile = NA,
    row.names = NULL
  )
  for (flag_col in flag_cols) {
    out[[paste0("any_", flag_col)]] <- NA
  }

  for (i in seq_len(nrow(keys))) {
    rows <- support[
      support$outcome_scale == keys$outcome_scale[i] &
        support$profile_id == keys$profile_id[i],
      ,
      drop = FALSE
    ]
    flagged <- rows$method[rows$flag_count > 0L]
    out$n_methods[i] <- nrow(rows)
    out$flagged_methods[i] <- paste(flagged, collapse = ",")
    out$any_flag[i] <- length(flagged) > 0L
    out$usable_profile[i] <- !out$any_flag[i]
    for (flag_col in flag_cols) {
      out[[paste0("any_", flag_col)]][i] <- any(rows[[flag_col]])
    }
  }

  out
}

empty_sipp_401k_frame <- function(cols) {
  out <- as.data.frame(setNames(rep(list(logical()), length(cols)), cols))
  out[0L, , drop = FALSE]
}

clip_propensity <- function(e_hat, clip) {
  list(
    e_hat = pmin(pmax(e_hat, clip), 1 - clip),
    raw_range = range(e_hat)
  )
}

fit_propensity_fullsample <- function(data, w_cols, num_trees, seed) {
  forest <- grf::regression_forest(
    X = as.matrix(data[, w_cols, drop = FALSE]),
    Y = data$D,
    num.trees = num_trees,
    seed = seed
  )
  as.numeric(predict(forest)$predictions)
}

fit_outcome_arm_fullsample <- function(data, w_cols, y_col, arm, num_trees,
                                       seed) {
  in_arm <- data$D == arm
  if (!any(in_arm)) {
    stop(sprintf("No observations with D = %d.", arm), call. = FALSE)
  }

  w_mat <- as.matrix(data[, w_cols, drop = FALSE])
  forest <- grf::regression_forest(
    X = w_mat[in_arm, , drop = FALSE],
    Y = data[in_arm, y_col],
    num.trees = num_trees,
    seed = seed
  )

  pred <- numeric(nrow(data))
  pred[in_arm] <- as.numeric(predict(forest)$predictions)
  pred[!in_arm] <- as.numeric(
    predict(forest, newdata = w_mat[!in_arm, , drop = FALSE])$predictions
  )
  pred
}

fit_nuisances_crossfit <- function(data, w_cols, y_cols, n_folds, num_trees,
                                   seed) {
  n <- nrow(data)
  fold_ids <- sample(rep_len(seq_len(n_folds), n))
  w_mat <- as.matrix(data[, w_cols, drop = FALSE])

  e_hat <- numeric(n)
  m0_hat <- lapply(y_cols, function(y) numeric(n))
  m1_hat <- lapply(y_cols, function(y) numeric(n))

  for (fold_k in seq_len(n_folds)) {
    te <- fold_ids == fold_k
    tr <- !te
    fold_seed <- seed + fold_k

    e_forest <- grf::regression_forest(
      X = w_mat[tr, , drop = FALSE],
      Y = data$D[tr],
      num.trees = num_trees,
      seed = fold_seed
    )
    e_hat[te] <- as.numeric(
      predict(e_forest, newdata = w_mat[te, , drop = FALSE])$predictions
    )

    for (scale_name in names(y_cols)) {
      y_col <- y_cols[[scale_name]]
      for (arm in c(0L, 1L)) {
        tr_arm <- tr & data$D == arm
        if (!any(tr_arm)) {
          stop(
            sprintf("Fold %d has no training observations with D = %d.", fold_k, arm),
            call. = FALSE
          )
        }
        forest <- grf::regression_forest(
          X = w_mat[tr_arm, , drop = FALSE],
          Y = data[tr_arm, y_col],
          num.trees = num_trees,
          seed = fold_seed + 100L * (arm + 1L)
        )
        pred_te <- as.numeric(
          predict(forest, newdata = w_mat[te, , drop = FALSE])$predictions
        )
        if (arm == 0L) {
          m0_hat[[scale_name]][te] <- pred_te
        } else {
          m1_hat[[scale_name]][te] <- pred_te
        }
      }
    }
  }

  list(e_hat = e_hat, m0_hat = m0_hat, m1_hat = m1_hat, fold_ids = fold_ids)
}

aipw_pseudo_outcome <- function(y, d, e_hat, m0_hat, m1_hat) {
  m1_hat - m0_hat +
    d * (y - m1_hat) / e_hat -
    (1 - d) * (y - m0_hat) / (1 - e_hat)
}

ci_bounds <- function(tau, se, alpha) {
  z <- stats::qnorm(1 - alpha / 2)
  list(lo = tau - z * se, hi = tau + z * se)
}

make_estimate_rows <- function(profiles, method, tau, se, alpha,
                               ess = NA_real_, extra = list()) {
  ci <- ci_bounds(tau, se, alpha)
  rows <- data.frame(
    profile_id = profiles$profile_id,
    method = method,
    tau = as.numeric(tau),
    se = as.numeric(se),
    ci_lo = ci$lo,
    ci_hi = ci$hi,
    ess = as.numeric(ess),
    row.names = NULL
  )
  for (nm in names(extra)) {
    rows[[nm]] <- extra[[nm]]
  }
  rows
}

fit_dnn_arm <- function(data, sipp, y_col, nuisance, tuning_mode, k,
                        s = NULL, tuning = NULL, se = TRUE,
                        jackknife = "fast") {
  dat <- data
  dat$e_hat <- nuisance$e_hat
  dat$m0_hat <- nuisance$m0_hat
  dat$m1_hat <- nuisance$m1_hat

  args <- list(
    data = dat,
    x_cols = sipp$x_cols,
    y_col = y_col,
    t_col = sipp$t_col,
    e_col = "e_hat",
    m0_col = "m0_hat",
    m1_col = "m1_hat",
    X0 = sipp$X0,
    method = "exact",
    k = k,
    tuning_mode = tuning_mode,
    se = se,
    jackknife = jackknife,
    return = c("W", "ess", "psi", "dx")
  )
  if (!is.null(s)) {
    args$s <- s
  }
  if (!is.null(tuning)) {
    args$tuning <- tuning
  }

  do.call(RdistNN::cate_dnn, args)
}

fit_causal_forest_arm <- function(data, sipp, y_col, nuisance, num_trees,
                                  seed) {
  x_mat <- as.matrix(data[, sipp$x_cols, drop = FALSE])
  y_hat <- nuisance$e_hat * nuisance$m1_hat +
    (1 - nuisance$e_hat) * nuisance$m0_hat

  cf <- grf::causal_forest(
    X = x_mat,
    Y = data[, y_col],
    W = data$D,
    Y.hat = y_hat,
    W.hat = nuisance$e_hat,
    num.trees = num_trees,
    seed = seed
  )
  pred <- predict(cf, newdata = sipp$X0, estimate.variance = TRUE)
  list(tau = pred$predictions, se = sqrt(pred$variance.estimates))
}

fit_sipp_401k_legacy_profiles <- function(sipp, config, run_info = list()) {
  data <- sipp$data
  n <- nrow(data)
  w_cols <- sipp$baseline_w_cols
  y_cols <- as.list(sipp$outcome_cols)
  source <- sipp$source
  source$n_fit <- nrow(data)
  if (is.null(source$fit_sample_rows)) {
    source$fit_sample_rows <- if ("row_id" %in% names(data)) data$row_id else seq_len(nrow(data))
  }
  source$outcome_scales <- names(y_cols)

  message("Fitting full-sample propensity forest ...")
  e_fs_raw <- fit_propensity_fullsample(
    data, w_cols,
    num_trees = config$rf_num_trees, seed = config$seed
  )
  e_fs <- clip_propensity(e_fs_raw, config$propensity_clip)

  message("Fitting cross-fitted nuisances (", config$n_folds, " folds) ...")
  cf_nuis <- fit_nuisances_crossfit(
    data, w_cols, y_cols,
    n_folds = config$n_folds,
    num_trees = config$rf_num_trees,
    seed = config$seed
  )
  e_cf <- clip_propensity(cf_nuis$e_hat, config$propensity_clip)

  x_mat <- as.matrix(data[, sipp$x_cols, drop = FALSE])

  estimates <- list()
  tuning_grids <- list()
  dnn_weights <- list()
  profile_diagnostics <- list()
  dnn_tuning <- RdistNN::tuning_dnn_paper(
    lambda = config$dnn_lambda,
    lambda_grid = config$dnn_lambda_grid
  )

  for (scale_name in names(y_cols)) {
    y_col <- y_cols[[scale_name]]
    message("Outcome scale `", scale_name, "` (column `", y_col, "`):")

    message("  full-sample outcome forests ...")
    nuis_fs <- list(
      e_hat = e_fs$e_hat,
      m0_hat = fit_outcome_arm_fullsample(
        data, w_cols, y_col,
        arm = 0L,
        num_trees = config$rf_num_trees, seed = config$seed + 11L
      ),
      m1_hat = fit_outcome_arm_fullsample(
        data, w_cols, y_col,
        arm = 1L,
        num_trees = config$rf_num_trees, seed = config$seed + 12L
      )
    )
    nuis_cf <- list(
      e_hat = e_cf$e_hat,
      m0_hat = cf_nuis$m0_hat[[scale_name]],
      m1_hat = cf_nuis$m1_hat[[scale_name]]
    )

    message("  dnn_fullsample (tuning grid) ...")
    fit_fs <- fit_dnn_arm(
      data, sipp, y_col,
      nuisance = nuis_fs,
      tuning_mode = "grid", k = config$dnn_k, tuning = dnn_tuning
    )

    message("  dnn_crossfit ...")
    fit_cf <- fit_dnn_arm(
      data, sipp, y_col,
      nuisance = nuis_cf,
      tuning_mode = "default", k = config$dnn_k, tuning = dnn_tuning
    )

    message("  causal_forest ...")
    cforest <- fit_causal_forest_arm(
      data, sipp, y_col,
      nuisance = nuis_cf,
      num_trees = config$cf_num_trees, seed = config$seed + 21L
    )

    scale_rows <- rbind(
      make_estimate_rows(
        sipp$profiles, "dnn_fullsample", fit_fs$tau, fit_fs$se,
        alpha = config$alpha, ess = fit_fs$ess,
        extra = list(s = fit_fs$meta$s, n_local = NA_integer_)
      ),
      make_estimate_rows(
        sipp$profiles, "dnn_crossfit", fit_cf$tau, fit_cf$se,
        alpha = config$alpha, ess = fit_cf$ess,
        extra = list(s = fit_cf$meta$s, n_local = NA_integer_)
      ),
      make_estimate_rows(
        sipp$profiles, "causal_forest", cforest$tau, cforest$se,
        alpha = config$alpha,
        extra = list(s = NA_integer_, n_local = NA_integer_)
      )
    )
    scale_rows$outcome_scale <- scale_name
    estimates[[scale_name]] <- scale_rows

    tuning_grids[[scale_name]] <- fit_fs$tuning_grid
    dnn_weights[[scale_name]] <- list(
      fullsample = fit_fs$W,
      crossfit = fit_cf$W,
      pseudo_outcome_fullsample = fit_fs$psi,
      pseudo_outcome_crossfit = fit_cf$psi
    )
    profile_diagnostics[[scale_name]] <- rbind(
      profile_weight_diagnostics(
        weights = fit_fs$W,
        d = data$D,
        e_hat = nuis_fs$e_hat,
        profiles = sipp$profiles,
        method = "dnn_fullsample",
        outcome_scale = scale_name,
        x_mat = x_mat,
        x0_mat = sipp$X0,
        raw_data = data,
        pseudo_outcome = fit_fs$psi,
        tau_hat = fit_fs$tau
      ),
      profile_weight_diagnostics(
        weights = fit_cf$W,
        d = data$D,
        e_hat = nuis_cf$e_hat,
        profiles = sipp$profiles,
        method = "dnn_crossfit",
        outcome_scale = scale_name,
        x_mat = x_mat,
        x0_mat = sipp$X0,
        raw_data = data,
        pseudo_outcome = fit_cf$psi,
        tau_hat = fit_cf$tau
      )
    )
  }

  estimates <- do.call(rbind, estimates)
  rownames(estimates) <- NULL
  estimates <- estimates[
    , c(
      "outcome_scale", "profile_id", "method", "tau", "se",
      "ci_lo", "ci_hi", "ess", "s", "n_local"
    )
  ]
  profile_diagnostics <- do.call(rbind, profile_diagnostics)
  rownames(profile_diagnostics) <- NULL

  list(
    estimates = estimates,
    method_metadata = sipp_401k_method_metadata("legacy"),
    profile_diagnostics = profile_diagnostics,
    tuning_grids = tuning_grids,
    dnn_weights = dnn_weights,
    analysis_data = make_analysis_data_payload(sipp),
    nuisances = list(
      fullsample = list(
        e_hat = e_fs$e_hat,
        e_hat_raw_range = e_fs$raw_range
      ),
      crossfit = list(
        e_hat = e_cf$e_hat,
        e_hat_raw_range = e_cf$raw_range,
        fold_ids = cf_nuis$fold_ids
      )
    ),
    profiles = sipp$profiles,
    config = config,
    run_info = run_info,
    source = source
  )
}

source(file.path(.sipp_helper_dir, "prediction_helpers.R"), local = TRUE)
source(file.path(.sipp_helper_dir, "frozen_diagnostics.R"), local = TRUE)
source(file.path(.sipp_helper_dir, "run_store.R"), local = TRUE)
source(file.path(.sipp_helper_dir, "affine_export.R"), local = TRUE)
fit_sipp_401k_profiles <- function(sipp, config, run_info = list()) {
  if (is.null(config$method_family) || config$method_family == "legacy")
    fit_sipp_401k_legacy_profiles(sipp, config, run_info)
  else if (config$method_family == "affine") fit_sipp_401k_affine_profiles(sipp, config, run_info)
  else stop("Unknown method family.", call. = FALSE)
}
