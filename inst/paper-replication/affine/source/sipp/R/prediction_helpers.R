# Controlled prediction paths. All-tree and OOB paths share forest objects.
validate_sipp_401k_prediction_inputs <- function(sipp, config) {
  d <- sipp$data
  if (config$n_folds < 2L || config$n_folds > nrow(d) ||
      any(table(factor(d$D, levels = 0:1)) < 2L))
    stop("Each treatment arm needs at least two observations and valid folds.", call. = FALSE)
  if (length(config$propensity_clip) != 1L || !is.finite(config$propensity_clip) ||
      config$propensity_clip <= 0 || config$propensity_clip >= .5)
    stop("Propensity clip must be between zero and one half.", call. = FALSE)
  if (length(config$alpha) != 1L || !is.finite(config$alpha) || config$alpha <= 0 || config$alpha >= 1)
    stop("Alpha must be between zero and one.", call. = FALSE)
  needed <- unique(c(sipp$baseline_w_cols, unname(sipp$outcome_cols), "D", "row_id", sipp$x_cols))
  if (!all(needed %in% names(d)) || any(!is.finite(as.matrix(d[, needed]))) ||
      anyDuplicated(d$row_id) || any(!d$D %in% 0:1))
    stop("Fitting requires finite aligned rows, binary treatment and unique row IDs.", call. = FALSE)
  orders <- pmax(1, floor(config$dnn_lambda_grid * nrow(d)^((length(sipp$x_cols) + 1) /
                                                         (length(sipp$x_cols) + 2))))
  if (any(!is.finite(orders)) || any(orders >= nrow(d)) || config$dnn_k != 1L ||
      anyDuplicated(config$dnn_lambda_grid))
    stop("Declared order grid must support distinct candidates and k=1 retained deletion.", call. = FALSE)
  invisible(TRUE)
}

sipp_401k_regression_forest <- function(X, Y, config, seed) {
  grf::regression_forest(X = X, Y = Y, num.trees = config$rf_num_trees,
                         seed = seed, num.threads = config$num_threads)
}

fit_sipp_401k_prediction_paths <- function(sipp, config) {
  validate_sipp_401k_prediction_inputs(sipp, config)
  d <- sipp$data; n <- nrow(d); X <- as.matrix(d[, sipp$baseline_w_cols, drop = FALSE])
  set.seed(config$seed)
  folds <- sample(rep_len(seq_len(config$n_folds), n))
  for (k in seq_len(config$n_folds))
    if (any(table(factor(d$D[folds != k], levels = 0:1)) < 2L))
      stop("A cross-fit training fold has inadequate treatment-arm support.", call. = FALSE)
  message("Fitting shared full/OOB propensity forest ...")
  ef <- sipp_401k_regression_forest(X, d$D, config, config$seed)
  paths <- list(
    fullsample = list(e_raw = as.numeric(predict(ef, newdata = X)$predictions)),
    oob = list(e_raw = as.numeric(predict(ef)$predictions)),
    crossfit = list(e_raw = numeric(n), fold_ids = folds, folds = list()))
  models <- list(fullsample = list(propensity = ef, outcomes = list()), crossfit = list())
  for (mechanism in names(paths)) {
    paths[[mechanism]]$m0_hat <- paths[[mechanism]]$m1_hat <-
      setNames(lapply(sipp$outcome_cols, function(x) numeric(n)), names(sipp$outcome_cols))
  }
  for (scale in names(sipp$outcome_cols)) {
    message("Fitting shared full/OOB outcome forests: ", scale)
    for (arm in 0:1) {
      in_arm <- d$D == arm
      forest <- sipp_401k_regression_forest(X[in_arm, , drop = FALSE],
        d[[sipp$outcome_cols[[scale]]]][in_arm], config, config$seed + 11L + arm)
      full <- as.numeric(predict(forest, newdata = X)$predictions)
      oob <- full; oob[in_arm] <- as.numeric(predict(forest)$predictions)
      field <- paste0("m", arm, "_hat")
      paths$fullsample[[field]][[scale]] <- full
      paths$oob[[field]][[scale]] <- oob
      models$fullsample$outcomes[[scale]][[paste0("m", arm)]] <- forest
    }
  }
  for (k in seq_len(config$n_folds)) {
    message("Fitting held-out nuisance fold ", k, "/", config$n_folds)
    te <- folds == k; tr <- !te
    ef <- sipp_401k_regression_forest(X[tr, , drop = FALSE], d$D[tr], config, config$seed + k)
    paths$crossfit$e_raw[te] <- as.numeric(predict(ef, newdata = X[te, , drop = FALSE])$predictions)
    paths$crossfit$folds[[k]] <- list(training_row_ids = d$row_id[tr],
      heldout_row_ids = d$row_id[te], seed = config$seed + k,
      training_arm_sizes = table(factor(d$D[tr], levels = 0:1)))
    models$crossfit[[k]] <- list(propensity = ef, outcomes = list())
    for (scale in names(sipp$outcome_cols)) for (arm in 0:1) {
      train_arm <- tr & d$D == arm
      forest <- sipp_401k_regression_forest(X[train_arm, , drop = FALSE],
        d[[sipp$outcome_cols[[scale]]]][train_arm], config, config$seed + k + 100L * (arm + 1L))
      field <- paste0("m", arm, "_hat")
      paths$crossfit[[field]][[scale]][te] <-
        as.numeric(predict(forest, newdata = X[te, , drop = FALSE])$predictions)
      models$crossfit[[k]]$outcomes[[scale]][[paste0("m", arm)]] <- forest
    }
  }
  for (mechanism in names(paths)) {
    a <- paths[[mechanism]]
    a$e_hat <- clip_propensity(a$e_raw, config$propensity_clip)$e_hat
    a$e_hat_raw_range <- range(a$e_raw)
    a$clipped <- a$e_hat != a$e_raw
    a$row_id <- d$row_id
    for (scale in names(sipp$outcome_cols)) {
      bounds <- sipp$specification$outcomes[[scale]]$bounds
      if (!is.null(bounds)) for (field in c("m0_hat", "m1_hat"))
        a[[field]][[scale]] <- pmax(bounds[1], pmin(bounds[2], a[[field]][[scale]]))
    }
    if (any(!is.finite(c(a$e_hat, unlist(a$m0_hat), unlist(a$m1_hat)))))
      stop("Nuisance predictions contain nonfinite values; increase trees or inspect support.", call. = FALSE)
    paths[[mechanism]] <- a
  }
  list(nuisances = paths, models = models)
}

make_sipp_401k_signal <- function(data, y_col, nuisance, scale) {
  e <- nuisance$e_hat; m0 <- nuisance$m0_hat[[scale]]; m1 <- nuisance$m1_hat[[scale]]
  contrast <- m1 - m0
  eligible <- data$D * (data[[y_col]] - m1) / e
  ineligible <- -(1 - data$D) * (data[[y_col]] - m0) / (1 - e)
  data.frame(row_id = data$row_id, e_raw = nuisance$e_raw, e_hat = e,
             clipped = nuisance$clipped, m0_hat = m0, m1_hat = m1,
             contrast = contrast, eligible_residual = eligible,
             ineligible_residual = ineligible, phi = contrast + eligible + ineligible)
}

fit_sipp_401k_affine_profiles <- function(sipp, config, run_info = list()) {
  pred <- fit_sipp_401k_prediction_paths(sipp, config)
  tuning <- RdistNN::tuning_dnn_paper(lambda = config$dnn_lambda,
                                     lambda_grid = config$dnn_lambda_grid)
  estimates <- diagnostics <- list()
  grids <- weights <- signals <- list()
  for (scale in names(sipp$outcome_cols)) for (mechanism in names(pred$nuisances)) {
    nuisance <- pred$nuisances[[mechanism]]
    signal <- make_sipp_401k_signal(sipp$data, sipp$outcome_cols[[scale]], nuisance, scale)
    scalar_nuisance <- list(e_hat = nuisance$e_hat, m0_hat = nuisance$m0_hat[[scale]],
                            m1_hat = nuisance$m1_hat[[scale]])
    fit <- fit_dnn_arm(sipp$data, sipp, sipp$outcome_cols[[scale]], scalar_nuisance,
      tuning_mode = if (mechanism == "fullsample") "grid" else "default",
      k = config$dnn_k, tuning = tuning)
    method <- paste0("dnn_", mechanism)
    rows <- make_estimate_rows(sipp$profiles, method, fit$tau, fit$se, config$alpha,
      ess = fit$ess, extra = list(s = fit$meta$s, n_local = NA_integer_))
    rows$outcome_scale <- scale
    estimates[[length(estimates) + 1L]] <- rows
    diagnostics[[length(diagnostics) + 1L]] <- profile_weight_diagnostics(
      weights = fit$W, d = sipp$data$D, e_hat = nuisance$e_hat,
      profiles = sipp$profiles, method = method, outcome_scale = scale,
      x_mat = as.matrix(sipp$data[, sipp$x_cols, drop = FALSE]), x0_mat = sipp$X0,
      raw_data = sipp$data, pseudo_outcome = signal$phi, tau_hat = fit$tau)
    index <- length(diagnostics)
    diagnostics[[index]]$raw_propensity_min <- vapply(seq_len(ncol(fit$W)),
      function(j) min(nuisance$e_raw[fit$W[,j] > 1e-12]), numeric(1))
    diagnostics[[index]]$raw_propensity_max <- vapply(seq_len(ncol(fit$W)),
      function(j) max(nuisance$e_raw[fit$W[,j] > 1e-12]), numeric(1))
    weights[[scale]][[mechanism]] <- fit$W
    weights[[scale]][[paste0("pseudo_outcome_", mechanism)]] <- signal$phi
    signals[[scale]][[mechanism]] <- signal
    if (mechanism == "fullsample") grids[[scale]] <- fit$tuning_grid
  }
  source <- sipp$source
  source$n_fit <- nrow(sipp$data); source$fit_sample_rows <- sipp$data$row_id
  list(schema_version = 2L, estimates = do.call(rbind, estimates),
       method_metadata = sipp_401k_method_metadata("affine"),
       profile_diagnostics = do.call(rbind, diagnostics), tuning_grids = grids,
       dnn_weights = weights, analysis_data = make_analysis_data_payload(sipp),
       nuisances = pred$nuisances, signals = signals, models = if (config$save_models) pred$models else NULL,
       profiles = sipp$profiles, specification = sipp$specification,
       config = config, run_info = run_info, source = source)
}
