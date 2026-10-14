# Exact fitted-array analysis; no model fitting or RdistNN calls.
validate_sipp_401k_saved_grid <- function(fit, scale, n, kx) {
  grid <- fit$tuning_grids[[scale]]
  if (is.null(grid) || !is.data.frame(grid$candidates) || !is.data.frame(grid$estimates))
    stop("Missing complete saved tuning grid.", call. = FALSE)
  candidates <- grid$candidates; rows <- grid$estimates
  if (!all(c("candidate_id", "lambda", "s", "is_default") %in% names(candidates)) ||
      !all(c("candidate_id", "s", "target_id", "estimate", "se") %in% names(rows)))
    stop("Saved tuning grid schema is incomplete.", call. = FALSE)
  expected <- fit$config$dnn_lambda_grid
  if (length(expected) == 0L || anyDuplicated(expected) || anyDuplicated(candidates$candidate_id) ||
      nrow(candidates) != length(expected) || !setequal(candidates$lambda, expected))
    stop("Saved grid candidates do not cover the declared configuration.", call. = FALSE)
  orders <- pmax(1, floor(candidates$lambda * n^((kx + 1) / (kx + 2))))
  if (any(!is.finite(candidates$s)) || any(orders != candidates$s) ||
      any(orders >= n) || anyNA(candidates$is_default) || sum(candidates$is_default) != 1L ||
      candidates$lambda[candidates$is_default] != fit$config$dnn_lambda)
    stop("Invalid or inconsistent saved grid orders/default.", call. = FALSE)
  keys <- paste(rows$candidate_id, rows$target_id)
  expected_keys <- as.vector(outer(candidates$candidate_id, seq_len(nrow(fit$profiles)), paste))
  if (anyDuplicated(keys) || !setequal(keys, expected_keys) ||
      any(!is.finite(c(rows$s, rows$estimate, rows$se))) || any(rows$se < 0) ||
      any(rows$s != candidates$s[match(rows$candidate_id, candidates$candidate_id)]))
    stop("Missing, duplicate or nonfinite grid profile coverage.", call. = FALSE)
  grid
}

validate_sipp_401k_saved_signal <- function(fit, scale, mechanism, phi) {
  d <- fit$analysis_data; signal <- fit$signals[[scale]][[mechanism]]
  nuisance <- fit$nuisances[[mechanism]]
  fields <- c("row_id", "e_raw", "e_hat", "clipped", "m0_hat", "m1_hat",
              "contrast", "eligible_residual", "ineligible_residual", "phi")
  same <- function(a,b) length(a) == length(b) && all(is.finite(c(a,b))) &&
    max(abs(a-b)) <= 1e-10 * max(1, abs(a), abs(b))
  if (is.null(signal) || !all(fields %in% names(signal)) ||
      !identical(signal$row_id, d$row_id) || !identical(nuisance$row_id, d$row_id))
    stop("Saved nuisance/signal row IDs are inconsistent.", call. = FALSE)
  for (field in c("e_raw", "e_hat"))
    if (!same(signal[[field]], nuisance[[field]]))
      stop("Saved nuisance predictions disagree with signal predictions.", call. = FALSE)
  for (field in c("m0_hat", "m1_hat"))
    if (!same(signal[[field]], nuisance[[field]][[scale]]))
      stop("Saved outcome predictions disagree with signal components.", call. = FALSE)
  clipped <- pmax(fit$config$propensity_clip, pmin(1-fit$config$propensity_clip, signal$e_raw))
  if (!same(signal$e_hat, clipped) ||
      !identical(signal$clipped, signal$e_raw != clipped) ||
      !identical(nuisance$clipped, signal$clipped))
    stop("Saved propensity clipping predictions are inconsistent.", call. = FALSE)
  y_col <- fit$specification$outcomes[[scale]]$column
  if (is.null(y_col) || !y_col %in% names(d))
    stop("Saved component outcome definition is missing.", call. = FALSE)
  expected <- list(contrast = signal$m1_hat - signal$m0_hat,
    eligible_residual = d$D * (d[[y_col]] - signal$m1_hat) / signal$e_hat,
    ineligible_residual = -(1-d$D) * (d[[y_col]] - signal$m0_hat) / (1-signal$e_hat))
  for (field in names(expected))
    if (!same(signal[[field]], expected[[field]]))
      stop("Saved AIPW component identity is inconsistent.", call. = FALSE)
  if (!same(signal$phi, phi) || !same(rowSums(signal[,names(expected)]), phi))
    stop("Saved total signal and component sum are inconsistent.", call. = FALSE)
  invisible(TRUE)
}

sipp_401k_rank_weights <- function(n, s) {
  if (length(s) != 1L || !is.finite(s) || s < 1L || s > n || s != floor(s))
    stop("Invalid exact DNN order.", call. = FALSE)
  w <- numeric(n); w[1L] <- s / n
  if (n > 1L) for (r in seq_len(n - 1L))
    w[r + 1L] <- if (r <= n - s) w[r] * (n - r - s + 1) / (n - r) else 0
  w
}

sipp_401k_frozen_rank <- function(phi, s) {
  phi <- as.matrix(phi); n <- nrow(phi)
  if (n < 2L || s >= n || any(!is.finite(phi)))
    stop("Frozen deletion requires finite signals and order smaller than sample size.", call. = FALSE)
  w <- sipp_401k_rank_weights(n, s); wm <- sipp_401k_rank_weights(n - 1L, s)
  tau <- colSums(w * phi); loo <- matrix(0, n, ncol(phi))
  for (k in seq_len(ncol(phi))) {
    before <- c(wm * phi[seq_len(n - 1L), k], 0)
    after <- c(0, wm * phi[2:n, k])
    loo[, k] <- c(0, cumsum(before))[seq_len(n)] + c(rev(cumsum(rev(after)))[-1L], 0)
  }
  pv <- sweep(-(n - 1L) * loo, 2L, n * tau, "+")
  centered <- sweep(pv, 2L, colMeans(pv), "-")
  list(tau = tau, weights = w, loo = loo, pseudo_values = pv, centered = centered,
       covariance = crossprod(centered) / (n * (n - 1L)))
}

analyze_sipp_401k_frozen <- function(fit) {
  required <- c("analysis_data", "dnn_weights", "profiles", "estimates")
  if (!all(required %in% names(fit)))
    stop("Frozen analysis requires saved data, weights, profiles and estimates.", call. = FALSE)
  if (!is.null(fit$config$dnn_k) && fit$config$dnn_k != 1L)
    stop("This independent frozen analysis supports exact k=1 only.", call. = FALSE)
  d <- fit$analysis_data; n <- nrow(d); profiles <- fit$profiles
  if (!all(c("row_id", "D") %in% names(d)) || anyDuplicated(d$row_id) ||
      any(!d$D %in% 0:1)) stop("Invalid or unaligned analysis row identifiers.", call. = FALSE)
  cols <- fit$specification$metric$columns
  if (is.null(cols)) cols <- c("age_z", "asinh_inc_z")
  X <- as.matrix(d[, cols, drop = FALSE]); X0 <- as.matrix(profiles[, cols, drop = FALSE])
  if (any(!is.finite(c(X, X0)))) stop("Reporting geometry must be finite.", call. = FALSE)
  has_components <- identical(fit$schema_version, 2L)
  info <- routes <- paired <- tuning <- ccov <- cmeans <- list()
  ratio <- function(a, b) if (b > 0) a / b else NA_real_
  for (scale in names(fit$dnn_weights)) {
    grid <- validate_sipp_401k_saved_grid(fit, scale, n, length(cols))
    payload <- fit$dnn_weights[[scale]]
    mechanisms <- names(payload)[!startsWith(names(payload), "pseudo_outcome_")]
    for (j in seq_len(nrow(profiles))) {
      distance <- sqrt(rowSums(sweep(X, 2L, X0[j, ], "-")^2))
      ord <- order(distance, seq_len(n)); centered <- list()
      for (mechanism in mechanisms) {
        method <- paste0("dnn_", mechanism)
        w <- payload[[mechanism]]; phi <- payload[[paste0("pseudo_outcome_", mechanism)]]
        if (!is.matrix(w) || !identical(dim(w), c(n, nrow(profiles))) ||
            length(phi) != n || any(!is.finite(c(w, phi))))
          stop("Saved weights and signals must be finite and aligned.", call. = FALSE)
        row <- fit$estimates[fit$estimates$outcome_scale == scale &
          fit$estimates$profile_id == profiles$profile_id[j] & fit$estimates$method == method, ]
        if (nrow(row) != 1L) stop("Missing or duplicate estimate row.", call. = FALSE)
        z <- sipp_401k_frozen_rank(phi[ord], row$s)
        if (max(abs(z$weights - w[ord, j])) > 1e-12 ||
            abs(z$tau - row$tau) > 1e-7 ||
            abs(sqrt(z$covariance[1, 1]) - row$se) > 1e-5)
          stop("Saved weights or inference disagree with declared geometry/order.", call. = FALSE)
        if (has_components) validate_sipp_401k_saved_signal(fit, scale, mechanism, phi)
        e <- fit$nuisances[[mechanism]]$e_hat
        if (length(e) != n || any(!is.finite(e)) || any(e <= 0 | e >= 1))
          stop("Saved propensity arrays must be aligned and strictly inside (0,1).", call. = FALSE)
        pc <- z$centered[, 1L]; total <- sum(pc^2)
        h <- if (total > 0) pc^2 / total else rep(NA_real_, n)
        top <- if (total > 0) which.max(h) else 1L
        raw <- fit$nuisances[[mechanism]]$e_raw
        raw_available <- !is.null(raw)
        if (raw_available && (length(raw) != n || any(!is.finite(raw))))
          stop("Raw propensity array is not aligned.", call. = FALSE)
        a1 <- w[, j] * d$D / e; a0 <- w[, j] * (1 - d$D) / (1 - e)
        sc <- w[, j] * (phi - row$tau)
        info[[length(info) + 1L]] <- data.frame(outcome_scale = scale,
          profile_id = profiles$profile_id[j], method = method,
          prediction_mechanism = if (!has_components && mechanism == "fullsample")
            "historical single-fit/OOB" else mechanism,
          tau = z$tau, se = sqrt(z$covariance[1, 1]), s = row$s,
          ess = 1 / sum(w[, j]^2), treated_mass = sum(w[, j] * d$D),
          ipw_treated_kish = sipp_401k_kish_effective_size(a1),
          ipw_control_kish = sipp_401k_kish_effective_size(a0),
          jk_top1_sq_share = if (total > 0) max(h) else NA_real_,
          jk_top5_sq_share = if (total > 0) sum(head(sort(h, decreasing = TRUE), 5L)) else NA_real_,
          leverage_index = sum(h^2),
          top5_abs_score_share = ratio(sum(head(sort(abs(sc), decreasing = TRUE), 5L)), sum(abs(sc))),
          top_row_id = d$row_id[ord[top]], top_D = d$D[ord[top]],
          top_Y = if ("Y" %in% names(d)) d$Y[ord[top]] else NA_real_,
          top_propensity = e[ord[top]], top_rank = top, top_distance = distance[ord[top]],
          top_weight = w[ord[top], j], frozen_delete_tau = z$loo[top, 1L],
          raw_propensity_available = raw_available,
          local_clipped_mass = if (raw_available) sum(w[, j] * (raw != e)) else NA_real_,
          raw_propensity_min = if (raw_available) min(raw[w[, j] > 1e-12]) else NA_real_,
          raw_propensity_max = if (raw_available) max(raw[w[, j] > 1e-12]) else NA_real_,
          degenerate_variance = total == 0)
        routes[[length(routes) + 1L]] <- data.frame(outcome_scale = scale,
          profile_id = profiles$profile_id[j], method = method,
          variance_jackknife = z$covariance[1, 1],
          variance_squared_weights = sum(sc^2),
          variance_local_moment = sum(w[, j] * (phi - row$tau)^2) * sum(w[, j]^2))
        centered[[mechanism]] <- list(pc = pc, phi = phi, tau = row$tau,
                                     variance = z$covariance[1,1], w = w[,j])
        if (has_components) {
          signal <- fit$signals[[scale]][[mechanism]]
          fields <- c("contrast", "eligible_residual", "ineligible_residual")
          if (is.null(signal) || !identical(signal$row_id, d$row_id) ||
              !all(fields %in% names(signal)) ||
              max(abs(rowSums(signal[, fields]) - phi)) > 1e-8 * max(1, abs(phi)))
            stop("Saved signal components or row IDs are inconsistent.", call. = FALSE)
          cz <- sipp_401k_frozen_rank(as.matrix(signal[ord, fields]), row$s)
          if (abs(sum(cz$covariance) - z$covariance[1,1]) > 1e-6 * max(1, z$covariance[1,1]))
            stop("Component covariance does not reconstruct total variance.", call. = FALSE)
          cmeans[[length(cmeans) + 1L]] <- data.frame(outcome_scale = scale,
            profile_id = profiles$profile_id[j], method = method, component = fields, mean = cz$tau)
          ij <- expand.grid(component_1 = fields, component_2 = fields, stringsAsFactors = FALSE)
          ccov[[length(ccov) + 1L]] <- data.frame(outcome_scale = scale,
            profile_id = profiles$profile_id[j], method = method, ij,
            covariance = as.vector(cz$covariance))
        }
        if (row$s != grid$candidates$s[grid$candidates$is_default])
          stop("Saved estimate uses a different default grid order.", call. = FALSE)
        orders <- unique(grid$candidates$s)
        for (s in orders[is.finite(orders)]) {
          tz <- sipp_401k_frozen_rank(phi[ord], s)
          th <- tz$centered[,1]^2 / sum(tz$centered[,1]^2)
          tuning[[length(tuning) + 1L]] <- data.frame(outcome_scale = scale,
            profile_id = profiles$profile_id[j], method = method, s = s, tau = tz$tau,
            se = sqrt(tz$covariance[1,1]), ess = 1 / sum(tz$weights^2),
            radius_90_z = distance[ord][which(cumsum(tz$weights) >= .9)[1L]],
            jk_top1_sq_share = if (any(is.finite(th))) max(th) else NA_real_,
            top_row_id = if (any(is.finite(th))) d$row_id[ord[which.max(th)]] else NA_integer_)
          if (mechanism == "fullsample") {
            expected <- fit$tuning_grids[[scale]]$estimates
            expected <- expected[expected$s == s & expected$target_id == j, ]
            if (any(abs(expected$estimate - tz$tau) > 1e-7) ||
                any(abs(expected$se - sqrt(tz$covariance[1,1])) > 1e-5))
              stop("Stored tuning inference disagrees with frozen reconstruction.", call. = FALSE)
          }
        }
      }
      if (length(centered) > 1L) for (pair in combn(names(centered), 2L, simplify = FALSE)) {
        a <- centered[[pair[1L]]]; b <- centered[[pair[2L]]]
        if (max(abs(a$w - b$w)) > 1e-12)
          stop("Paired comparisons require common final weights.", call. = FALSE)
        covariance <- sum(a$pc * b$pc) / (n * (n - 1L))
        movement <- a$w * (b$phi - a$phi)
        paired[[length(paired) + 1L]] <- data.frame(outcome_scale = scale,
          profile_id = profiles$profile_id[j], method_1 = paste0("dnn_", pair[1L]),
          method_2 = paste0("dnn_", pair[2L]), gap = b$tau - a$tau,
          absolute_signal_movement = sum(abs(movement)),
          treated_gap = sum(movement[d$D == 1]), control_gap = sum(movement[d$D == 0]),
          covariance = covariance, correlation = ratio(covariance, sqrt(a$variance * b$variance)),
          paired_gap_se = sqrt(max(0, a$variance + b$variance - 2 * covariance)))
      }
    }
  }
  bind <- function(x) if (length(x)) do.call(rbind, x) else data.frame()
  list(information = bind(info), covariance_routes = bind(routes), paired = bind(paired),
       tuning = bind(tuning), component_means = bind(cmeans), component_covariance = bind(ccov),
       availability = list(component_decomposition = has_components,
         raw_propensity = has_components),
       scope = "Frozen fitted-array diagnostics; pointwise reference, no coverage or causal-identification certification.")
}
