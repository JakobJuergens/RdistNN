# Structural and contract tests for the specialized plotting panels: plqr/gmm
# autoplot (P-1, G-1), Gram conditioning (P-2), quantile calibration
# (P-3), CATE influence (C-5), TDNN scale decomposition (NP-5), and
# kernel-order sensitivity (NP-6).

.plot_special_plqr_fit <- function(n = 90L, m = 5L, seed = 101L, se = FALSE) {
  demo <- plqr_demo_data(n = n, p = 2L, p_d = 2L, tau = 0.5, seed = seed)
  df <- demo$data
  plqr_dnn(
    X = as.matrix(df[paste0("X", 1:2)]),
    D = as.matrix(df[paste0("D", 1:2)]),
    g = as.numeric(df$g_oracle),
    mD = as.matrix(df[paste0("mD", 1:2)]),
    Y = as.numeric(df$Y),
    theta_eval = demo$truth$theta_true, tau = 0.5,
    X0 = cbind(seq(-0.5, 0.5, length.out = m), 0),
    k = 1L, s = 10L, se = se
  )
}

.plot_special_cate_fit <- function(n = 100L, m = 5L, seed = 102L) {
  set.seed(seed)
  X <- matrix(rnorm(n * 2L), n, 2L,
              dimnames = list(NULL, c("X1", "X2")))
  psi <- rnorm(n, mean = X[, 1L])
  cate_dnn(
    X = X, X0 = cbind(seq(-1, 1, length.out = m), 0),
    pseudo_outcome = psi, k = 1L, s = 15L, se = FALSE,
    return = c("W", "dx")
  )
}

test_that("plqr and gmm autoplot methods carry class-specific dressing", {
  skip_if_not_installed("ggplot2")
  fit_plqr <- .plot_special_plqr_fit(se = TRUE)
  p <- ggplot2::autoplot(fit_plqr)

  # PLQR moments get the zero reference line and the beta/tau subtitle.
  expect_true("GeomHline" %in% .plot_layer_geoms(p))
  expect_match(p$labels$subtitle, "PLQR moments at theta_eval")
  expect_match(p$labels$subtitle, "tau = 0.50", fixed = TRUE)
  expect_identical(class(p$facet)[1L], "FacetWrap")

  dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 100L, kx = 2L, q = 2L,
                     seed = 103L)
  model <- gmm_model_dr_iv(
    z_cols = paste0("Z", seq_len(dgp$q)),
    pi_cols = paste0("piZ", seq_len(dgp$q), "_oracle"),
    muY_col = "muY_oracle",
    muD_col = "muD_oracle"
  )
  fit_gmm <- gmm_dnn(
    data = dgp$data, x_cols = c("X1", "X2"), x0 = dgp$x0,
    model = model, k = 1L, s = 12L
  )
  p_gmm <- ggplot2::autoplot(fit_gmm)
  expect_match(p_gmm$labels$subtitle, "GMM model")
  # Coefficient estimates carry no zero line.
  expect_false("GeomHline" %in% .plot_layer_geoms(p_gmm))
  expect_false("GeomVline" %in% .plot_layer_geoms(p_gmm))

  # plot() methods print invisibly.
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  ret <- withVisible(plot(fit_plqr))
  expect_false(ret$visible)
  expect_s3_class(ret$value, "ggplot")
})

test_that("plot_conditioning facets eigenvalue and log condition number", {
  skip_if_not_installed("ggplot2")
  fit <- .plot_special_plqr_fit()
  dd <- plot_data_diagnostics(fit)
  p <- plot_conditioning(fit)

  expect_identical(levels(p$data$metric),
                   c("smallest eigenvalue", "log10 condition number"))
  # Values reproduce the diagnostics columns (log10 for the condition
  # number facet).
  avail <- dd$gram_available
  expect_equal(
    p$data$value[p$data$metric == "smallest eigenvalue"],
    dd$gram_min_eig[avail]
  )
  expect_equal(
    p$data$value[p$data$metric == "log10 condition number"],
    log10(dd$gram_condition_number[avail])
  )
  # Guides: zero for the eigenvalue facet, log10 cond_limits otherwise.
  hline <- p$layers[[which(.plot_layer_geoms(p) == "GeomHline")[1L]]]
  expect_equal(sort(hline$data$yintercept), sort(c(0, 3, 6)))
  expect_match(p$labels$caption, "ill-conditioned")

  expect_error(plot_conditioning(fit, cond_limits = c(1e6, 1e3)),
               "increasing")
  expect_error(plot_conditioning(.plot_special_cate_fit()), "plqr_result")
})

test_that("plot_calibration reproduces the weighted indicator share", {
  skip_if_not_installed("ggplot2")
  fit <- .plot_special_plqr_fit(m = 5L)
  p <- plot_calibration(fit)

  # Manual recomputation from stored components.
  below <- fit$.model$below_threshold
  share_ref <- as.numeric(crossprod(below, fit$W))
  expect_equal(p$data$share, share_ref)

  # Reference band: tau +/- z sqrt(tau(1-tau)/neff).
  neff <- 1 / colSums(fit$W^2)
  half <- qnorm(0.975) * sqrt(0.5 * 0.5 / neff)
  expect_equal(p$data$band_upper - p$data$band_lower, 2 * half)
  geoms <- .plot_layer_geoms(p)
  expect_true(all(c("GeomRibbon", "GeomHline", "GeomLine") %in% geoms))
  expect_match(p$labels$caption, "reference band")

  # Argument path: residuals reproduce the stored indicator.
  fit_old <- fit
  fit_old$.model$below_threshold <- NULL
  expect_error(plot_calibration(fit_old), "residuals")
  demo <- plqr_demo_data(n = 90L, p = 2L, p_d = 2L, tau = 0.5, seed = 101L)
  df <- demo$data
  r <- as.numeric(df$Y) -
    as.vector(as.matrix(df[paste0("D", 1:2)]) %*% demo$truth$theta_true) -
    as.numeric(df$g_oracle)
  p_arg <- plot_calibration(fit_old, residuals = r)
  expect_equal(p_arg$data$share, share_ref)

  expect_error(plot_calibration(.plot_special_cate_fit()), "plqr_result")
})

test_that("plot_influence computes contributions and labels extremes", {
  skip_if_not_installed("ggplot2")
  fit <- .plot_special_cate_fit(m = 5L)
  p <- plot_influence(fit, targets = 3L)

  # Contributions equal w * (psi - tau_hat) for the shown target.
  psi <- fit$.model$psi
  contrib_ref <- fit$W[, 3L] * (psi - fit$tau[3L])
  expect_equal(sort(p$data$contribution), sort(contrib_ref))
  # Contributions sum to (approximately) zero at the fitted root.
  expect_equal(sum(p$data$contribution), 0, tolerance = 1e-10)

  geoms <- .plot_layer_geoms(p)
  expect_true("GeomText" %in% geoms)
  txt <- p$layers[[which(geoms == "GeomText")[1L]]]$data
  expect_identical(nrow(txt), 3L)
  # Labelled rows are the largest absolute contributions.
  expect_equal(
    sort(abs(txt$contribution)),
    sort(sort(abs(contrib_ref), decreasing = TRUE)[1:3])
  )
  # Distance is stored, so the x axis is distance.
  expect_identical(p$labels$x, "Distance")
  expect_match(p$labels$caption, "top-1 weight share")

  expect_error(plot_influence(.plot_special_plqr_fit()), "cate_result")
  fit_nopsi <- fit
  fit_nopsi$.model$psi <- NULL
  expect_error(plot_influence(fit_nopsi), "pseudo_outcome")
  p_arg <- plot_influence(fit_nopsi, targets = 3L, pseudo_outcome = psi)
  expect_equal(sort(p_arg$data$contribution), sort(contrib_ref))
})

test_that("plot_tdnn_scales overlays refits that match manual estimates", {
  skip_if_not_installed("ggplot2")
  set.seed(104)
  n <- 90L
  x_mat <- matrix(rnorm(n * 2L), n, 2L,
                  dimnames = list(NULL, c("X1", "X2")))
  y_vec <- x_mat[, 1L] + rnorm(n, sd = 0.3)
  X0 <- cbind(seq(-1, 1, length.out = 5L), 0)
  fit <- np_tdnn(X = x_mat, X0 = X0, Y = y_vec, s1 = 5L, s2 = 10L,
                 se = FALSE)
  p <- plot_tdnn_scales(fit)

  expect_identical(levels(p$data$estimator),
                   c("DNN(s = 5)", "DNN(s = 10)", "TDNN(5, 10)"))
  expect_identical(nrow(p$data), 3L * 5L)

  # The overlay curves equal manual refits and the stored TDNN estimate.
  ref_s1 <- np_dnn(X = x_mat, X0 = X0, Y = y_vec, k = 1L, s = 5L,
                   se = FALSE)
  expect_equal(
    p$data$estimate[p$data$estimator == "DNN(s = 5)"],
    as.numeric(ref_s1$mu)
  )
  expect_equal(
    p$data$estimate[p$data$estimator == "TDNN(5, 10)"],
    as.numeric(fit$mu)
  )

  # Wrong class and missing stored data error with remedies.
  expect_error(plot_tdnn_scales(.plot_special_cate_fit()), "tdnn_np_result")
  fit_old <- fit
  fit_old$.model <- NULL
  expect_error(plot_tdnn_scales(fit_old), "supply both")
  p_args <- plot_tdnn_scales(fit_old, X = x_mat, Y = y_vec)
  expect_equal(p_args$data$estimate, p$data$estimate)
})

test_that("plot_tuning evaluates candidate grids and marks the fit", {
  skip_if_not_installed("ggplot2")
  set.seed(105)
  n <- 80L
  x_mat <- matrix(rnorm(n * 2L), n, 2L)
  y_vec <- x_mat[, 1L] + rnorm(n, sd = 0.3)
  fit <- np_dnn(X = x_mat, x0 = c(0, 0), Y = y_vec, k = 1L, se = FALSE)

  p <- plot_tuning(fit, se = FALSE)
  s_fit <- fit$meta$s
  expect_true(s_fit %in% p$data$s)
  expect_identical(sum(p$data$selected), 1L)
  expect_identical(p$data$s[p$data$selected], s_fit)

  # Candidate estimates equal scalar fits without requiring plot_tuning()
  # to call public estimators internally.
  s_other <- setdiff(unique(p$data$s), s_fit)[1L]
  ref <- np_dnn(X = x_mat, x0 = c(0, 0), Y = y_vec, k = 1L, s = s_other,
                se = FALSE)
  expect_equal(p$data$estimate[p$data$s == s_other],
               as.numeric(ref$mu))
  # The fitted point reproduces the original estimate.
  expect_equal(p$data$estimate[p$data$selected], as.numeric(fit$mu))

  # The vline marks the fitted s; caption names provenance.
  geoms <- .plot_layer_geoms(p)
  vl <- p$layers[[which(geoms == "GeomVline")[1L]]]
  expect_identical(vl$data$xintercept, s_fit)
  expect_match(p$labels$caption, "paper_default")

  # cate method, with se: pointranges drawn.
  fit_cate <- .plot_special_cate_fit(m = 2L)
  p_cate <- plot_tuning(fit_cate, se = TRUE)
  expect_true("GeomPointrange" %in% .plot_layer_geoms(p_cate))
  expect_identical(class(p_cate$facet)[1L], "FacetWrap")

  p_plqr <- plot_tuning(.plot_special_plqr_fit(), se = FALSE)
  expect_s3_class(p_plqr, "ggplot")
  expect_identical(class(p_plqr$facet)[1L], "FacetWrap")
  bad_state <- fit
  bad_state$.model$x_mat <- NULL
  expect_error(plot_tuning(bad_state), "required fitted training data.*refit")
})

test_that("PLQR plot quantity controls zero lines and plotted values", {
  skip_if_not_installed("ggplot2")
  fit <- .plot_special_plqr_fit()
  geoms_of <- function(p) vapply(p$layers, function(l) class(l$geom)[1L], character(1L))

  # Default: localized moments with a zero reference line.
  p_moment <- plot_estimates(fit)
  expect_true(any(c("GeomHline", "GeomVline") %in% geoms_of(p_moment)))
  expect_equal(
    sort(unique(p_moment$data$estimate)),
    sort(unique(as.numeric(as.matrix(fit$moment)))),
    tolerance = 1e-12
  )

  # Supplied-theta plots: repeated theta_eval, no automatic zero line.
  p_theta <- plot_estimates(fit, quantity = "theta")
  expect_false(any(c("GeomHline", "GeomVline") %in% geoms_of(p_theta)))
  expect_equal(
    sort(unique(p_theta$data$estimate)),
    sort(unique(as.numeric(as.matrix(fit$theta)))),
    tolerance = 1e-12
  )
})

test_that("candidate-grid PLQR plot data stays in moment coordinates", {
  skip_if_not_installed("ggplot2")
  demo <- plqr_demo_data(n = 70, p = 2, p_d = 2, tau = 0.5, seed = 21)
  df <- demo$data
  fit <- plqr_dnn(
    X = as.matrix(df[paste0("X", 1:2)]),
    D = as.matrix(df[paste0("D", 1:2)]),
    g = as.numeric(df$g_oracle),
    mD = as.matrix(df[paste0("mD", 1:2)]),
    Y = as.numeric(df$Y),
    theta_eval = demo$truth$theta_true, tau = 0.5,
    X0 = rbind(c(0, 0), c(0.4, -0.4)), k = 1L,
    tuning_mode = "grid", se = FALSE,
    return = c("W", "ranks")
  )

  p <- plot_estimates(fit, candidates = "all")
  grid_est <- fit$tuning_grid$estimates
  expect_true(all(p$data$estimate %in% grid_est$estimate))
  # The default candidate rows are the fit's localized moments, not theta.
  default_id <- fit$tuning_grid$candidates$candidate_id[
    fit$tuning_grid$candidates$is_default
  ][1L]
  default_rows <- p$data[p$data$candidate_id == default_id, ]
  expect_equal(
    sort(default_rows$estimate),
    sort(as.numeric(as.matrix(fit$moment))),
    tolerance = 1e-12
  )
})
