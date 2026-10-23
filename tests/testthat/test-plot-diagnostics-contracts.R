# Deeper contract tests for the CATE diagnostic plotting panels. Complements the
# structural checks in test-plot-diagnostics.R with:
#   - built-layer identities for the overlap mass panel (masses recomputed
#     independently from W and t, band/guide geometry, pinned palette),
#   - the user-supplied-treatment path verified numerically,
#   - overlap display-mode boundaries and input validation,
#   - propensity density facet counts, weight mapping, and subsetting,
#   - diagnostics facet values against the extractor and the panel-level
#     restriction of the 0.9n guide,
#   - dashboard panel counts, argument passthrough, single-panel
#     degradation, and a full print smoke test.

.plot_diag_contract_cate_df <- function(n = 120L, seed = 81L) {
  set.seed(seed)
  df <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  df$e <- plogis(0.5 * df$X1)
  df$T <- as.integer(runif(n) < df$e)
  df$m0 <- df$X1
  df$m1 <- df$X1 + 1
  df$Y <- ifelse(df$T == 1L, df$m1, df$m0) + rnorm(n, sd = 0.3)
  df
}

.plot_diag_contract_cate_fit <- function(n = 120L, m = 5L, seed = 81L, se = FALSE,
                          s = 15L) {
  df <- .plot_diag_contract_cate_df(n = n, seed = seed)
  cate_dnn(
    data = df, x_cols = c("X1", "X2"), y_col = "Y", t_col = "T",
    e_col = "e", m0_col = "m0", m1_col = "m1",
    X0 = cbind(seq(-1, 1, length.out = m), 0),
    k = 1L, s = s, se = se
  )
}

.plot_diag_weight_corruptions <- function(W) {
  nonfinite <- W
  nonfinite[1L] <- Inf
  list(observations = W[-1L, , drop = FALSE], targets = W[, -1L, drop = FALSE],
    vector = as.numeric(W), complex = W + 0i, nonfinite = nonfinite)
}

test_that("weight diagnostics and extraction require aligned finite retained weights", {
  cate <- .plot_diag_contract_cate_fit(n = 24L, m = 2L, s = 3L)
  common <- list(X = cate$.model$x_mat, X0 = cate$x0, Y = cate$.model$psi, se = FALSE)
  fits <- list(cate = cate, np = do.call(np_dnn, c(common, list(s = 3L))),
    tdnn = do.call(np_tdnn, c(common, list(s1 = 2L, s2 = 4L))))
  for (fit in fits) {
    expect_no_error(diagnostics(fit))
    expect_identical(nrow(plot_data_weights(fit)), fit$meta$n * fit$meta$m)
    for (weights in .plot_diag_weight_corruptions(fit$W)) {
      bad <- fit
      bad$W <- weights
      before <- serialize(bad, NULL)
      for (route in list(diagnostics, plot_data_weights)) {
        expect_error(route(bad), "object\\$W.*[Rr]efit")
      }
      .expect_identical_fast(serialize(bad, NULL), before)
    }
  }
})

test_that("overlap and influence panels share retained weight validation", {
  skip_if_not_installed("ggplot2")
  df <- .plot_diag_contract_cate_df(n = 24L)
  cate <- .plot_diag_contract_cate_fit(n = 24L, m = 2L, s = 3L)
  clate <- clate_dnn(X = as.matrix(df[c("X1", "X2")]), X0 = cate$x0,
    Y = df$Y, D = df$T, A = df$T, e = df$e, m_y0 = df$m0, m_y1 = df$m1,
    m_d0 = rep(0, nrow(df)), m_d1 = rep(1, nrow(df)), s = 3L, se = FALSE)
  routes <- list(function(x) plot_overlap(x, type = "mass"),
    function(x) plot_overlap(x, type = "propensity"), plot_influence)
  for (fit in list(cate, clate)) for (weights in .plot_diag_weight_corruptions(fit$W)) {
    bad <- fit
    bad$W <- weights
    for (route in routes) expect_error(route(bad), "object\\$W.*[Rr]efit")
  }
})

test_that("overlap mass built layers reproduce independent computations", {
  skip_if_not_installed("ggplot2")
  fit <- .plot_diag_contract_cate_fit(m = 5L)
  limits <- c(0.2, 0.8)
  p <- plot_overlap(fit, type = "mass", share_limits = limits)
  built <- ggplot2::ggplot_build(p)
  geoms <- .plot_layer_geoms(p)

  # Independent mass computation from stored components.
  t_num <- fit$.model$t_vec
  treated_ref <- as.numeric(crossprod(t_num, fit$W))
  control_ref <- as.numeric(crossprod(1 - t_num, fit$W))

  line_df <- built$data[[which(geoms == "GeomLine")[1L]]]
  cols <- split(line_df$y, line_df$colour)
  expect_setequal(names(cols), c("#0072B2", "#D55E00"))
  # Okabe-Ito assignment: treated is blue, control vermillion; lines are
  # ordered along the profile axis, which is increasing in X1 here.
  expect_equal(sort(cols[["#0072B2"]]), sort(treated_ref))
  expect_equal(sort(cols[["#D55E00"]]), sort(control_ref))

  # Warning-band rectangles sit exactly at the limits.
  rect_idx <- which(geoms == "GeomRect")
  expect_equal(built$data[[rect_idx[1L]]]$ymax, limits[1L])
  expect_equal(built$data[[rect_idx[2L]]]$ymin, limits[2L])
  hline <- p$layers[[which(geoms == "GeomHline")[1L]]]
  expect_equal(sort(hline$data$yintercept), limits)
})

test_that("user-supplied treatment reproduces the stored-treatment panel", {
  skip_if_not_installed("ggplot2")
  df <- .plot_diag_contract_cate_df(seed = 82L)
  fit_data <- cate_dnn(
    data = df, x_cols = c("X1", "X2"), y_col = "Y", t_col = "T",
    e_col = "e", m0_col = "m0", m1_col = "m1",
    X0 = cbind(seq(-1, 1, length.out = 5L), 0),
    k = 1L, s = 15L, se = FALSE
  )
  fit_mat <- cate_dnn(
    X = as.matrix(df[, c("X1", "X2")]),
    X0 = cbind(seq(-1, 1, length.out = 5L), 0),
    pseudo_outcome = rnorm(nrow(df)),
    k = 1L, s = 15L, se = FALSE
  )

  p_stored <- plot_overlap(fit_data, type = "mass")
  p_arg <- plot_overlap(fit_mat, type = "mass", treatment = df$T)
  # Same X and targets: identical weights, hence identical masses.
  expect_equal(p_arg$data$mass, p_stored$data$mass)
  expect_identical(p_arg$data$arm, p_stored$data$arm)
})

test_that("overlap display mode and inputs are validated at the boundaries", {
  skip_if_not_installed("ggplot2")
  # Profile lines at m = 4, dodged columns at m = 3.
  fit4 <- .plot_diag_contract_cate_fit(m = 4L, seed = 83L)
  expect_true("GeomLine" %in% .plot_layer_geoms(plot_overlap(fit4, type = "mass")))
  fit3 <- .plot_diag_contract_cate_fit(m = 3L, seed = 84L)
  expect_true("GeomCol" %in% .plot_layer_geoms(plot_overlap(fit3, type = "mass")))

  # share_limits validation.
  expect_error(plot_overlap(fit4, share_limits = c(0.9, 0.1)), "increasing")
  expect_error(plot_overlap(fit4, share_limits = 0.5), "increasing|two")
  expect_error(plot_overlap(fit4, share_limits = c(0.1, Inf)), "finite")

  # Overlap panels need the stored weight matrix.
  df <- .plot_diag_contract_cate_df(seed = 85L)
  fit_now <- cate_dnn(
    data = df, x_cols = c("X1", "X2"), y_col = "Y", t_col = "T",
    e_col = "e", m0_col = "m0", m1_col = "m1",
    x0 = c(0, 0), k = 1L, s = 15L, se = FALSE, return = c("ess")
  )
  expect_error(plot_overlap(fit_now, type = "mass"), class = "rdistnn_extra_unavailable")

  # Treatment argument must be binary.
  fit_mat <- cate_dnn(
    X = as.matrix(df[, c("X1", "X2")]), x0 = c(0, 0),
    pseudo_outcome = rnorm(nrow(df)), k = 1L, s = 15L, se = FALSE
  )
  expect_error(
    plot_overlap(fit_mat, type = "mass", treatment = runif(nrow(df))),
    "0|binary"
  )
})

test_that("propensity densities map weights and respect target subsetting", {
  skip_if_not_installed("ggplot2")
  fit <- .plot_diag_contract_cate_fit(m = 10L, seed = 86L)

  # Default subsetting above eight targets: at most three facets.
  p <- plot_overlap(fit, type = "propensity")
  built <- ggplot2::ggplot_build(p)
  expect_lte(length(unique(built$layout$layout$PANEL)), 3L)

  # Explicit targets: exactly those facets, in order.
  p2 <- plot_overlap(fit, type = "propensity", targets = c(2L, 7L))
  expect_identical(levels(p2$data$target), c("target 2", "target 7"))
  expect_identical(nrow(p2$data), 2L * fit$meta$n)

  # The density layer maps the per-target normalized weights.
  dens_layer <- p2$layers[[which(.plot_layer_geoms(p2) == "GeomDensity")[1L]]]
  expect_true("weight" %in% names(dens_layer$mapping))
  w2 <- p2$data$w[p2$data$target_id == 2L]
  expect_equal(sum(w2), 1)
  expect_equal(w2, fit$W[, 2L] / sum(fit$W[, 2L]))

  # Missing retained AIPW propensities require refitting for either view.
  fit_no_e <- fit
  fit_no_e$.model$e <- NULL
  expect_error(plot_overlap(fit_no_e, type = "propensity"), "object\\$\\.model\\$e.*[Rr]efit")
  expect_error(plot_overlap(fit_no_e, type = "mass"), "object\\$\\.model\\$e.*[Rr]efit")
})

test_that("diagnostics facets reproduce the extractor and isolate the guide", {
  skip_if_not_installed("ggplot2")
  fit <- .plot_diag_contract_cate_fit(m = 5L, seed = 87L)
  dd <- plot_data_diagnostics(fit)
  p <- ggplot2::autoplot(diagnostics(fit))
  built <- ggplot2::ggplot_build(p)
  geoms <- .plot_layer_geoms(p)

  # Facet values equal the extractor columns, target by target.
  long <- p$data
  expect_equal(long$value[long$metric == "n_eff"], dd$n_eff)
  expect_equal(long$value[long$metric == "max weight"], dd$max_weight)
  expect_equal(long$value[long$metric == "top-5 weight share"],
               dd$top5_weight_share)

  # The 0.9n guide is drawn in exactly one panel: the n_eff facet.
  hline_built <- built$data[[which(geoms == "GeomHline")[1L]]]
  expect_identical(length(unique(hline_built$PANEL)), 1L)
  layout <- built$layout$layout
  guide_panel <- layout$metric[layout$PANEL == unique(hline_built$PANEL)]
  expect_identical(as.character(guide_panel), "n_eff")

  # Works across estimator classes, including TDNN via np inheritance.
  set.seed(88)
  n <- 60L
  Xm <- matrix(rnorm(n * 2L), n, 2L)
  fit_tdnn <- np_tdnn(X = Xm, x0 = c(0, 0), Y = rnorm(n),
                      s1 = 4L, s2 = 8L, se = FALSE)
  p_tdnn <- ggplot2::autoplot(diagnostics(fit_tdnn))
  expect_s3_class(p_tdnn, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p_tdnn))
})

test_that("dashboards carry the expected panel counts and pass arguments", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  # Full-input cate dashboard: four panels (estimates, overlap,
  # propensity, weights).
  fit <- .plot_diag_contract_cate_fit(m = 5L, seed = 89L, se = TRUE)
  dash <- plot_dashboard(fit)
  expect_identical(length(dash$patches$plots) + 1L, 4L)

  # Matrix-interface fit: overlap and propensity unavailable, diagnostics
  # substituted -> three panels.
  df <- .plot_diag_contract_cate_df(seed = 90L)
  fit_mat <- cate_dnn(
    X = as.matrix(df[, c("X1", "X2")]),
    X0 = cbind(seq(-1, 1, length.out = 5L), 0),
    pseudo_outcome = rnorm(nrow(df)), k = 1L, s = 15L, se = FALSE
  )
  dash_mat <- plot_dashboard(fit_mat)
  expect_identical(length(dash_mat$patches$plots) + 1L, 3L)

  # Arguments pass through to the estimates panel (checked through the
  # degraded list, which returns the panels directly).
  testthat::local_mocked_bindings(
    .plot_pkg_installed = function(pkg) pkg != "patchwork",
    .package = "RdistNN"
  )
  suppressMessages(panels <- plot_dashboard(fit, level = 0.8))
  expect_identical(attr(panels$estimates$data, "level"), 0.8)
  expect_identical(names(panels),
                   c("estimates", "overlap", "propensity", "weights"))
})

test_that("dashboards degrade to a single panel and survive printing", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  # np fit without stored W: weight and diagnostics panels are dropped,
  # leaving a single-panel dashboard that still composes.
  set.seed(91)
  n <- 60L
  X <- matrix(rnorm(n * 2L), n, 2L)
  fit_now <- np_dnn(
    X = X, X0 = cbind(seq(-1, 1, length.out = 5L), 0),
    Y = X[, 1L] + rnorm(n, sd = 0.3), k = 1L, s = 10L,
    se = FALSE, return = c("ess")
  )
  dash <- plot_dashboard(fit_now)
  expect_s3_class(dash, "patchwork")
  expect_identical(length(dash$patches$plots) + 1L, 1L)

  # Full render smoke test: the composed dashboard prints without error.
  fit <- .plot_diag_contract_cate_fit(m = 5L, seed = 92L, se = TRUE)
  dash_full <- plot_dashboard(fit)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(print(dash_full))

  # paper_default provenance reaches the annotation.
  df <- .plot_diag_contract_cate_df(seed = 93L)
  fit_def <- cate_dnn(
    data = df, x_cols = c("X1", "X2"), y_col = "Y", t_col = "T",
    e_col = "e", m0_col = "m0", m1_col = "m1",
    X0 = cbind(seq(-1, 1, length.out = 5L), 0),
    k = 1L, se = FALSE
  )
  dash_def <- plot_dashboard(fit_def)
  expect_match(dash_def$patches$annotation$caption, "paper_default")
  expect_match(dash_def$patches$annotation$caption, "jackknife = none",
               fixed = TRUE)
})
