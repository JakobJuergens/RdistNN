# Structural and contract tests for the CATE diagnostic plotting panels: overlap
# (C-3/C-4), diagnostics facets (NP-4), dashboards (C-6 and the np
# analogue), and the meta$jackknife plumbing behind the dashboard
# annotation.

.plot_diag_cate_df <- function(n = 120L, seed = 71L) {
  .test_rng_scope({
    df <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
    df$e <- plogis(0.5 * df$X1)
    df$T <- as.integer(runif(n) < df$e)
    df$m0 <- df$X1
    df$m1 <- df$X1 + 1
    df$Y <- ifelse(df$T == 1L, df$m1, df$m0) + rnorm(n, sd = 0.3)
    df
  }, seed = seed)
}

.plot_diag_cate_fit <- function(n = 120L, m = 5L, seed = 71L, se = FALSE) {
  df <- .plot_diag_cate_df(n = n, seed = seed)
  cate_dnn(
    data = df, x_cols = c("X1", "X2"), y_col = "Y", t_col = "T",
    e_col = "e", m0_col = "m0", m1_col = "m1",
    X0 = cbind(seq(-1, 1, length.out = m), 0),
    k = 1L, s = 15L, se = se
  )
}

.plot_diag_matrix_cate_fit <- function(n = 100L, m = 5L, seed = 72L) {
  .test_rng_scope({
    X <- matrix(rnorm(n * 2L), n, 2L)
    cate_dnn(
      X = X, X0 = cbind(seq(-1, 1, length.out = m), 0),
      pseudo_outcome = rnorm(n, mean = X[, 1L]),
      k = 1L, s = 15L, se = FALSE
    )
  }, seed = seed)
}

test_that("meta stores the jackknife type across wrappers", {
  .test_rng_scope({
    fit_se <- .plot_diag_cate_fit(se = TRUE)
    expect_identical(fit_se$meta$jackknife, "fast")
    fit_nose <- .plot_diag_cate_fit(se = FALSE)
    expect_identical(fit_nose$meta$jackknife, "none")

    set.seed(73)
    n <- 50L
    X <- matrix(rnorm(n * 2L), n, 2L)
    y <- rnorm(n)
    fit_np <- np_dnn(X = X, x0 = c(0, 0), Y = y, k = 1L, s = 8L,
                     se = TRUE, jackknife = "loo")
    expect_identical(fit_np$meta$jackknife, "loo")
    fit_tdnn <- np_tdnn(X = X, x0 = c(0, 0), Y = y, s1 = 4L, s2 = 8L,
                        se = FALSE)
    expect_identical(fit_tdnn$meta$jackknife, "none")
  })
})

test_that("overlap mass panel draws arms, guides, bands, and violations", {
  skip_if_not_installed("ggplot2")
  fit <- .plot_diag_cate_fit(m = 5L)
  p <- plot_overlap(fit, type = "mass")
  geoms <- .plot_layer_geoms(p)

  # Two warning-band rects, dotted share guides, and profile lines.
  expect_identical(sum(geoms == "GeomRect"), 2L)
  expect_true("GeomHline" %in% geoms)
  expect_true("GeomLine" %in% geoms)
  expect_identical(p$labels$x, "X1")

  # Long data: both arms per target; masses partition unit mass.
  expect_identical(nrow(p$data), 2L * 5L)
  mass_sum <- tapply(p$data$mass, p$data$target_id, sum)
  expect_equal(as.numeric(mass_sum), rep(1, 5L), tolerance = 1e-10)

  # Masses agree with the diagnostics frame.
  dd <- plot_data_diagnostics(fit)
  treated <- p$data[p$data$arm == "treated", ]
  expect_equal(treated$mass[order(treated$target_id)], dd$treated_mass)

  # Tight limits force violations into the caption.
  p_tight <- plot_overlap(fit, type = "mass",
                          share_limits = c(0.49, 0.51))
  expect_match(p_tight$labels$caption, "treated share outside")

  # Few targets without a profile axis: dodged columns.
  fit2 <- .plot_diag_cate_fit(m = 2L, seed = 74L)
  p2 <- plot_overlap(fit2, type = "mass")
  expect_true("GeomCol" %in% .plot_layer_geoms(p2))
})

test_that("overlap inputs follow the precedence chain", {
  .test_rng_scope({
    skip_if_not_installed("ggplot2")
    fit_mat <- .plot_diag_matrix_cate_fit()

    # No stored treatment: actionable error naming both remedies.
    expect_error(plot_overlap(fit_mat, type = "mass"), "t_col")
    expect_error(plot_overlap(fit_mat, type = "mass"), "treatment")

    # User-supplied treatment vector unlocks the panel.
    set.seed(72)
    X <- matrix(rnorm(100L * 2L), 100L, 2L) # same draw as the fit helper
    t_vec <- as.integer(runif(100L) < 0.5)
    p <- plot_overlap(fit_mat, type = "mass", treatment = t_vec)
    expect_s3_class(p, "ggplot")
    expect_error(
      plot_overlap(fit_mat, type = "mass", treatment = t_vec[-1L]),
      "length"
    )

    # Propensity view: stored e gives densities; matrix fit errors with the
    # e_col remedy; a supplied propensity vector unlocks densities.
    fit <- .plot_diag_cate_fit(m = 5L)
    p_dens <- plot_overlap(fit, type = "propensity", targets = c(1L, 5L))
    expect_true("GeomDensity" %in% .plot_layer_geoms(p_dens))
    expect_identical(class(p_dens$facet)[1L], "FacetWrap")

    expect_error(plot_overlap(fit_mat, type = "propensity"), "e_col")
    p_dens2 <- plot_overlap(fit_mat, type = "propensity",
                            propensity = plogis(rnorm(100L)))
    expect_true("GeomDensity" %in% .plot_layer_geoms(p_dens2))
    expect_error(
      plot_overlap(fit_mat, type = "propensity", propensity = rep(2, 100L)),
      "in \\(0, 1\\)"
    )

    # Non-cate objects are rejected.
    set.seed(75)
    fit_np <- np_dnn(X = X, x0 = c(0, 0), Y = rnorm(100L), k = 1L, s = 8L,
                     se = FALSE)
    expect_error(plot_overlap(fit_np), "rdistnn_cate_result")
  })
})

test_that("diagnostics autoplot facets three metrics with the 0.9n guide", {
  skip_if_not_installed("ggplot2")
  fit <- .plot_diag_cate_fit(m = 5L)
  d <- diagnostics(fit)
  p <- ggplot2::autoplot(d)

  expect_s3_class(p, "ggplot")
  expect_identical(class(p$facet)[1L], "FacetWrap")
  expect_identical(levels(p$data$metric),
                   c("n_eff", "max weight", "top-5 weight share"))
  expect_identical(nrow(p$data), 3L * 5L)

  geoms <- .plot_layer_geoms(p)
  hline <- p$layers[[which(geoms == "GeomHline")[1L]]]
  expect_equal(hline$data$yintercept, 0.9 * fit$meta$n)
  expect_identical(as.character(hline$data$metric), "n_eff")
  expect_match(p$labels$caption, "0.9 n")

  # plot() prints and returns invisibly.
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  ret <- withVisible(plot(d))
  expect_false(ret$visible)
  expect_s3_class(ret$value, "ggplot")
})

test_that("cate dashboard composes with the metadata annotation contract", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  fit <- .plot_diag_cate_fit(m = 5L, se = TRUE)
  dash <- plot_dashboard(fit)
  expect_s3_class(dash, "patchwork")

  # Annotation contract: caption carries the meta fields verbatim.
  cap <- dash$patches$annotation$caption
  meta <- fit$meta
  expect_match(cap, sprintf("n = %d", meta$n), fixed = TRUE)
  expect_match(cap, sprintf("m = %d", meta$m), fixed = TRUE)
  expect_match(cap, sprintf("k = %d", meta$k), fixed = TRUE)
  expect_match(cap, sprintf("s = %d (manual)", meta$s), fixed = TRUE)
  expect_match(cap, sprintf("distance = %s", meta$distance), fixed = TRUE)
  expect_match(cap, sprintf("jackknife = %s", meta$jackknife), fixed = TRUE)
})

test_that("dashboards substitute and degrade per availability", {
  .test_rng_scope({
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("patchwork")

    # Matrix-interface cate fit: overlap unavailable, diagnostics facets
    # substituted; the dashboard still composes.
    fit_mat <- .plot_diag_matrix_cate_fit()
    dash <- plot_dashboard(fit_mat)
    expect_s3_class(dash, "patchwork")

    # np dashboard composes from estimates, diagnostics, and weight views.
    set.seed(76)
    n <- 80L
    X <- matrix(rnorm(n * 2L), n, 2L)
    fit_np <- np_dnn(
      X = X, X0 = cbind(seq(-1, 1, length.out = 5L), 0),
      Y = X[, 1L] + rnorm(n, sd = 0.3), k = 1L, s = 12L,
      return = c("W", "dx")
    )
    dash_np <- plot_dashboard(fit_np)
    expect_s3_class(dash_np, "patchwork")

    # Unsupported classes are rejected.
    expect_error(plot_dashboard(lm(y ~ x, data.frame(y = 1:4, x = 1:4))),
                 "not implemented")
  })
})

test_that("dashboard returns a named panel list when patchwork is absent", {
  skip_if_not_installed("ggplot2")
  fit <- .plot_diag_cate_fit(m = 5L)
  testthat::local_mocked_bindings(
    .plot_pkg_installed = function(pkg) pkg != "patchwork",
    .package = "RdistNN"
  )
  expect_message(
    panels <- plot_dashboard(fit),
    "patchwork is not installed"
  )
  expect_true(is.list(panels))
  expect_true(all(vapply(panels, inherits, logical(1L), what = "ggplot")))
  expect_true(all(c("estimates", "overlap") %in% names(panels)))
})
