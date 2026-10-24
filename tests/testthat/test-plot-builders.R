# Structural tests for the Phase 1 plot builders (plan panels NP-1, C-1,
# C-2, NP-2, NP-3) and the autoplot()/plot() methods. No rendering: tests
# assert on ggplot object structure, layer geoms, mapped data, and labels.

.pb_np_profile_fit <- function(n = 80L, m = 6L, se = TRUE,
                               return = c("W"), seed = 51L,
                               k = 1L, s = 12L) {
  .test_rng_scope({
    x_mat <- matrix(rnorm(n * 2L), n, 2L,
                    dimnames = list(NULL, c("X1", "X2")))
    y_vec <- x_mat[, 1L] + rnorm(n, sd = 0.3)
    X0 <- cbind(seq(-1, 1, length.out = m), 0)
    np_dnn(X = x_mat, X0 = X0, Y = y_vec, k = k, s = s,
           se = se, return = return)
  }, seed = seed)
}

.pb_cate_fit <- function(n = 80L, X0 = rbind(c(0, 0), c(1, 1), c(-1, 2)),
                         se = TRUE, seed = 52L) {
  .test_rng_scope({
    X <- matrix(rnorm(n * 2L), n, 2L,
                dimnames = list(NULL, c("X1", "X2")))
    psi <- rnorm(n, mean = X[, 1L])
    cate_dnn(X = X, X0 = X0, pseudo_outcome = psi, k = 1L, s = 12L, se = se)
  }, seed = seed)
}

test_that("missing ggplot2 produces an actionable guard error", {
  fit <- .pb_np_profile_fit(se = FALSE)
  testthat::local_mocked_bindings(
    .plot_pkg_installed = function(pkg) FALSE,
    .package = "RdistNN"
  )
  expect_error(plot_estimates(fit), "ggplot2")
  expect_error(plot_weights(fit), "install.packages", fixed = TRUE)
})

test_that("profile mode draws ribbon + line with the profile axis", {
  skip_if_not_installed("ggplot2")
  fit <- .pb_np_profile_fit(se = TRUE)
  p <- plot_estimates(fit)

  expect_s3_class(p, "ggplot")
  geoms <- .plot_layer_geoms(p)
  expect_true("GeomRibbon" %in% geoms)
  expect_true("GeomLine" %in% geoms)
  expect_false("GeomHline" %in% geoms) # no zero line for NP means
  expect_identical(p$labels$x, "X1")
  expect_identical(p$labels$y, "Estimate")
  expect_match(p$labels$caption, "n_eff")
  expect_match(p$labels$caption, "s = 12.*s/n =.*\\(manual\\)")
})

test_that("se = FALSE drops the band and says so; x_var overrides the axis", {
  skip_if_not_installed("ggplot2")
  fit <- .pb_np_profile_fit(se = FALSE)
  p <- plot_estimates(fit)
  geoms <- .plot_layer_geoms(p)
  expect_false("GeomRibbon" %in% geoms)
  expect_match(p$labels$caption, "no standard errors stored")

  p2 <- plot_estimates(fit, x_var = "X2")
  expect_identical(p2$labels$x, "X2")
  p3 <- plot_estimates(fit, x_var = 2)
  expect_identical(p3$labels$x, "X2")
  expect_error(plot_estimates(fit, x_var = "nope"), "x_var")
})

test_that("cate panels draw the zero line; small target sets go forest", {
  skip_if_not_installed("ggplot2")
  fit <- .pb_cate_fit()
  p <- plot_estimates(fit)
  geoms <- .plot_layer_geoms(p)

  # Non-profile with 3 targets: forest layout, zero vline, pointranges.
  expect_true("GeomVline" %in% geoms)
  expect_true("GeomSegment" %in% geoms)
  expect_identical(levels(p$data$target_label)[1L], "X1 = -1, X2 = 2")

  # Explicit pointrange style over target ids uses hline instead.
  p2 <- plot_estimates(fit, style = "pointrange")
  geoms2 <- .plot_layer_geoms(p2)
  expect_true("GeomHline" %in% geoms2)
  expect_false("GeomVline" %in% geoms2)

  # Profile style without a profile axis errors.
  expect_error(plot_estimates(fit, style = "profile"), "profile axis")
})

test_that("shade_significant adds a rug only where CIs exclude zero", {
  .test_rng_scope({
    skip_if_not_installed("ggplot2")
    set.seed(53)
    n <- 100L
    X <- matrix(rnorm(n * 2L), n, 2L)
    psi <- rnorm(n, mean = 4 * X[, 1L], sd = 0.5)
    X0 <- cbind(seq(-1.5, 1.5, length.out = 8L), 0)
    fit <- cate_dnn(X = X, X0 = X0, pseudo_outcome = psi, k = 1L, s = 15L,
                    se = TRUE)

    p0 <- plot_estimates(fit)
    expect_false("GeomRug" %in% .plot_layer_geoms(p0))
    p1 <- plot_estimates(fit, shade_significant = TRUE)
    geoms <- .plot_layer_geoms(p1)
    if ("GeomRug" %in% geoms) {
      rug_layer <- p1$layers[[which(geoms == "GeomRug")[1L]]]
      sig <- rug_layer$data
      expect_true(all(sign(sig$ci_lower) == sign(sig$ci_upper)))
    }
    # Zero hline present in profile mode for cate.
    expect_true("GeomHline" %in% geoms)
  })
})

test_that("autoplot and plot methods delegate to plot_estimates", {
  skip_if_not_installed("ggplot2")
  fit_np <- .pb_np_profile_fit(se = FALSE)
  fit_cate <- .pb_cate_fit(se = FALSE)

  p_np <- ggplot2::autoplot(fit_np)
  expect_s3_class(p_np, "ggplot")
  p_cate <- ggplot2::autoplot(fit_cate, style = "pointrange")
  expect_s3_class(p_cate, "ggplot")

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  ret <- withVisible(plot(fit_np))
  expect_false(ret$visible)
  expect_s3_class(ret$value, "ggplot")
})

test_that("weight profile and cumulative views carry guides and captions", {
  skip_if_not_installed("ggplot2")
  n <- 80L
  fit <- .pb_np_profile_fit(n = n, m = 2L, se = FALSE,
                            return = c("W", "dx"))
  p <- plot_weights(fit, type = "profile")
  geoms <- .plot_layer_geoms(p)
  expect_true("GeomStep" %in% geoms)
  vl <- p$layers[[which(geoms == "GeomVline")[1L]]]
  expect_identical(vl$data$xintercept, n - 12L + 1L)
  expect_identical(p$labels$x, "Rank")

  p_log <- plot_weights(fit, type = "profile", log_weights = TRUE)
  expect_match(p_log$labels$caption, "non-positive weights dropped")
  log_layer_data <- p_log$data
  expect_true(all(log_layer_data$weight > 0))

  p_cum <- plot_weights(fit, type = "cumulative")
  geoms_cum <- .plot_layer_geoms(p_cum)
  expect_true("GeomLine" %in% geoms_cum)
  expect_true("GeomHline" %in% geoms_cum)
  expect_true("GeomPoint" %in% geoms_cum)
  expect_identical(p_cum$labels$y, "Cumulative weight share")
})

test_that("weight profile guide uses exact DNN support edge for k > 1", {
  skip_if_not_installed("ggplot2")
  n <- 80L
  k <- 3L
  s <- 12L
  fit <- .pb_np_profile_fit(
    n = n, m = 2L, se = FALSE, return = c("W", "dx"), k = k, s = s
  )
  p <- plot_weights(fit, type = "profile")
  geoms <- .plot_layer_geoms(p)
  vl <- p$layers[[which(geoms == "GeomVline")[1L]]]

  expect_identical(vl$data$xintercept, n - s + k)
  expect_false(identical(vl$data$xintercept, n - s + 1L))
})

test_that("kernel view requires distances and marks the rank-s radius", {
  skip_if_not_installed("ggplot2")
  fit_nodx <- .pb_np_profile_fit(m = 1L, se = FALSE, return = c("W"))
  expect_error(plot_weights(fit_nodx, type = "kernel"), "dx")

  fit <- .pb_np_profile_fit(m = 1L, se = FALSE, return = c("W", "dx"))
  p <- plot_weights(fit, type = "kernel")
  geoms <- .plot_layer_geoms(p)
  expect_true("GeomStep" %in% geoms)
  expect_identical(p$labels$x, "Distance")
  vl <- p$layers[[which(geoms == "GeomVline")[1L]]]
  dw <- plot_data_weights(fit)
  expect_identical(vl$data$xintercept, dw$distance[fit$meta$s])
})

test_that("TDNN weight panels add the zero line; many targets are subset", {
  .test_rng_scope({
    skip_if_not_installed("ggplot2")
    set.seed(54)
    n <- 70L
    Xm <- matrix(rnorm(n * 2L), n, 2L)
    fit_tdnn <- np_tdnn(X = Xm, x0 = c(0, 0), Y = rnorm(n),
                        s1 = 4L, s2 = 8L, se = FALSE,
                        return = c("W", "dx"))
    p <- plot_weights(fit_tdnn, type = "kernel")
    geoms <- .plot_layer_geoms(p)
    expect_true("GeomHline" %in% geoms)
    expect_match(p$labels$caption, "signed")

    fit_many <- .pb_np_profile_fit(m = 10L, se = FALSE, return = c("W"))
    p_many <- plot_weights(fit_many, type = "profile")
    expect_lte(length(unique(p_many$data$target_id)), 3L)
    p_pick <- plot_weights(fit_many, type = "profile", targets = c(2L, 4L))
    expect_identical(sort(unique(p_pick$data$target_id)), c(2L, 4L))
    expect_error(plot_weights(fit_many, targets = 99L), "target ids")
  })
})
