# Deeper contract tests for the Phase 1 plot builders. Complements the
# structural checks in test-plot-builders.R with:
#   - full ggplot_build() smoke coverage over the mode/type matrix,
#   - built-layer data identities against the plot_data_* extractors,
#   - level propagation into ribbon bounds,
#   - display-mode boundary rules (m = 3/4, m = 12/13),
#   - facet behavior for multi-coefficient (plqr/gmm) fits,
#   - forest ordering and x_var-forced profiles,
#   - rug/significance consistency with the extractor,
#   - kernel guides for multiple targets and class-agnostic weight panels,
#   - guard ordering when both ggplot2 and W are missing.

.pbc_np_fit <- function(n = 80L, m = 6L, se = TRUE, return = c("W"),
                        seed = 61L, s = 12L) {
  set.seed(seed)
  x_mat <- matrix(rnorm(n * 2L), n, 2L,
                  dimnames = list(NULL, c("X1", "X2")))
  y_vec <- x_mat[, 1L] + rnorm(n, sd = 0.3)
  X0 <- cbind(seq(-1, 1, length.out = m), 0)
  np_dnn(X = x_mat, X0 = X0, Y = y_vec, k = 1L, s = s, se = se,
         return = return)
}

.pbc_cate_targets <- function(m, seed = 62L) {
  # Non-profile target sets: both coordinates vary.
  set.seed(seed)
  cbind(seq(-1, 1, length.out = m), seq(1, -1, length.out = m))
}

.pbc_cate_fit <- function(n = 90L, m = 3L, se = TRUE, seed = 63L,
                          psi_scale = 1) {
  set.seed(seed)
  X <- matrix(rnorm(n * 2L), n, 2L,
              dimnames = list(NULL, c("X1", "X2")))
  psi <- rnorm(n, mean = psi_scale * X[, 1L])
  cate_dnn(X = X, X0 = .pbc_cate_targets(m), pseudo_outcome = psi,
           k = 1L, s = 12L, se = se)
}

test_that("every estimate mode and weight view survives ggplot_build", {
  skip_if_not_installed("ggplot2")
  fit_prof <- .pbc_np_fit(m = 6L, se = TRUE, return = c("W", "dx"))
  fit_prof_nose <- .pbc_np_fit(m = 6L, se = FALSE, return = c("W", "dx"))
  fit_small <- .pbc_cate_fit(m = 3L, se = TRUE)

  plots <- list(
    plot_estimates(fit_prof),
    plot_estimates(fit_prof, style = "pointrange"),
    plot_estimates(fit_prof_nose),
    plot_estimates(fit_small),
    plot_estimates(fit_small, style = "pointrange"),
    plot_estimates(fit_small, se = TRUE, style = "forest"),
    plot_weights(fit_prof, type = "profile"),
    plot_weights(fit_prof, type = "profile", log_weights = TRUE),
    plot_weights(fit_prof, type = "cumulative"),
    plot_weights(fit_prof, type = "kernel"),
    plot_weights(fit_prof, type = "profile", targets = 2L)
  )
  for (p in plots) {
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
  }
})

test_that("built ribbon and line layers reproduce the extractor values", {
  skip_if_not_installed("ggplot2")
  fit <- .pbc_np_fit(m = 6L, se = TRUE)
  level <- 0.8
  p <- plot_estimates(fit, level = level)
  df <- plot_data_estimates(fit, level = level)
  built <- ggplot2::ggplot_build(p)

  geoms <- vapply(p$layers, function(l) class(l$geom)[1L], character(1L))
  ribbon <- built$data[[which(geoms == "GeomRibbon")[1L]]]
  line <- built$data[[which(geoms == "GeomLine")[1L]]]

  # X0 is increasing in X1, so built rows align with extractor rows.
  expect_equal(ribbon$ymin, df$ci_lower)
  expect_equal(ribbon$ymax, df$ci_upper)
  expect_equal(ribbon$x, df$x0_1)
  expect_equal(line$y, df$estimate)

  # The plot-level data is the extractor frame itself.
  expect_identical(p$data$estimate, df$estimate)
  expect_identical(attr(p$data, "level"), level)
})

test_that("display-mode boundaries follow the documented rules", {
  skip_if_not_installed("ggplot2")
  geoms_of <- function(p) {
    vapply(p$layers, function(l) class(l$geom)[1L], character(1L))
  }

  # Profile axis with m = 4: profile band; with m = 3: pointrange.
  fit4 <- .pbc_np_fit(m = 4L, se = TRUE)
  expect_true("GeomRibbon" %in% geoms_of(plot_estimates(fit4)))
  fit3 <- .pbc_np_fit(m = 3L, se = TRUE)
  g3 <- geoms_of(plot_estimates(fit3))
  expect_false("GeomRibbon" %in% g3)
  expect_true("GeomPointrange" %in% g3)

  # Non-profile targets: forest through m = 12, pointrange at m = 13.
  fit12 <- .pbc_cate_fit(m = 12L, se = FALSE)
  p12 <- plot_estimates(fit12)
  expect_true("target_label" %in% names(p12$data))
  fit13 <- .pbc_cate_fit(m = 13L, se = FALSE)
  p13 <- plot_estimates(fit13)
  expect_false("target_label" %in% names(p13$data))
  expect_identical(p13$labels$x, "target")

  # m = 1: single pointrange (or point without se).
  fit1 <- .pbc_np_fit(m = 1L, se = TRUE)
  g1 <- geoms_of(plot_estimates(fit1))
  expect_true("GeomPointrange" %in% g1)
})

test_that("forest ordering puts the first target on top", {
  skip_if_not_installed("ggplot2")
  fit <- .pbc_cate_fit(m = 3L, se = TRUE)
  p <- plot_estimates(fit, style = "forest")
  lv <- levels(p$data$target_label)
  labs_map <- p$data$target_label[order(p$data$target_id)]
  # Reversed levels: the first target's label is the last (top) level.
  expect_identical(lv[length(lv)], as.character(labs_map[1L]))
})

test_that("x_var forces a profile when several coordinates vary", {
  skip_if_not_installed("ggplot2")
  fit <- .pbc_cate_fit(m = 6L, se = TRUE)

  # Auto mode: no unique varying coordinate -> forest (m <= 12).
  p_auto <- plot_estimates(fit)
  expect_true("target_label" %in% names(p_auto$data))

  # x_var override: profile over the named coordinate.
  p_x1 <- plot_estimates(fit, x_var = "X1")
  geoms <- vapply(p_x1$layers, function(l) class(l$geom)[1L], character(1L))
  expect_true("GeomLine" %in% geoms)
  expect_identical(p_x1$labels$x, "X1")
  # CATE keeps its zero hline in the forced-profile mode.
  expect_true("GeomHline" %in% geoms)
})

test_that("multi-coefficient fits facet by coefficient with free scales", {
  skip_if_not_installed("ggplot2")
  demo <- plqr_demo_data(n = 80L, p = 2L, p_d = 2L, tau = 0.5, seed = 12)
  df_p <- demo$data
  fit_plqr <- plqr_dnn(
    X = as.matrix(df_p[paste0("X", 1:2)]),
    D = as.matrix(df_p[paste0("D", 1:2)]),
    g = as.numeric(df_p$g_oracle),
    mD = as.matrix(df_p[paste0("mD", 1:2)]),
    Y = as.numeric(df_p$Y),
    theta_eval = demo$truth$theta_true, tau = 0.5,
    X0 = rbind(c(0, 0), c(0.5, 0.5)), k = 1L, s = 8L, se = TRUE
  )
  p <- plot_estimates(fit_plqr)
  expect_identical(class(p$facet)[1L], "FacetWrap")
  built <- ggplot2::ggplot_build(p)
  expect_identical(length(built$layout$layout$PANEL), 2L)

  # Single-coefficient fits stay unfaceted.
  p_np <- plot_estimates(.pbc_np_fit(m = 2L, se = FALSE))
  expect_identical(class(p_np$facet)[1L], "FacetNull")
})

test_that("significance rug agrees with the extractor's intervals", {
  skip_if_not_installed("ggplot2")
  set.seed(64)
  n <- 100L
  X <- matrix(rnorm(n * 2L), n, 2L,
              dimnames = list(NULL, c("X1", "X2")))
  psi <- rnorm(n, mean = 4 * X[, 1L], sd = 0.5)
  X0 <- cbind(seq(-1.5, 1.5, length.out = 8L), 0)
  fit <- cate_dnn(X = X, X0 = X0, pseudo_outcome = psi, k = 1L, s = 15L,
                  se = TRUE)

  df <- plot_data_estimates(fit)
  sig <- df[!is.na(df$ci_lower) &
              sign(df$ci_lower) == sign(df$ci_upper) &
              df$ci_lower != 0, , drop = FALSE]

  p <- plot_estimates(fit, shade_significant = TRUE)
  geoms <- vapply(p$layers, function(l) class(l$geom)[1L], character(1L))
  if (nrow(sig) > 0L) {
    expect_true("GeomRug" %in% geoms)
    rug_idx <- which(geoms == "GeomRug")[1L]
    built <- ggplot2::ggplot_build(p)
    expect_equal(sort(built$data[[rug_idx]]$x), sort(sig$x0_1))
  } else {
    expect_false("GeomRug" %in% geoms)
  }
})

test_that("built weight layers match the weights frame across views", {
  skip_if_not_installed("ggplot2")
  fit <- .pbc_np_fit(m = 2L, se = FALSE, return = c("W", "dx"))
  dw <- plot_data_weights(fit)

  # Kernel view: two rank-s radius guides, one per target, equal to the
  # rank-s distances in the frame.
  p_k <- plot_weights(fit, type = "kernel")
  geoms <- vapply(p_k$layers, function(l) class(l$geom)[1L], character(1L))
  vl <- p_k$layers[[which(geoms == "GeomVline")[1L]]]
  radii <- vapply(
    1:2, function(j) dw$distance[dw$target_id == j][fit$meta$s],
    numeric(1L)
  )
  expect_equal(sort(vl$data$xintercept), sort(radii))

  # Cumulative view: marker points sit on the cumulative curve at the
  # top-share ranks.
  p_c <- plot_weights(fit, type = "cumulative")
  geoms_c <- vapply(p_c$layers, function(l) class(l$geom)[1L], character(1L))
  pts <- p_c$layers[[which(geoms_c == "GeomPoint")[1L]]]$data
  expect_true(all(pts$rank %in% c(1L, 5L, 10L)))
  merged <- merge(pts, dw, by = c("target_id", "rank"),
                  suffixes = c("_layer", ""))
  expect_equal(merged$cum_share_layer, merged$cum_share)

  # Profile view built data reproduces the frame in rank order.
  p_p <- plot_weights(fit, type = "profile", targets = 1L)
  built <- ggplot2::ggplot_build(p_p)
  step <- built$data[[1L]]
  d1 <- dw[dw$target_id == 1L, ]
  expect_equal(step$x, d1$rank)
  expect_equal(step$y, d1$weight)
})

test_that("weight panels are class-agnostic and honest about provenance", {
  skip_if_not_installed("ggplot2")
  fit_cate <- .pbc_cate_fit(m = 2L, se = FALSE)
  p <- plot_weights(fit_cate, type = "profile")
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$caption, "s = 12 \\(manual\\)")

  # Default-s fits report paper_default provenance in the caption.
  set.seed(65)
  n <- 60L
  X <- matrix(rnorm(n * 2L), n, 2L)
  fit_def <- np_dnn(X = X, x0 = c(0, 0), Y = rnorm(n), k = 1L, se = FALSE)
  p_def <- plot_estimates(fit_def)
  expect_match(p_def$labels$caption, "paper_default")
})

test_that("default weight subsetting uses first, median-n_eff, and last", {
  skip_if_not_installed("ggplot2")
  fit_many <- .pbc_np_fit(m = 10L, se = FALSE)
  dw <- plot_data_weights(fit_many)
  all_targets <- sort(unique(dw$target_id))
  neff <- vapply(
    all_targets,
    function(j) 1 / sum(dw$weight[dw$target_id == j]^2),
    numeric(1L)
  )
  med <- all_targets[which.min(abs(neff - stats::median(neff)))]
  expected <- unique(c(all_targets[1L], med, all_targets[length(all_targets)]))

  p <- plot_weights(fit_many, type = "profile")
  expect_setequal(unique(p$data$target_id), expected)
  expect_true(all(c(1L, 10L) %in% unique(p$data$target_id)))
  expect_lte(length(unique(p$data$target_id)), 3L)
})

test_that("the ggplot2 guard fires before the missing-W error", {
  fit_now <- .pbc_np_fit(m = 1L, se = FALSE, return = c("ess"))
  testthat::local_mocked_bindings(
    .plot_pkg_installed = function(pkg) FALSE,
    .package = "RdistNN"
  )
  # Guard first: the error names ggplot2, not the weight matrix.
  expect_error(plot_weights(fit_now), "ggplot2")
})

test_that("autoplot passes arguments through to plot_estimates", {
  skip_if_not_installed("ggplot2")
  fit <- .pbc_np_fit(m = 6L, se = TRUE)
  p80 <- ggplot2::autoplot(fit, level = 0.8)
  expect_identical(attr(p80$data, "level"), 0.8)
  p_pr <- ggplot2::autoplot(fit, style = "pointrange")
  geoms <- vapply(p_pr$layers, function(l) class(l$geom)[1L], character(1L))
  expect_true("GeomPointrange" %in% geoms)

  # TDNN fits dispatch through the np_result autoplot method.
  set.seed(66)
  n <- 60L
  Xm <- matrix(rnorm(n * 2L), n, 2L)
  fit_tdnn <- np_tdnn(X = Xm, x0 = c(0, 0), Y = rnorm(n),
                      s1 = 4L, s2 = 8L, se = FALSE)
  p_tdnn <- ggplot2::autoplot(fit_tdnn)
  expect_s3_class(p_tdnn, "ggplot")
  expect_match(p_tdnn$labels$caption, "s1 = 4, s2 = 8")
})
