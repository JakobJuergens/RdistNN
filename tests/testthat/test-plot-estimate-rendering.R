.render_count_points <- function(grob) {
  if (inherits(grob, "points")) return(length(grob$x))
  children <- if (inherits(grob, "gtable")) {
    grob$grobs[grepl("^panel", grob$layout$name)]
  } else c(grob$grobs, as.list(grob$children))
  sum(vapply(children, .render_count_points, integer(1L)))
}

.render_ratio_fit <- function(grid = FALSE, profile = FALSE,
                              numerator = c(2, 2, 2, 2, 3, 6, 4, 9)) {
  args <- list(X = matrix(1:8, ncol = 1L),
    X0 = matrix(if (profile) c(0, 2, 4, 6, 8) else c(0, 8), ncol = 1L),
    numerator = numerator, denominator = rep(1, 8L),
    s = 7L, return = character())
  if (grid) {
    args$s <- NULL
    args$tuning_mode <- "grid"
    args$tuning <- new_tuning("rendering", eval = function(...) list(s = 7L),
      eval_grid = function(...) data.frame(s = c(6L, 7L)))
  }
  do.call(affine_ratio_dnn, args)
}

test_that("finite ratio estimates draw even when their Wald interval is unavailable", {
  .with_partial_result_warning({
  skip_if_not_installed("ggplot2")
  for (style in c("forest", "pointrange")) {
    p <- plot_estimates(.render_ratio_fit(), style = style)
    expect_true(any(is.na(p$data$ci_lower)))
    expect_true(any(is.finite(p$data$ci_lower)))
    expect_warning(g <- ggplot2::ggplotGrob(p), NA)
    expect_equal(.render_count_points(g), 2L)
    p <- plot_estimates(.render_ratio_fit(grid = TRUE), style = style, candidates = "all")
    expect_warning(g <- ggplot2::ggplotGrob(p), NA)
    expect_equal(.render_count_points(g), 4L)
  }
  })
})

test_that("plotting reports an actionable error for an old optional dependency", {
  .with_partial_result_warning({
  skip_if_not_installed("ggplot2")
  testthat::local_mocked_bindings(
    .plot_pkg_version = function(pkg) numeric_version("3.4.0"), .package = "RdistNN")
  expect_error(plot_estimates(.render_ratio_fit()), "ggplot2 >= 3.4.1.*3.4.0")
  })
})

test_that("profile estimates remain visible at targets without intervals", {
  .with_partial_result_warning({
  .with_distance_tie_warning({
  skip_if_not_installed("ggplot2")
  p <- plot_estimates(.render_ratio_fit(profile = TRUE), style = "profile")
  expect_warning(g <- ggplot2::ggplotGrob(p), NA)
  expect_equal(.render_count_points(g), 5L)
  })
  })
})

test_that("missing candidate intervals retain shared point positions", {
  .with_partial_result_warning({
  skip_if_not_installed("ggplot2")
  for (style in c("forest", "pointrange")) {
    p <- plot_estimates(.render_ratio_fit(grid = TRUE), style = style, candidates = "all")
    built <- ggplot2::ggplot_build(p)
    geoms <- vapply(p$layers, function(layer) class(layer$geom)[1L], character(1L))
    intervals <- built$data[[which(geoms == "GeomSegment")]]
    points <- built$data[[which(geoms == "GeomPoint")]]
    expect_equal(nrow(intervals), 2L)
    expect_equal(nrow(points), 4L)
    coordinate <- if (style == "forest") "y" else "x"
    expect_true(all(intervals[[coordinate]] %in% points[[coordinate]]))
  }
  })
})

test_that("profile ribbons do not bridge an unavailable interior interval", {
  .with_partial_result_warning({
  .with_distance_tie_warning({
  skip_if_not_installed("ggplot2")
  f <- .render_ratio_fit(profile = TRUE, numerator = c(2, 2, 2, 2, 4, 4, 4, 9))
  expect_identical(f$status$wald_available[c(3L, 4L, 5L)], c(TRUE, FALSE, TRUE))
  p <- plot_estimates(f, style = "profile")
  ribbon <- p$layers[[which(.plot_layer_geoms(p) == "GeomRibbon")]]$data
  for (group in unique(ribbon$interval_group)) {
    x <- ribbon$x0_1[ribbon$interval_group == group]
    expect_false(min(x) < 6 && max(x) > 6)
  }
  expect_warning(ggplot2::ggplotGrob(p), NA)
  })
  })
})
test_that("forest plots retain separate identities for repeated reporting targets", {
  skip_if_not_installed("ggplot2")
  X <- cbind(1:8, (1:8)^2)
  X0 <- rbind(c(.1, .2), c(3.1, 4.2), c(.1, .2))
  A <- rep(0:1, 4)
  tuning <- new_tuning("repeated_forest_targets", eval = function(...) list(s = 3L),
    eval_grid = function(...) data.frame(s = 2:4))
  constructors <- list(
    np = function(...) np_dnn(Y = 1:8, ...),
    cate = function(...) cate_dnn(pseudo_outcome = 1:8, ...),
    ratio = function(...) affine_ratio_dnn(numerator = 1:8, denominator = rep(1, 8), ...),
    clate = function(...) clate_dnn(Y = 2 * A + (1:8) / 10, D = A, A = A,
      e = rep(.5, 8), m_y0 = rep(0, 8), m_y1 = rep(2, 8),
      m_d0 = rep(0, 8), m_d1 = rep(1, 8), ...)
  )
  for (make in constructors) {
    fit <- make(X = X, X0 = X0, tuning = tuning, tuning_mode = "grid", se = FALSE)
    for (candidates in c("default", "all")) {
      p <- plot_estimates(fit, candidates = candidates)
      expect_identical(length(unique(p$data$target_label)), 3L)
      expect_match(as.character(p$data$target_label[p$data$target_id == 1L][1L]), "target 1")
      expect_match(as.character(p$data$target_label[p$data$target_id == 3L][1L]), "target 3")
      expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
      expect_no_error(ggplot2::autoplot(fit, candidates = candidates))
    }
  }
})
