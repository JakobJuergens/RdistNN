.render_count_points <- function(grob) {
  if (inherits(grob, "points")) return(length(grob$x))
  children <- if (inherits(grob, "gtable")) {
    grob$grobs[grepl("^panel", grob$layout$name)]
  } else c(grob$grobs, as.list(grob$children))
  sum(vapply(children, .render_count_points, integer(1L)))
}

.render_ratio_fit <- function(grid = FALSE, profile = FALSE) {
  args <- list(X = matrix(1:8, ncol = 1L),
    X0 = matrix(if (profile) c(0, 2, 4, 6, 8) else c(0, 8), ncol = 1L),
    numerator = c(2, 2, 2, 2, 3, 6, 4, 9), denominator = rep(1, 8L),
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

test_that("plotting reports an actionable error for an old optional dependency", {
  skip_if_not_installed("ggplot2")
  testthat::local_mocked_bindings(
    .plot_pkg_version = function(pkg) numeric_version("3.4.0"), .package = "RdistNN")
  expect_error(plot_estimates(.render_ratio_fit()), "ggplot2 >= 3.4.1.*3.4.0")
})

test_that("profile estimates remain visible at targets without intervals", {
  skip_if_not_installed("ggplot2")
  p <- plot_estimates(.render_ratio_fit(profile = TRUE), style = "profile")
  expect_warning(g <- ggplot2::ggplotGrob(p), NA)
  expect_equal(.render_count_points(g), 5L)
})

test_that("missing candidate intervals retain shared point positions", {
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

test_that("profile ribbons do not bridge an unavailable interior interval", {
  skip_if_not_installed("ggplot2")
  f <- .render_ratio_fit(profile = TRUE)
  f$status$wald_available[4L] <- FALSE
  p <- plot_estimates(f, style = "profile")
  ribbon <- p$layers[[which(.plot_layer_geoms(p) == "GeomRibbon")]]$data
  for (group in unique(ribbon$interval_group)) {
    x <- ribbon$x0_1[ribbon$interval_group == group]
    expect_false(min(x) < 6 && max(x) > 6)
  }
  expect_warning(ggplot2::ggplotGrob(p), NA)
})
