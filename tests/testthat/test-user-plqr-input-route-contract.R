.plqr_input_route_fixture <- function() {
  n <- 8L
  data <- data.frame(x = seq_len(n), d = seq_len(n), g = 0, md = 0,
    y = seq_len(n), density = .5)
  list(
    data = list(data = data, x_cols = "x", D_cols = "d", g_col = "g",
      mD_cols = "md", y_col = "y", theta_eval = .5, x0 = 0, s = 3L,
      se = FALSE),
    matrix = list(X = matrix(data$x), D = matrix(data$d), g = data$g,
      mD = matrix(data$md), Y = data$y, theta_eval = .5, x0 = 0, s = 3L,
      se = FALSE)
  )
}

test_that("PLQR rejects ignored matrix arguments with data inputs", {
  args <- .plqr_input_route_fixture()$data
  for (arg in c("X", "D", "g", "mD", "Y")) {
    conflicting <- args
    conflicting[[arg]] <- "invalid direct input"
    expect_error(do.call(plqr_dnn, conflicting),
      paste0("`", arg, "`.*cannot be combined with `data`"))
  }
})

test_that("PLQR requires data for its column selectors", {
  args <- .plqr_input_route_fixture()$matrix
  for (arg in c("x_cols", "D_cols", "g_col", "mD_cols", "y_col")) {
    conflicting <- args
    conflicting[[arg]] <- "missing_column"
    expect_error(do.call(plqr_dnn, conflicting), "data.*column selectors")
  }
})

test_that("PLQR route validation preserves direct densities with data inputs", {
  fixture <- .plqr_input_route_fixture()
  tabular <- do.call(plqr_dnn, fixture$data)
  numeric <- do.call(plqr_dnn, fixture$matrix)
  expect_equal(unname(tabular$moment), unname(numeric$moment), tolerance = 1e-12)

  args <- fixture$data
  args$se <- TRUE
  direct <- do.call(plqr_dnn, c(args, list(density = rep(.5, 8L))))
  column <- do.call(plqr_dnn, c(args, list(density_col = "density")))
  expect_equal(direct$moment, column$moment, tolerance = 1e-12)
  expect_equal(direct$theta_vcov, column$theta_vcov, tolerance = 1e-12)
  expect_true(all(is.finite(direct$theta_se)))
})
