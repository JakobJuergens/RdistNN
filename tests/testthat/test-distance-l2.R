test_that("L2 matches basic expectations", {
  X <- gen_X()
  x0 <- gen_x0(p = ncol(X))
  d <- dist_l2(X, x0)
  ds <- dist_l2(X, x0, squared = TRUE)
  expect_true(all(d >= 0))
  expect_true(all.equal(d^2, ds, tolerance = 1e-12) == TRUE)
})

test_that("L2 order is invariant to squared flag", {
  X <- gen_X()
  x0 <- gen_x0(p = ncol(X))
  o1 <- .distance_ranking(X, x0, distance_l2(squared = FALSE))$ord
  o2 <- .distance_ranking(X, x0, distance_l2(squared = TRUE))$ord
  expect_identical(o1, o2)
})

# Battery: unsquared
run_distance_battery(
  metric_name = "l2",
  mk_spec = function(...) distance_l2(squared = FALSE),
  args = list(),
  expect_sq = FALSE
)

# Battery: squared
run_distance_battery(
  metric_name = "l2_sq",
  mk_spec = function(...) distance_l2(squared = TRUE),
  args = list(),
  expect_sq = TRUE
)
