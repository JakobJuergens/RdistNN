test_that("Scaled L2 equals unscaled when scale=1 (both flags)", {
  X <- gen_X()
  x0 <- gen_x0(p = ncol(X))
  inv1 <- rep(1, ncol(X))
  expect_equal(dist_l2_scaled(X, x0, inv1, FALSE), dist_l2(X, x0, FALSE))
  expect_equal(dist_l2_scaled(X, x0, inv1, TRUE), dist_l2(X, x0, TRUE))
})

test_that("Scaled L2 handles eps-floor (finite output)", {
  X0 <- gen_X()
  mod <- gen_constcol_and_inv(X0, eps = 1e-12)
  X <- mod$X
  inv <- mod$inv
  x0 <- gen_x0(p = ncol(X))
  d <- dist_l2_scaled(X, x0, inv, squared = FALSE)
  ds <- dist_l2_scaled(X, x0, inv, squared = TRUE)
  expect_true(all(is.finite(d)) && all(is.finite(ds)))
  expect_true(all(d >= 0) && all(ds >= 0))
  expect_true(all.equal(d^2, ds, tolerance = 1e-12) == TRUE)
})

# Battery: unsquared (data-driven scaling)
run_distance_battery(
  metric_name = "l2_scaled",
  mk_spec = function(...) distance_l2_scaled(rule = "sd", squared = FALSE),
  orderer = function(X, x0, ...) order_by_l2_scaled(X, x0, rule = "sd", squared = FALSE),
  args = list(),
  expect_sq = FALSE
)

# Battery: squared
run_distance_battery(
  metric_name = "l2_scaled_sq",
  mk_spec = function(...) distance_l2_scaled(rule = "sd", squared = TRUE),
  orderer = function(X, x0, ...) order_by_l2_scaled(X, x0, rule = "sd", squared = TRUE),
  args = list(),
  expect_sq = TRUE
)

# Optional: dedicated C++ vs wrapper (when available)
test_that("Scaled L2: C++ fast path equals wrapper path (if available)", {
  skip_if_not(exists("rustat_dist_l2_scaled", mode = "function"))
  X <- gen_X(seed = 21)
  x0 <- gen_x0(p = ncol(X), seed = 22)
  inv <- rep(3, ncol(X))
  d_cpp <- rustat_dist_l2_scaled(as.matrix(X), as.numeric(x0), as.numeric(inv), FALSE)
  d_wr <- dist_l2_scaled(X, x0, inv, FALSE)
  ds_cpp <- rustat_dist_l2_scaled(as.matrix(X), as.numeric(x0), as.numeric(inv), TRUE)
  ds_wr <- dist_l2_scaled(X, x0, inv, TRUE)
  expect_equal(d_cpp, d_wr)
  expect_equal(ds_cpp, ds_wr)
})
