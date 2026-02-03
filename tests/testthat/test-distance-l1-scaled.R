test_that("Scaled L1 equals unscaled when scale=1", {
  X <- gen_X()
  x0 <- gen_x0(p = ncol(X))
  inv1 <- rep(1, ncol(X))
  expect_equal(dist_l1_scaled(X, x0, inv1), dist_l1(X, x0))
})

test_that("Scaled L1 handles eps-floor (finite output)", {
  X0 <- gen_X()
  mod <- gen_constcol_and_inv(X0, eps = 1e-12)
  X <- mod$X
  inv <- mod$inv
  x0 <- gen_x0(p = ncol(X))
  d <- dist_l1_scaled(X, x0, inv)
  expect_true(all(is.finite(d)))
  expect_true(all(d >= 0))
})

# Battery: use data-driven scaling rule to avoid dim mismatch
run_distance_battery(
  metric_name = "l1_scaled",
  mk_spec = function(...) distance_l1_scaled(rule = "sd"),
  orderer = function(X, x0, ...) order_by_l1_scaled(X, x0, rule = "sd"),
  args = list(),
  expect_sq = NULL
)

# If the dedicated C++ symbol exists, check C++ vs wrapper
test_that("Scaled L1: C++ fast path equals wrapper path (if available)", {
  skip_if_not(exists("rustat_dist_l1_scaled", mode = "function"))
  X <- gen_X(seed = 11)
  x0 <- gen_x0(p = ncol(X), seed = 12)
  inv <- rep(2, ncol(X))
  d_cpp <- rustat_dist_l1_scaled(as.matrix(X), as.numeric(x0), as.numeric(inv))
  d_wr <- dist_l1_scaled(X, x0, inv)
  expect_equal(d_cpp, d_wr)
})
