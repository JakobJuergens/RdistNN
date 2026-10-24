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

test_that("Scaled L1 rejects non-finite and non-numeric scale inputs", {
  X <- gen_X()
  x0 <- gen_x0(p = ncol(X))
  inv <- rep(1, ncol(X))

  inv_bad <- inv
  inv_bad[1L] <- Inf
  expect_error(dist_l1_scaled(X, x0, inv_bad), "inv_s.*non-finite", ignore.case = TRUE)
  expect_error(dist_l1_scaled(X, x0, as.character(inv)), "inv_s.*numeric", ignore.case = TRUE)

  expect_error(
    order_by_l1_scaled(X, x0, rule = "custom", scale_vec = replace(inv, 1L, Inf)),
    "scale_vec.*non-finite",
    ignore.case = TRUE
  )
  expect_error(
    order_by_l1_scaled(X, x0, rule = "custom", scale_vec = as.character(inv)),
    "scale_vec.*numeric",
    ignore.case = TRUE
  )
  expect_error(
    order_by_l1_scaled(X, x0, rule = "custom", scale_vec = inv, eps = 0),
    "eps.*positive finite",
    ignore.case = TRUE
  )

  spec <- distance_l1_scaled(rule = "sd")
  bad_ctx <- list(inv_s = replace(inv, 1L, NaN))
  expect_error(spec$eval(X, x0, bad_ctx), "x0/inv_s.*non-finite", ignore.case = TRUE)
  expect_error(spec$eval_batch(X, matrix(x0, nrow = 1L), bad_ctx), "inv_s.*non-finite", ignore.case = TRUE)
})

# Battery: use data-driven scaling rule to avoid dim mismatch
run_distance_battery(
  metric_name = "l1_scaled",
  mk_spec = function(...) distance_l1_scaled(rule = "sd"),
  orderer = function(X, x0, ...) order_by_l1_scaled(X, x0, rule = "sd"),
  args = list(),
  expect_sq = NULL
)

test_that("Scaled L1: native and forced R paths agree", {
  withr::local_options(RdistNN.use_compiled = FALSE)
  X <- gen_X(seed = 11)
  x0 <- gen_x0(p = ncol(X), seed = 12)
  inv <- rep(2, ncol(X))
  d_cpp <- rustat_dist_l1_scaled(as.matrix(X), as.numeric(x0), as.numeric(inv))
  d_wr <- dist_l1_scaled(X, x0, inv)
  expect_equal(d_cpp, d_wr)
})
