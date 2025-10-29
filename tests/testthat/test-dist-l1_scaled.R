# tests/testthat/test-dist-l1-scaled.R
# -------------------------------------------------------------------
# Scaled L1: shapes, compiled parity, batch, and custom scales
# Tags:
#   [engine]   shape/finite
#   [compiled] parity vs R fallback
#   [batch]    batch equals columnwise calls
#   [custom]   custom scale_vec; parity to unscaled when scale==1
# -------------------------------------------------------------------

test_that("[engine] distance_l1_scaled works and returns finite distances", {
  set.seed(1)
  X <- matrix(rnorm(30), 10, 3)
  x0 <- rnorm(3)
  dspec <- distance_l1_scaled("sd")
  ctx <- distance_setup(dspec, X)

  d <- distance_eval(dspec, X, x0, ctx)
  expect_length(d, nrow(X))
  expect_true(all(is.finite(d)))
})

test_that("[compiled] rustat_dist_l1_scaled parity with R fallback", {
  set.seed(42)
  X <- matrix(rnorm(40), 20, 2)
  x0 <- rnorm(2)
  dspec <- distance_l1_scaled("mad")
  ctx <- distance_setup(dspec, X)

  # R fallback reference
  ref <- {
    dif <- sweep(X, 2L, x0, "-")
    dif <- sweep(dif, 2L, ctx$inv_s, "*")
    rowSums(abs(dif))
  }

  skip_if_not(exists("rustat_dist_l1_scaled", mode = "function"))
  out <- distance_eval(dspec, X, x0, ctx)
  expect_equal(out, ref, tolerance = 0)
})

test_that("[batch] batch equals columnwise calls", {
  set.seed(7)
  X <- matrix(rnorm(24), 8, 3)
  X0 <- matrix(rnorm(6), 2, 3)
  dspec <- distance_l1_scaled("iqr")
  ctx <- distance_setup(dspec, X)

  M <- distance_eval_batch(dspec, X, X0, ctx)
  M2 <- sapply(seq_len(nrow(X0)), function(j) distance_eval(dspec, X, X0[j, ], ctx))
  expect_equal(M, M2, tolerance = 0)
})

test_that("[custom] custom scale_vec works; all-ones matches unscaled L1", {
  set.seed(3)
  X <- matrix(rnorm(20), 10, 2)
  x0 <- rnorm(2)

  dspec_un <- distance_l1() # unscaled L1
  dspec_one <- distance_l1_scaled("custom", scale_vec = c(1, 1))
  ctx_one <- distance_setup(dspec_one, X)

  d_un <- distance_eval(dspec_un, X, x0)
  d_one <- distance_eval(dspec_one, X, x0, ctx_one)

  expect_equal(d_one, d_un, tolerance = 0)

  # length mismatch should error on setup
  expect_error(
    distance_l1_scaled("custom", scale_vec = 1) |> distance_setup(X),
    "length p = ncol\\(data\\)"
  )
})
