# tests/testthat/test-dist-l2-scaled.R
# -------------------------------------------------------------------
# Scaled L2: shapes, parity, batch, registry-agnostic use
# Tags:
#   [engine]  shape/finite & squared toggle
#   [compiled] parity vs R fallback
#   [batch]   batch equals columnwise calls
#   [custom]  custom scale_vec; parity to unscaled when scale==1
# -------------------------------------------------------------------

test_that("[engine] distance_l2_scaled works; squared toggle", {
  set.seed(1)
  X <- matrix(rnorm(30), 10, 3)
  x0 <- rnorm(3)
  dspec <- distance_l2_scaled("sd", squared = TRUE)
  ctx <- distance_setup(dspec, X)

  d_sq <- distance_eval(dspec, X, x0, ctx)
  d <- {
    dspec2 <- distance_l2_scaled("sd", squared = FALSE)
    ctx2 <- distance_setup(dspec2, X)
    distance_eval(dspec2, X, x0, ctx2)
  }
  expect_length(d_sq, nrow(X))
  expect_true(all(is.finite(d_sq)))
  expect_true(all(is.finite(d)))
  expect_equal(d_sq, d^2, tolerance = 1e-12)
})

test_that("[compiled] rustat_dist_l2_scaled parity with R fallback", {
  set.seed(42)
  X <- matrix(rnorm(40), 20, 2)
  x0 <- rnorm(2)
  dspec <- distance_l2_scaled("mad", squared = TRUE)
  ctx <- distance_setup(dspec, X)

  # Reference R fallback
  ref_sq <- {
    dif <- sweep(X, 2L, x0, "-")
    dif <- sweep(dif, 2L, ctx$inv_s, "*")
    rowSums(dif * dif)
  }
  skip_if_not(exists("rustat_dist_l2_scaled", mode = "function"))
  out_sq <- distance_eval(dspec, X, x0, ctx)
  expect_equal(out_sq, ref_sq, tolerance = 1e-12)
})

test_that("[batch] batch equals columnwise calls", {
  set.seed(7)
  X <- matrix(rnorm(24), 8, 3)
  X0 <- matrix(rnorm(6), 2, 3)
  dspec <- distance_l2_scaled("iqr", squared = TRUE)
  ctx <- distance_setup(dspec, X)

  M <- distance_eval_batch(dspec, X, X0, ctx)
  M2 <- sapply(seq_len(nrow(X0)), function(j) distance_eval(dspec, X, X0[j, ], ctx))
  expect_equal(M, M2, tolerance = 1e-12)
})

test_that("[custom] custom scale_vec works; scale==1 matches unscaled L2", {
  set.seed(3)
  X <- matrix(rnorm(20), 10, 2)
  x0 <- rnorm(2)

  # custom all-ones -> same as unscaled squared L2
  dspec_un <- distance_l2(squared = TRUE)
  dspec_one <- distance_l2_scaled("custom", scale_vec = c(1, 1), squared = TRUE)
  ctx_one <- distance_setup(dspec_one, X)

  d_un <- distance_eval(dspec_un, X, x0)
  d_one <- distance_eval(dspec_one, X, x0, ctx_one)
  expect_equal(d_one, d_un, tolerance = 1e-12)

  # mismatch length should error
  expect_error(
    distance_l2_scaled("custom", scale_vec = 1, squared = TRUE) |> distance_setup(X),
    "length p = ncol\\(data\\)"
  )
})
