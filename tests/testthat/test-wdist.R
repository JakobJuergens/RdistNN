# tests/testthat/test-wdist.R
# -------------------------------------------------------------------
# Distance-weight scaffold (R/wdist.R)
# Tags:
#   [smoke]   Construction & printing
#   [single]  eval_dx path (one target)
#   [batch]   eval_dx_batch path (many targets)
#   [setup]   Context setup usage
#   [errors]  Input validation & API guards
#   [conv]    Convenience helpers & predicates
# -------------------------------------------------------------------

test_that("[smoke] can construct a minimal R dist_weight and print it", {
  softmax_w <- function(dx, ctx = NULL) {
    z <- exp(-as.numeric(dx))
    z / sum(z)
  }
  dw <- new_dist_weight(name = "softmax", eval_dx = softmax_w)

  expect_s3_class(dw, "dist_weight")
  expect_identical(dw$name, "softmax")
  expect_output(print(dw), "<dist_weight: softmax>")
  expect_output(print(dw), "R eval_dx")
})

test_that("[single] eval_dx returns numeric(n) nonnegative with unit sum", {
  dw <- new_dist_weight(
    name    = "softmax",
    eval_dx = function(dx, ctx = NULL) { z <- exp(-dx); z / sum(z) }
  )
  dx <- c(1, 3, 2, 0)
  w  <- dist_weight_eval(dw, dx)
  expect_type(w, "double")
  expect_length(w, length(dx))
  expect_true(all(w >= 0))
  expect_equal(sum(w), 1, tolerance = 1e-12)
})

test_that("[batch] eval_dx_batch returns numeric(n x m) with column sums 1", {
  # Batch-only implementation
  dwb <- new_dist_weight(
    name = "softmax_batch",
    eval_dx_batch = function(DX, ctx = NULL) {
      Z <- exp(-DX)
      sweep(Z, 2L, colSums(Z), "/")
    }
  )
  DX <- matrix(c(0, 1, 2, 3,
                 3, 2, 1, 0), nrow = 4, ncol = 2)
  W <- dist_weight_eval_batch(dwb, DX)
  expect_true(is.matrix(W))
  expect_equal(dim(W), dim(DX))
  expect_true(all(W >= 0))
  expect_equal(colSums(W), rep(1, ncol(DX)), tolerance = 1e-12)

  # Calling single-target path on a batch-only dw should error
  expect_error(dist_weight_eval(dwb, DX[,1]), "No R `eval_dx` provided")
})

test_that("[setup] dist_weight_setup provides ctx used inside eval", {
  dw <- new_dist_weight(
    name = "scaled_gauss",
    setup = function(data, ...) list(h = 2.0),
    eval_dx = function(dx, ctx = NULL) {
      h <- ctx$h
      w <- exp(-dx / (2*h*h))
      w / sum(w)
    }
  )
  ctx <- dist_weight_setup(dw, matrix(0, 5, 2))
  expect_equal(ctx$h, 2.0, tolerance = 0)

  dx <- c(0, 1, 4)
  w  <- dist_weight_eval(dw, dx, ctx = ctx)
  # Compare to explicit formula with h=2
  w_star <- {
    z <- exp(-dx / 8)
    z / sum(z)
  }
  expect_equal(w, w_star, tolerance = 1e-15)
})

test_that("[errors] input validation for dx/DX and API guards", {
  dw <- new_dist_weight(
    name    = "softmax",
    eval_dx = function(dx, ctx = NULL) { z <- exp(-dx); z / sum(z) }
  )
  expect_error(dist_weight_eval(dw, numeric(0)), "positive length")
  expect_error(dist_weight_eval(dw, c(1, NA)), "cannot contain NA")
  expect_error(dist_weight_eval(dw, c(1, Inf)), "must be finite")

  # NEW: single-target path rejects matrices
  expect_error(
    dist_weight_eval(dw, matrix(1, 2, 2)),
    "dx must be a numeric vector; use dist_weight_eval_batch\\(\\) for matrices"
  )

  dwb <- new_dist_weight(
    name = "softmax_batch",
    eval_dx_batch = function(DX, ctx = NULL) {
      Z <- exp(-DX); sweep(Z, 2L, colSums(Z), "/")
    }
  )
  expect_error(dist_weight_eval_batch(dwb, matrix(NA_real_, 2, 2)), "must be finite")
  expect_error(dist_weight_eval_batch(dwb, matrix(1, 0, 2)), "positive nrow and ncol")
  expect_error(dist_weight_eval(dw, matrix(1, 2, 2)), "dx must be a numeric vector")  # duplicate guard ok

  # Returned shape/content must be checked
  bad <- new_dist_weight(
    name = "bad_single",
    eval_dx = function(dx, ctx = NULL) c(1, 2) # wrong length
  )
  expect_error(dist_weight_eval(bad, c(1, 2, 3)), "must return numeric\\(length\\(dx\\)\\)")

  badb <- new_dist_weight(
    name = "bad_batch",
    eval_dx_batch = function(DX, ctx = NULL) matrix(1, 3, 1) # wrong dims
  )
  expect_error(dist_weight_eval_batch(badb, matrix(1, 2, 2)), "must return a finite numeric matrix")
})


test_that("[conv] as_dist_weight wraps a plain function; predicate works; randomized flag prints", {
  dw <- as_dist_weight(function(dx, ctx = NULL) rep(1/length(dx), length(dx)), name = "uniform")
  expect_s3_class(dw, "dist_weight")
  expect_true(is_dist_weight(dw))
  expect_output(print(dw), "<dist_weight: uniform>")
  expect_output(print(dw), "R eval_dx")

  # randomized printing
  dwr <- new_dist_weight(
    name = "rand",
    eval_dx = function(dx, ctx = NULL) rep(1/length(dx), length(dx)),
    randomized = TRUE
  )
  expect_output(print(dwr), "(randomized)")
})
