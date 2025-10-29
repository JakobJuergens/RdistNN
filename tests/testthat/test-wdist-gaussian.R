# tests/testthat/test-wdist-gaussian-l2.R
# -------------------------------------------------------------------
# Gaussian kernel distance-weight over L2 distances
# Tags:
#   [smoke]   Construction & printing
#   [single]  Monotonicity and normalization (squared)
#   [batch]   Column sums, and normalize=FALSE parity
#   [equiv]   Equivalence: squared=TRUE vs squared=FALSE with sqrt inputs
#   [underflow]  Extreme small h -> 1-NN fallback
#   [api]     Integration with weights_matrix() using distance_l2
# -------------------------------------------------------------------

test_that("[smoke] construct and print gaussian_l2 dist_weight", {
  dw <- wdist_gaussian_l2(h = 0.5, squared = TRUE)
  expect_s3_class(dw, "dist_weight")
  expect_output(print(dw), "<dist_weight: gaussian_l2")
  expect_output(print(dw), "R eval_dx_batch")
})

test_that("[single] weights are monotone decreasing in (squared) distance and sum to 1", {
  dw <- wdist_gaussian_l2(h = 1, squared = TRUE, normalize = TRUE)
  dx2 <- c(0, 1, 4, 9) # squared L2 distances
  w <- dist_weight_eval(dw, dx2)

  expect_length(w, length(dx2))
  expect_true(all(w >= 0))
  expect_equal(sum(w), 1, tolerance = 1e-12)

  # Monotone: closer -> larger weight
  expect_true(w[1] > w[2] && w[2] > w[3] && w[3] > w[4])
})

test_that("[batch] column sums are 1; normalize=FALSE matches raw exp(-DX/(2h^2))", {
  h <- 0.7
  dw <- wdist_gaussian_l2(h = h, squared = TRUE, normalize = TRUE)

  DX <- matrix(c(
    0, 1, 4, 9,
    3, 2, 1, 0
  ), nrow = 4, ncol = 2, byrow = FALSE)
  W <- dist_weight_eval_batch(dw, DX)
  expect_equal(dim(W), dim(DX))
  expect_true(all(W >= 0))
  expect_equal(colSums(W), rep(1, ncol(DX)), tolerance = 1e-12)

  # normalize = FALSE parity check
  dwn <- wdist_gaussian_l2(h = h, squared = TRUE, normalize = FALSE)
  Wn <- dist_weight_eval_batch(dwn, DX)
  Wstar <- exp(-DX / (2 * h * h))
  expect_equal(Wn, Wstar, tolerance = 1e-15)
})

test_that("[equiv] squared=TRUE vs squared=FALSE with matched inputs", {
  h <- 0.9
  dx <- c(0.0, 0.5, 1.25, 2.0) # plain L2 distances

  dw_sq <- wdist_gaussian_l2(h = h, squared = TRUE, normalize = TRUE)
  dw_ln <- wdist_gaussian_l2(h = h, squared = FALSE, normalize = TRUE)

  w_sq <- dist_weight_eval(dw_sq, dx^2) # pass squared distances
  w_ln <- dist_weight_eval(dw_ln, dx) # pass plain distances

  expect_equal(w_sq, w_ln, tolerance = 1e-15)

  # batch version
  DX <- cbind(dx, rev(dx))
  Wsq <- dist_weight_eval_batch(dw_sq, DX^2)
  Wln <- dist_weight_eval_batch(dw_ln, DX)
  expect_equal(Wsq, Wln, tolerance = 1e-15)
})

test_that("[underflow] extreme small h triggers 1-NN fallback (no zero distances)", {
  # Choose all distances >= 1 and tiny bandwidth so exp terms underflow to 0
  dx2 <- c(1, 4, 9, 16) # squared distances; min at index 1
  h <- 1e-3 # exp(-dx2/(2h^2)) ~ exp(-O(1e6)) -> 0
  dw <- wdist_gaussian_l2(h = h, squared = TRUE, normalize = TRUE)

  w <- dist_weight_eval(dw, dx2)
  expect_equal(w, c(1, 0, 0, 0))

  # Batch: each column’s min should get all mass
  DX <- cbind(
    c(4, 9, 16, 25), # min at row 1
    c(9, 1, 16, 25), # min at row 2
    c(25, 9, 4, 36) # min at row 3
  )
  W <- dist_weight_eval_batch(dw, DX)
  expect_equal(W[, 1], c(1, 0, 0, 0))
  expect_equal(W[, 2], c(0, 1, 0, 0))
  expect_equal(W[, 3], c(0, 0, 1, 0))
})

test_that("[api] integrates with weights_matrix using distance_l2(squared=TRUE)", {
  set.seed(1)
  X <- matrix(rnorm(40), 20, 2)
  X0 <- matrix(rnorm(6), 3, 2)

  dx_spec <- distance_l2(squared = TRUE) # produces DX
  dw <- wdist_gaussian_l2(h = 0.5, squared = TRUE, normalize = TRUE)

  out <- weights_matrix(
    X, X0,
    ranking = dx_spec, # compute distances via distance_l2
    weight = dw, # distance-based weighting rule
    return = c("W", "dx")
  )
  expect_s3_class(out, "weights_result")
  expect_equal(dim(out$W), c(nrow(X), nrow(X0)))
  expect_true(all(out$W >= 0))
  expect_equal(colSums(out$W), rep(1, nrow(X0)), tolerance = 1e-12)
  expect_true(is.matrix(out$dx))
})
