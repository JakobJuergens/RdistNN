# tests/testthat/test-rank-api.R
# ------------------------------------------------------------------------------

# Small helper to build a random ranks matrix (each col a permutation 1..n)
.rand_ranks_mat <- function(n, m) {
  R <- matrix(NA_integer_, n, m)
  set.seed(123)
  for (j in seq_len(m)) R[, j] <- sample.int(n, n, replace = FALSE)
  R
}

# Small helper to build a random ranks matrix (each col a permutation 1..n)
.rand_ranks_mat <- function(n, m) {
  R <- matrix(NA_integer_, n, m)
  set.seed(123)
  for (j in seq_len(m)) R[, j] <- sample.int(n, n, replace = FALSE)
  R
}

test_that("is_rank_weight() and as_rank_weight() basics", {
  skip_if_not(exists("is_rank_weight"), "is_rank_weight() not found")
  skip_if_not(exists("as_rank_weight"), "as_rank_weight() not found")

  # Simple harmonic rank weight w_i ∝ 1/r_i
  rw <- as_rank_weight(function(r, ctx = NULL) {
    w <- 1 / as.numeric(r)
    w / sum(w)
  }, name = "inv_rank")

  expect_true(is_rank_weight(rw))
  expect_equal(rw$name, "inv_rank")

  # Single-column evaluation
  n <- 7L
  r <- sample.int(n)
  w <- rw$eval(r)
  expect_type(w, "double")
  expect_length(w, n)
  expect_true(all(is.finite(w)))
  expect_equal(sum(w), 1, tolerance = 1e-12)
})

test_that("rank_weight_eval() and rank_weight_eval_batch() work (eval only path)", {
  skip_if_not(exists("as_rank_weight"), "as_rank_weight() not found")
  skip_if_not(exists("rank_weight_eval"), "rank_weight_eval() not found")
  skip_if_not(exists("rank_weight_eval_batch"), "rank_weight_eval_batch() not found")

  # eval-only spec (no eval_batch): batch should fall back to apply(eval, 2, ...)
  rw_eval_only <- as_rank_weight(function(r, ctx = NULL) {
    w <- 1 / as.numeric(r)
    w / sum(w)
  }, name = "eval_only")

  n <- 8L; m <- 3L
  R <- .rand_ranks_mat(n, m)

  Wb <- rank_weight_eval_batch(rw_eval_only, R)
  expect_equal(dim(Wb), c(n, m))
  expect_true(all(is.finite(Wb)))
  expect_equal(unname(colSums(Wb)), rep(1, m), tolerance = 1e-12)

  # Single-target matches the batch first column
  w1 <- rank_weight_eval(rw_eval_only, R[, 1L])
  expect_equal(w1, Wb[, 1L])
})

test_that("rank_weight_eval() falls back to eval_batch when eval is missing", {
  skip_if_not(exists("rank_weight_eval"), "rank_weight_eval() not found")
  skip_if_not(exists("rank_weight_eval_batch"), "rank_weight_eval_batch() not found")

  # batch-only spec (no eval): single-target should route through eval_batch
  rw_batch_only <- list(
    name = "batch_only",
    eval_batch = function(R, ctx = NULL) {
      # w_i ∝ 1 / r_i^2
      W <- 1 / (R^2)
      W / rep(colSums(W), each = nrow(W))
    }
  )
  rw_batch_only$.__kind <- "rank_weight"
  class(rw_batch_only) <- c("rank_weight", class(rw_batch_only))

  n <- 9L; m <- 4L
  R <- .rand_ranks_mat(n, m)

  # Batch
  Wb <- rank_weight_eval_batch(rw_batch_only, R)
  expect_equal(dim(Wb), c(n, m))
  expect_equal(unname(colSums(Wb)), rep(1, m), tolerance = 1e-12)

  # Single target path
  w1 <- rank_weight_eval(rw_batch_only, R[, 2L])
  expect_equal(w1, Wb[, 2L])
})

test_that("ctx is forwarded to eval / eval_batch", {
  skip_if_not(exists("as_rank_weight"), "as_rank_weight() not found")

  # ctx: power alpha; w_i ∝ 1 / r_i^alpha
  rw_ctx <- as_rank_weight(function(r, ctx = NULL) {
    alpha <- if (is.null(ctx) || is.null(ctx$alpha)) 1 else ctx$alpha
    w <- 1 / (as.numeric(r)^alpha)
    w / sum(w)
  }, name = "ctx_power")

  n <- 6L; m <- 2L
  R <- .rand_ranks_mat(n, m)

  W1 <- rank_weight_eval_batch(rw_ctx, R, ctx = list(alpha = 1))
  W2 <- rank_weight_eval_batch(rw_ctx, R, ctx = list(alpha = 2))

  expect_equal(dim(W1), c(n, m))
  expect_equal(dim(W2), c(n, m))
  # Different alpha should change weights
  expect_false(isTRUE(all.equal(W1, W2)))
  # Still normalized
  expect_equal(unname(colSums(W1)), rep(1, m), tolerance = 1e-12)
  expect_equal(unname(colSums(W2)), rep(1, m), tolerance = 1e-12)
})

test_that("register_rank(), get_rank(), rank_names()/list_ranks() work (new + legacy)", {
  skip_if_not(exists("register_rank"), "register_rank() not found")
  skip_if_not(exists("get_rank"), "get_rank() not found")
  skip_if_not(exists("rank_names"), "rank_names() not found")
  skip_if_not(exists("list_ranks"), "list_ranks() not found")

  # New form: (name, spec function)
  nm_new <- paste0("rw_new_", as.integer(runif(1, 1e6, 2e6)))
  spec_new <- function(r, ctx = NULL) {
    w <- 1 / as.numeric(r)
    w / sum(w)
  }
  register_rank(nm_new, spec_new)
  got_new <- get_rank(nm_new)
  expect_true(is_rank_weight(got_new))
  expect_equal(got_new$name, nm_new)
  expect_true(nm_new %in% rank_names())
  expect_true(nm_new %in% list_ranks())

  # Legacy form: register_rank(spec) where spec has $name
  nm_legacy <- paste0("rw_legacy_", as.integer(runif(1, 2e6, 3e6)))
  spec_legacy <- as_rank_weight(function(r, ctx = NULL) {
    w <- 1 / (as.numeric(r)^2)
    w / sum(w)
  }, name = nm_legacy)
  register_rank(spec_legacy)  # legacy one-arg
  got_legacy <- get_rank(nm_legacy)
  expect_true(is_rank_weight(got_legacy))
  expect_equal(got_legacy$name, nm_legacy)

  # Unknown name errors list known entries
  expect_error(get_rank("__no_such_rank__"), "Unknown rank_weight")
})

test_that("input validation: bad ranks and bad specs produce clear errors", {
  skip_if_not(exists("rank_weight_eval"), "rank_weight_eval() not found")
  skip_if_not(exists("rank_weight_eval_batch"), "rank_weight_eval_batch() not found")
  skip_if_not(exists("as_rank_weight"), "as_rank_weight() not found")

  rw <- as_rank_weight(function(r, ctx = NULL) {
    w <- 1 / as.numeric(r)
    w / sum(w)
  }, name = "inv_rank")

  # Non-permutation
  expect_error(rank_weight_eval(rw, c(1L, 1L, 2L)), "permutation", ignore.case = TRUE)

  # NA
  expect_error(rank_weight_eval(rw, c(1L, NA_integer_, 3L)), "cannot contain NA", ignore.case = TRUE)

  # Zero length
  expect_error(rank_weight_eval(rw, integer(0)), "positive length", ignore.case = TRUE)

  # Bad spec: neither eval nor eval_batch
  bad_spec <- list(name = "bad")
  bad_spec$.__kind <- "rank_weight"
  class(bad_spec) <- c("rank_weight", class(bad_spec))
  expect_error(rank_weight_eval(bad_spec, c(1L, 2L)), "neither eval nor eval_batch")

  # Bad returns: wrong length in eval
  bad_len <- as_rank_weight(function(r, ctx = NULL) 1, name = "bad_len")
  expect_error(rank_weight_eval(bad_len, c(1L, 2L, 3L)), "length\\(ranks\\)")

  # Bad returns: non-finite
  bad_inf <- as_rank_weight(function(r, ctx = NULL) rep(Inf, length(r)), name = "bad_inf")
  expect_error(rank_weight_eval(bad_inf, c(1L, 2L, 3L)), "finite", ignore.case = TRUE)

  # Batch validator through matrix path: wrong shape returned
  rw_bad_batch_shape <- list(
    name = "bad_batch_shape",
    eval_batch = function(R, ctx = NULL) matrix(0, 1, 1)
  )
  rw_bad_batch_shape$.__kind <- "rank_weight"
  class(rw_bad_batch_shape) <- c("rank_weight", class(rw_bad_batch_shape))
  R <- .rand_ranks_mat(5L, 2L)
  expect_error(rank_weight_eval_batch(rw_bad_batch_shape, R), "nrow\\(R\\) x ncol\\(R\\)")
})
