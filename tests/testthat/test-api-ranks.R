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

  n <- 8L
  m <- 3L
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

  n <- 9L
  m <- 4L
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

  n <- 6L
  m <- 2L
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
  register_rank(spec_legacy) # legacy one-arg
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

test_that("as_rank_weight handles functions without ctx (and ignores provided ctx)", {
  skip_if_not(exists("as_rank_weight"), "as_rank_weight() not found")
  skip_if_not(exists("rank_weight_eval"), "rank_weight_eval() not found")

  # fun has signature fun(r) with no ctx
  rw_nocontext <- as_rank_weight(function(r) {
    w <- 1 / as.numeric(r)
    w / sum(w)
  }, name = "no_ctx")

  n <- 7L
  r <- sample.int(n)
  # Providing ctx must not error and must be ignored by adapter
  w <- rank_weight_eval(rw_nocontext, r, ctx = list(alpha = 99))
  expect_length(w, n)
  expect_true(all(is.finite(w)))
  expect_equal(sum(w), 1, tolerance = 1e-12)
})

test_that("batch fallback produces an n x 1 matrix when m == 1 (eval-only spec)", {
  skip_if_not(exists("as_rank_weight"), "as_rank_weight() not found")
  skip_if_not(exists("rank_weight_eval_batch"), "rank_weight_eval_batch() not found")

  rw_eval_only <- as_rank_weight(function(r) {
    w <- 1 / as.numeric(r)
    w / sum(w)
  }, name = "eval_only_singlecol")

  n <- 5L
  m <- 1L
  R <- .rand_ranks_mat(n, m)

  W <- rank_weight_eval_batch(rw_eval_only, R)
  expect_true(is.matrix(W))
  expect_equal(dim(W), c(n, m))
  expect_equal(unname(colSums(W)), rep(1, m), tolerance = 1e-12)
})

test_that("when both eval and eval_batch exist, eval_batch is preferred (single and batch)", {
  skip_if_not(exists("rank_weight_eval"), "rank_weight_eval() not found")
  skip_if_not(exists("rank_weight_eval_batch"), "rank_weight_eval_batch() not found")

  calls <- new.env(parent = emptyenv())
  calls$eval <- 0L
  calls$batch <- 0L

  spec <- list(
    name = "both_paths",
    eval = function(r, ctx = NULL) {
      calls$eval <- calls$eval + 1L
      # Should not be called if eval_batch exists
      rep.int(NA_real_, length(r))
    },
    eval_batch = function(R, ctx = NULL) {
      calls$batch <- calls$batch + 1L
      # pmax() drops dims; rebuild the matrix shape explicitly
      W <- 1 / pmax(1, as.numeric(R))
      W <- matrix(W, nrow = nrow(R), ncol = ncol(R))
      # or: dim(W) <- dim(R)
      W <- sweep(W, 2, colSums(W), "/")
      W
    }
  )
  spec$.__kind <- "rank_weight"
  class(spec) <- c("rank_weight", class(spec))

  n <- 6L
  m <- 3L
  R <- .rand_ranks_mat(n, m)

  # Batch path should trigger eval_batch
  W <- rank_weight_eval_batch(spec, R)
  expect_equal(dim(W), c(n, m))
  expect_equal(calls$batch, 1L)
  expect_equal(calls$eval, 0L)

  # Single-target path should still prefer eval_batch
  w1 <- rank_weight_eval(spec, R[, 2L])
  expect_length(w1, n)
  expect_equal(calls$batch, 2L) # called again
  expect_equal(calls$eval, 0L)
})

test_that("register_rank overwrite is idempotent and replaces behavior", {
  skip_if_not(exists("register_rank") && exists("get_rank") && exists("rank_weight_eval"))

  nm <- paste0("rw_over_", as.integer(runif(1, 1e6, 2e6)))
  spec1 <- function(r) {
    w <- 1 / as.numeric(r)
    w / sum(w)
  }
  spec2 <- function(r) {
    w <- 1 / (as.numeric(r)^2)
    w / sum(w)
  }

  register_rank(nm, spec1)
  r <- sample.int(8L)
  w1 <- rank_weight_eval(get_rank(nm), r)

  # Overwrite with different behavior
  register_rank(nm, spec2)
  w2 <- rank_weight_eval(get_rank(nm), r)

  expect_false(isTRUE(all.equal(w1, w2)))
  expect_equal(sum(w2), 1, tolerance = 1e-12)
})

test_that("rank_weight_eval_batch validates columns are permutations", {
  skip_if_not(exists("rank_weight_eval_batch") && exists("as_rank_weight"))

  rw <- as_rank_weight(function(r) (1 / as.numeric(r)) / sum(1 / as.numeric(r)), "perm_check")
  n <- 6L
  m <- 2L
  R <- .rand_ranks_mat(n, m)
  R[, 2L] <- c(1L, 1L, 3L, 4L, 5L, 6L) # not a permutation
  expect_error(rank_weight_eval_batch(rw, R), "permutation", ignore.case = TRUE)
})

test_that("rank_weight_eval rejects NA or NaN returns from user functions", {
  skip_if_not(exists("as_rank_weight") && exists("rank_weight_eval"))

  rw_na <- as_rank_weight(function(r) rep(NA_real_, length(r)), "bad_na")
  expect_error(rank_weight_eval(rw_na, 1:4), "finite", ignore.case = TRUE)

  rw_nan <- as_rank_weight(function(r) rep(NaN, length(r)), "bad_nan")
  expect_error(rank_weight_eval(rw_nan, 1:4), "finite", ignore.case = TRUE)
})

test_that("integerish numeric ranks are accepted; non-integer permutations rejected", {
  skip_if_not(exists("rank_weight_eval") && exists("as_rank_weight"))

  rw <- as_rank_weight(function(r) (1 / as.numeric(r)) / sum(1 / as.numeric(r)), "intish")

  # numeric but integerish -> OK
  r_num <- as.numeric(1:5)
  expect_no_error(rank_weight_eval(rw, r_num))

  # non-integer values -> becomes non-permutation after as.integer()
  r_bad <- c(1, 2, 2.2, 4, 5)
  expect_error(rank_weight_eval(rw, r_bad), "permutation", ignore.case = TRUE)
})

test_that("n = 1 trivial case works (both eval and eval_batch paths)", {
  skip_if_not(exists("rank_weight_eval") && exists("rank_weight_eval_batch") && exists("as_rank_weight"))

  rw_eval <- as_rank_weight(function(r) 1, "trivial_eval")
  w <- rank_weight_eval(rw_eval, 1L)
  expect_equal(w, 1)

  rw_batch <- list(
    name = "trivial_batch",
    eval_batch = function(R, ctx = NULL) matrix(1, nrow = 1L, ncol = ncol(R))
  )
  rw_batch$.__kind <- "rank_weight"
  class(rw_batch) <- c("rank_weight", class(rw_batch))
  W <- rank_weight_eval_batch(rw_batch, matrix(1L, 1L, 1L))
  expect_equal(dim(W), c(1L, 1L))
  expect_equal(W[1, 1], 1)
})

test_that("is_rank_weight recognizes list specs flagged with __kind", {
  skip_if_not(exists("is_rank_weight"))

  rw_list <- list(
    name = "kind_list",
    eval = function(r, ctx = NULL) {
      w <- 1 / as.numeric(r)
      w / sum(w)
    }
  )
  rw_list$.__kind <- "rank_weight"
  class(rw_list) <- c("rank_weight", class(rw_list))
  expect_true(is_rank_weight(rw_list))
})
