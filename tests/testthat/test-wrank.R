# -------------------------------------------------------------------
# rank_weight (S3) — scaffolding tests for rank.R
# Tags:
#   [smoke]    Construction & printing
#   [validate] Validation errors
#   [setup]    Context plumbing
#   [engine]   Eval / Eval-batch shapes, finite checks, input-path errors
#   [inputs]   Rank vector/matrix validation via public APIs
#   [helpers]  as_rank_weight / is_rank_weight / ord_to_ranks
# -------------------------------------------------------------------

test_that("[smoke] can construct a minimal rank_weight and print it", {
  f <- function(ranks, ctx = NULL) {
    # simple normalized 1/r scheme
    r <- as.integer(ranks)
    out <- 1 / r
    out / sum(out)
  }
  rw <- new_rank_weight(name = "inv_rank", eval = f)
  expect_s3_class(rw, "rank_weight")
  expect_output(print(rw), "<rank_weight: ")
  expect_output(print(rw), "impl :\\s*R eval")
})

test_that("[validate] errors on empty name or missing implementation", {
  f <- function(ranks, ctx = NULL) rep(0, length(ranks))
  expect_error(
    new_rank_weight(name = "", eval = f),
    "rank-weight name must be non-empty"
  )
  expect_error(
    new_rank_weight(name = "noimpl"),
    "Provide at least one of: ptr, eval, eval_batch"
  )
})

test_that("[setup] context from rank_weight_setup is threaded into eval and eval_batch", {
  used_ctx_eval <- FALSE
  used_ctx_batch <- FALSE

  rw <- new_rank_weight(
    name = "ctx_rw",
    setup = function(data, scale = 2) list(scale = scale),
    eval = function(ranks, ctx = NULL) {
      used_ctx_eval <<- isTRUE(ctx$scale == 3)
      w <- rep(0, length(ranks))
      w[ranks <= 1] <- 1 # put mass on top-1
      w * ctx$scale # scale, will be renormalized here for testing purposes
    },
    eval_batch = function(R, ctx = NULL) {
      used_ctx_batch <<- isTRUE(ctx$scale == 3)
      W <- matrix(0, nrow(R), ncol(R))
      for (j in seq_len(ncol(R))) W[R[, j] == 1L, j] <- ctx$scale
      W
    }
  )

  n <- 5L
  ranks <- 1:n
  ctx <- rank_weight_setup(rw, data = NULL, scale = 3)
  w <- rank_weight_eval(rw, ranks, ctx)
  expect_true(used_ctx_eval)
  expect_equal(length(w), n)

  R <- cbind(1:n, n:1)
  W <- rank_weight_eval_batch(rw, R, ctx)
  expect_true(used_ctx_batch)
  expect_equal(dim(W), c(n, 2L))
})

test_that("[engine] rank_weight_eval enforces numeric(n) finite; rank_weight_eval_batch enforces size", {
  # good eval
  rw_ok <- new_rank_weight(
    name = "ok",
    eval = function(ranks, ctx = NULL) {
      w <- rep(0, length(ranks))
      w[ranks <= 2] <- 0.5
      w
    }
  )
  n <- 6L
  ranks <- sample.int(n)
  out <- rank_weight_eval(rw_ok, ranks)
  expect_true(is.numeric(out))
  expect_length(out, n)
  expect_true(all(is.finite(out)))

  # bad eval: wrong length
  rw_bad_len <- new_rank_weight(
    name = "bad_len",
    eval = function(ranks, ctx = NULL) rep(1, length(ranks) + 1L)
  )
  expect_error(
    rank_weight_eval(rw_bad_len, ranks),
    "must return numeric\\(length\\(ranks\\)\\) with finite entries"
  )

  # good eval_batch
  rw_b_ok <- new_rank_weight(
    name = "okb",
    eval_batch = function(R, ctx = NULL) {
      W <- matrix(0, nrow(R), ncol(R))
      for (j in seq_len(ncol(R))) W[R[, j] <= 2L, j] <- 0.5
      W
    }
  )
  R <- cbind(sample.int(n), sample.int(n))
  OB <- rank_weight_eval_batch(rw_b_ok, R)
  expect_true(is.matrix(OB))
  expect_equal(dim(OB), dim(R))
  expect_true(all(is.finite(OB)))

  # bad eval_batch: wrong shape
  rw_b_bad <- new_rank_weight(
    name = "badb",
    eval_batch = function(R, ctx = NULL) matrix(1, nrow(R) + 1L, ncol(R))
  )
  expect_error(
    rank_weight_eval_batch(rw_b_bad, R),
    "must return a finite numeric matrix of size nrow\\(R\\) x ncol\\(R\\)"
  )
})

test_that("[engine] path errors when respective R fallback is absent", {
  n <- 5L
  ranks <- sample.int(n)
  R <- cbind(sample.int(n), sample.int(n))

  rw_only_batch <- new_rank_weight(
    name = "only_batch",
    eval_batch = function(R, ctx = NULL) matrix(0, nrow(R), ncol(R))
  )
  expect_error(
    rank_weight_eval(rw_only_batch, ranks),
    "No R `eval` provided; use compiled `ptr` or R `eval_batch`\\."
  )

  rw_only_eval <- new_rank_weight(
    name = "only_eval",
    eval = function(ranks, ctx = NULL) rep(0, length(ranks))
  )
  expect_error(
    rank_weight_eval_batch(rw_only_eval, R),
    "No R `eval_batch` provided; use compiled `ptr` or R `eval`\\."
  )
})

test_that("[inputs] rank vector validation fires via public API", {
  rw <- new_rank_weight(
    name = "v",
    eval = function(ranks, ctx = NULL) rep(1 / length(ranks), length(ranks))
  )
  expect_error(rank_weight_eval(rw, integer(0)), "positive length")
  expect_error(rank_weight_eval(rw, c(1L, NA_integer_)), "cannot contain NA")
  expect_error(rank_weight_eval(rw, c(1L, 1L)), "permutation of 1..n")
})

test_that("[inputs] rank matrix validation fires via public API", {
  rw <- new_rank_weight(
    name = "m",
    eval_batch = function(R, ctx = NULL) matrix(1 / nrow(R), nrow(R), ncol(R))
  )
  expect_error(rank_weight_eval_batch(rw, matrix(integer(0), 0, 0)), "positive nrow and ncol")
  R_bad <- cbind(c(1L, 1L, 2L), c(1L, 2L, NA_integer_))
  expect_error(rank_weight_eval_batch(rw, R_bad), "cannot contain NA|permutation of 1..n")
})

test_that("[helpers] as_rank_weight wraps a function; is_rank_weight predicate; ord_to_ranks round-trip", {
  # as_rank_weight + predicate
  f <- function(ranks, ctx = NULL) rep(1 / length(ranks), length(ranks))
  rw <- as_rank_weight(f, name = "wrapped")
  expect_s3_class(rw, "rank_weight")
  expect_true(is_rank_weight(rw))

  # ord_to_ranks round-trip property: order(ranks) == ord
  set.seed(7)
  n <- 11L
  ord <- sample.int(n)
  ranks <- ord_to_ranks(ord)
  expect_true(all(ranks[ord] == seq_len(n)))
  expect_identical(order(ranks), ord)

  # ord_to_ranks errors
  expect_error(ord_to_ranks(c(1L, 1L, 3L)), "permutation of 1..n")
})
