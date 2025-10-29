# -------------------------------------------------------------------
# Distance Spec (S3) — core behavior tests
# Tags:
#   [smoke]   Construction & printing
#   [validate]Validation errors
#   [setup]   Setup/ctx threading
#   [engine]  Eval shapes/finite & input checks
#   [batch]   Batch eval parity (columnwise)
# -------------------------------------------------------------------

test_that("[smoke] can construct a minimal R distance and print it", {
  dfun <- function(X, x0, ctx = NULL) rowSums(abs(X - matrix(x0, nrow(X), ncol(X), TRUE)))
  d <- new_distance(name = "abs_sum", eval = dfun)
  expect_s3_class(d, "dnn_distance")
  expect_output(print(d), "<dnn_distance: abs_sum>")
})

test_that("[validate] new_distance errors on empty name or no implementation", {
  dfun <- function(X, x0, ctx = NULL) rep(0, nrow(X))
  expect_error(new_distance(name = "", eval = dfun), "distance name must be non-empty")
  expect_error(new_distance(name = "no_impl"), "Provide at least one of: ptr, eval, eval_batch")
})

test_that("[setup] context from distance_setup is threaded into eval", {
  called <- FALSE
  d <- new_distance(
    name = "ctx_check",
    setup = function(data, flag = TRUE) list(flag = flag),
    eval  = function(X, x0, ctx = NULL) {
      stopifnot(is.list(ctx), isTRUE(ctx$flag))
      called <<- TRUE
      rep(0, nrow(X))
    }
  )
  X  <- matrix(0, 3, 2); x0 <- c(1, 1)
  ctx <- distance_setup(d, X, flag = TRUE)
  out <- distance_eval(d, X, x0, ctx)
  expect_true(called)
  expect_equal(out, rep(0, 3))
})

test_that("[engine] distance_eval returns numeric(n) finite; input checks fire", {
  d <- new_distance(
    name = "eucl_R",
    eval = function(X, x0, ctx = NULL) {
      dx2 <- rowSums((X - matrix(x0, nrow(X), ncol(X), TRUE))^2)
      sqrt(dx2)
    }
  )
  X  <- matrix(rnorm(12), 6, 2)
  x0 <- rnorm(2)
  out <- distance_eval(d, X, x0)
  expect_type(out, "double")
  expect_length(out, nrow(X))
  expect_true(all(is.finite(out)))

  Xbad <- matrix(numeric(), 5, 0)
  expect_error(distance_eval(d, Xbad, numeric(0)), "at least one column")
  expect_error(distance_eval(d, X, x0[1]), "x0 length must equal ncol")
})

test_that("[batch] distance_eval_batch shape & parity with columnwise eval", {
  dev <- function(X, x0, ctx = NULL) rowSums(abs(X - matrix(x0, nrow(X), ncol(X), TRUE)))
  db  <- function(X, X0, ctx = NULL) {
    X  <- as.matrix(X); X0 <- as.matrix(X0)
    out <- matrix(NA_real_, nrow(X), nrow(X0))
    for (j in seq_len(nrow(X0))) out[, j] <- dev(X, X0[j, ], ctx)
    out
  }
  d <- new_distance(name = "abs_sum_batch", eval = dev, eval_batch = db)

  X  <- matrix(rnorm(20), 10, 2)
  X0 <- matrix(rnorm(6), 3, 2)
  M  <- distance_eval_batch(d, X, X0)

  expect_s3_class(d, "dnn_distance")
  expect_true(is.matrix(M))
  expect_equal(dim(M), c(nrow(X), nrow(X0)))
  # parity with columnwise eval
  M2 <- sapply(seq_len(nrow(X0)), function(j) distance_eval(d, X, X0[j, ]))
  expect_equal(M, M2, tolerance = 0)
  # dimension mismatch
  expect_error(distance_eval_batch(d, X, matrix(1, 1, 3)), "ncol\\(X0\\) must match ncol\\(X\\)")
})
