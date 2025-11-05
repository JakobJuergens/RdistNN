# tests/testthat/test-api-weights.R
# ------------------------------------------------------------------------------

test_that("helpers: deterministic ordering + ranks conversion round-trip", {
  skip_if_not(exists("ord_to_ranks"), "ord_to_ranks() not found in package")

  set.seed(1)
  n <- 6
  p <- 2
  m <- 3
  X  <- matrix(round(rnorm(n * p), 3), n, p)
  X0 <- matrix(round(rnorm(m * p), 3), m, p)

  # Build DX with intentional ties in column 2
  DX <- vapply(seq_len(m), function(j) {
    rowSums((X - matrix(X0[j, ], n, p, byrow = TRUE))^2)
  }, numeric(n))
  DX[, 2] <- c(1, 1, 2, 2, 3, 3) # create ties

  # What api_weights.R would compute from DX
  ORD <- apply(DX, 2L, .order_with_ties)
  RKS <- apply(ORD, 2L, ord_to_ranks)

  # Ties: equal distances order by row index
  expect_equal(ORD[, 2], order(DX[, 2], seq_along(DX[, 2])))

  # Each column is a permutation 1..n
  for (j in seq_len(m)) {
    expect_setequal(ORD[, j], seq_len(n))
    expect_setequal(RKS[, j], seq_len(n))
  }
})

test_that(".resolve_ranking works for dx, ord, ranks, and function outputs", {
  skip_if_not(exists(".resolve_ranking"), ".resolve_ranking() not found")

  n <- 5; p <- 2; m <- 2
  X  <- matrix(seq_len(n * p), n, p)
  X0 <- matrix(c(0, 0, 10, 10), m, p)

  DX <- vapply(seq_len(m), function(j) {
    rowSums((X - matrix(X0[j, ], n, p, byrow = TRUE))^2)
  }, numeric(n))
  ORD <- apply(DX, 2L, .order_with_ties)
  RKS <- apply(ORD, 2L, ord_to_ranks)

  # 1) list(dx=...)
  out_dx <- .resolve_ranking(list(dx = DX), X, X0, NULL)
  expect_equal(out_dx$dx, DX)
  expect_equal(out_dx$ord, ORD)
  expect_equal(out_dx$ranks, RKS)

  # 2) list(ord=...)
  out_ord <- .resolve_ranking(list(ord = ORD), X, X0, NULL)
  expect_null(out_ord$dx)
  expect_equal(out_ord$ord, ORD)
  expect_equal(out_ord$ranks, RKS)

  # 3) list(ranks=...)
  out_rk <- .resolve_ranking(list(ranks = RKS), X, X0, NULL)
  expect_null(out_rk$dx)
  expect_null(out_rk$ord)
  expect_equal(out_rk$ranks, RKS)

  # 4) function returning a list with dx
  ranker_fun <- function(X, X0, ctx) list(dx = DX)
  out_fn <- .resolve_ranking(ranker_fun, X, X0, NULL)
  expect_equal(out_fn$dx, DX)
  expect_equal(out_fn$ord, ORD)
  expect_equal(out_fn$ranks, RKS)
})

test_that(".resolve_ranking rejects malformed shapes", {
  skip_if_not(exists(".resolve_ranking"), ".resolve_ranking() not found")
  n <- 4; p <- 2; m <- 2
  X  <- matrix(0, n, p)
  X0 <- matrix(0, m, p)

  bad_ord <- matrix(1L, n - 1L, m) # wrong n
  expect_error(.resolve_ranking(list(ord = bad_ord), X, X0),
               "`ord` must be n x m", fixed = TRUE)

  bad_ranks <- matrix(1L, n, m)
  bad_ranks[, 1] <- c(1L, 1L, 2L, 3L) # not a permutation
  expect_error(.resolve_ranking(list(ranks = bad_ranks), X, X0),
               "Each column of `ranks` must be a permutation",
               fixed = TRUE)
})

test_that(".resolve_weight identifies rank_weight and zero-arg factories", {
  skip_if_not(exists(".resolve_weight"), ".resolve_weight() not found")
  skip_if_not(exists("as_rank_weight"), "as_rank_weight() not found")

  # A simple rank_weight: w_i ∝ 1 / rank_i
  rw <- as_rank_weight(function(r, ctx = NULL) {
    w <- 1 / as.numeric(r)
    w / sum(w)
  }, name = "inv_rank")

  got <- .resolve_weight(rw)
  expect_equal(got$kind, "rank")
  expect_equal(got$spec$name, "inv_rank")

  got2 <- .resolve_weight(function() rw) # zero-arg factory
  expect_equal(got2$kind, "rank")
  expect_equal(got2$spec$name, "inv_rank")
})

test_that("weights_matrix works with distance-based ranking + rank_weight", {
  skip_if_not(exists("as_rank_weight"), "as_rank_weight() not found")

  n <- 5; p <- 2; m <- 3
  set.seed(42)
  X  <- matrix(rnorm(n * p), n, p)
  X0 <- matrix(rnorm(m * p), m, p)

  # Precompute distances for deterministic test
  DX <- vapply(seq_len(m), function(j) {
    rowSums((X - matrix(X0[j, ], n, p, byrow = TRUE))^2)
  }, numeric(n))

  # Rank-weight: harmonic by ranks
  rw <- as_rank_weight(function(r, ctx = NULL) {
    w <- 1 / as.numeric(r)
    w / sum(w)
  }, name = "inv_rank")

  res <- weights_matrix(
    X = X, X0 = X0,
    ranking = list(dx = DX),  # bypass registry; use precomputed dx
    weight  = rw,
    return  = c("W", "ord", "ranks", "dx")
  )

  # Shapes
  expect_equal(dim(res$W), c(n, m))
  expect_equal(dim(res$ord), c(n, m))
  expect_equal(dim(res$ranks), c(n, m))
  expect_equal(dim(res$dx), c(n, m))

  # Column sums ~ 1
  expect_true(all(abs(colSums(res$W) - 1) < 1e-12))

  # Metadata
  expect_equal(res$meta$n, n)
  expect_equal(res$meta$m, m)
  expect_equal(res$meta$p, p)
  expect_equal(res$meta$weight_kind, "rank")
  expect_equal(res$meta$weight_name, "inv_rank")

  # Deterministic ties: if we create equal dx entries in a column,
  # ordering must break ties by row index.
  DX_tie <- DX
  DX_tie[1:2, 1] <- 7
  DX_tie[3:4, 1] <- 7
  res2 <- weights_matrix(
    X = X, X0 = X0,
    ranking = list(dx = DX_tie),
    weight  = rw,
    return  = "ord"
  )
  expect_equal(res2$ord[, 1], order(DX_tie[, 1], seq_along(DX_tie[, 1])))
})

test_that("weights_matrix handles x0 vs X0 and basic shape errors", {
  skip_if_not(exists("as_rank_weight"), "as_rank_weight() not found")

  n <- 4; p <- 2; m <- 1
  X  <- matrix(1, n, p)
  x0 <- c(0, 0)

  DX <- matrix(rowSums((X - matrix(x0, n, p, byrow = TRUE))^2), n, m)

  rw <- as_rank_weight(function(r, ctx = NULL) {
    w <- 1 / as.numeric(r)
    w / sum(w)
  }, name = "inv_rank")

  # x0 convenience path
  res <- weights_matrix(X = X, x0 = x0, ranking = list(dx = DX), weight = rw)
  expect_equal(dim(res$W), c(n, m))

  # can't pass both X0 and x0
  expect_error(
    weights_matrix(
      X = X, X0 = matrix(x0, 1), x0 = x0,
      ranking = list(dx = DX), weight = rw
    ),
    "Provide either `X0` or `x0`, not both."
  )

  # mismatched ncol
  expect_error(
    weights_matrix(
      X = X, X0 = matrix(0, 2, p + 1),
      ranking = list(dx = DX), weight = rw
    ),
    "ncol\\(X0\\) must equal ncol\\(X\\)"
  )
})

test_that("print.weights_result is informative", {
  n <- 3; m <- 2; p <- 1
  out <- list(
    W = matrix(0, n, m),
    meta = list(n = n, m = m, p = p, weight_kind = "rank", weight_name = "mock")
  )
  class(out) <- "weights_result"
  expect_output(print(out), "<weights_result>")
  expect_output(print(out), "data  : n=3, p=1, m=2")
  expect_output(print(out), "weight: rank \\(mock\\)")
})

test_that("rank-weight path works if available (DNN wrappers)", {
  # Optional integration test; skip if your DNN factory isn't present yet.
  skip_if_not(
    exists("weights_dnn_exact") &&
      exists("rank_weight_dnn_exact") &&
      exists("ord_to_ranks"),
    "DNN rank-weight API not found; skipping."
  )

  set.seed(9)
  n <- 6; p <- 2; m <- 2
  X  <- matrix(rnorm(n * p), n, p)
  X0 <- matrix(rnorm(m * p), m, p)

  res <- weights_dnn_exact(X, X0, distance = "l2", k = 2L, s = 3L, return = c("W", "ranks"))
  expect_equal(dim(res$W), c(n, m))
  expect_true(all(is.finite(res$W)))
  expect_equal(dim(res$ranks), c(n, m))
  expect_true(all(abs(colSums(res$W) - 1) < 1e-8))
})
