# tests/testthat/test-fast-order.R
# -------------------------------------------------------------------
# Fast ordering/ranking: rustat_argsort_dist(), rustat_ranks_from_dist(),
# rustat_topk_dist()
#
# Tags:
#   [parity]  Parity with base R on finite, nonnegative inputs
#   [ties]    Stability (ties broken by index)
#   [inverse] ranks/order inverse property
#   [topk]    Top-k equals head(order, k)
#   [edge]    Edge cases (n=0/1, k out of range)
#   [nf]      Non-finite policy: finite first; then -Inf < +Inf < NaN
#   [mono]    Monotone transform invariance
#   [det]     Determinism
# -------------------------------------------------------------------

test_that("[parity] argsort matches base::order(method='radix') on finite, nonnegative", {
  set.seed(123)
  n <- 5000L
  x <- runif(n) * 10
  ord_fast <- rustat_argsort_dist(x)
  ord_R <- order(x, method = "radix")
  expect_identical(ord_fast, ord_R)

  # include ties
  set.seed(42)
  x2 <- sample(0:50, n, replace = TRUE) / 10
  ord_fast2 <- rustat_argsort_dist(x2)
  ord_R2 <- order(x2, method = "radix")
  expect_identical(ord_fast2, ord_R2)
})

test_that("[ties] stability: equal keys keep increasing original index", {
  x <- rep(7, 100L)
  ord <- rustat_argsort_dist(x)
  expect_identical(ord, seq_along(x)) # stable
  rnk <- rustat_ranks_from_dist(x)
  expect_identical(rnk, seq_along(x)) # ties.method = 'first'
})

test_that("[inverse] ranks and argsort are inverses", {
  set.seed(1)
  x <- runif(1000)
  ord <- rustat_argsort_dist(x)
  rnk <- rustat_ranks_from_dist(x)

  # rnk[ ord[i] ] == i
  expect_true(all(rnk[ord] == seq_along(ord)))
  # ord[ rnk[i] ] == i
  expect_true(all(ord[rnk] == seq_along(rnk)))
})

test_that("[topk] top-k equals head(order, k) (with deterministic tiebreak)", {
  set.seed(99)
  n <- 3000L
  x <- sample(0:1000, n, replace = TRUE) / 10 # many ties
  for (k in c(0L, 1L, 2L, 5L, 17L, 100L, 999L, n)) {
    idx_k <- rustat_topk_dist(x, k)
    ord <- rustat_argsort_dist(x)
    expect_identical(idx_k, head(ord, k))
  }
})

test_that("[edge] length 0/1 and k out of bounds", {
  x0 <- numeric(0)
  expect_identical(rustat_argsort_dist(x0), integer(0))
  expect_identical(rustat_ranks_from_dist(x0), integer(0))
  expect_identical(rustat_topk_dist(x0, 0), integer(0))

  x1 <- 3.14
  expect_identical(rustat_argsort_dist(x1), 1L)
  expect_identical(rustat_ranks_from_dist(x1), 1L)
  expect_identical(rustat_topk_dist(x1, 5), 1L) # clamped to n
  expect_identical(rustat_topk_dist(x1, -3), integer(0)) # clamped to 0
})

test_that("[nf] non-finite policy: finite first; then -Inf < +Inf < NaN; ties by index", {
  x <- c(NaN, -Inf, Inf, 1, -1, 0, NaN, Inf, -Inf, 2, 2, 0)
  n <- length(x)

  # Expected order computed in pure R to match the C++ comparator semantics.
  finite <- is.finite(x)
  # Key1: finite first (0 for finite, 1 for non-finite)
  key1 <- as.integer(!finite)
  # Key2: numeric ascending for finite; 0 sentinel for non-finite (refined by key3)
  key2 <- ifelse(finite, x, 0)
  # Key3: for non-finite only: -Inf (0) < +Inf (1) < NaN (2); 0 for finite
  key3 <- ifelse(finite, 0L,
    ifelse(is.nan(x), 2L,
      ifelse(x > 0, 1L, 0L)
    )
  )
  # Key4: index tiebreak
  key4 <- seq_len(n)

  ord_expected <- order(key1, key2, key3, key4, method = "radix")
  ord_fast <- rustat_argsort_dist(x)

  expect_identical(ord_fast, ord_expected)

  # Ranks consistent with that order
  rnk <- rustat_ranks_from_dist(x)
  expect_true(all(rnk[ord_fast] == seq_along(ord_fast)))
})

test_that("[mono] monotone transforms preserve order (add & positive scale)", {
  set.seed(777)
  x <- abs(rnorm(5000)) # nonnegative
  ord <- rustat_argsort_dist(x)

  y1 <- x + 5
  y2 <- x * 7.3
  expect_identical(rustat_argsort_dist(y1), ord)
  expect_identical(rustat_argsort_dist(y2), ord)
})

test_that("[det] determinism across repeated runs", {
  set.seed(2024)
  x <- runif(10000)
  o1 <- rustat_argsort_dist(x)
  o2 <- rustat_argsort_dist(x)
  r1 <- rustat_ranks_from_dist(x)
  r2 <- rustat_ranks_from_dist(x)
  expect_identical(o1, o2)
  expect_identical(r1, r2)
})

test_that("signed zero is an input-order tie in compiled sorting", {
  cpp <- .rdistnn_compiled_fun("rustat_argsort_dist")
  expect_true(is.function(cpp))
  expect_identical(cpp(c(-0, 0, 1, -0)), c(1L, 2L, 4L, 3L))
})
