# tests/testthat/test-api-weights.R
# ------------------------------------------------------------------------------

test_that("helpers: deterministic ordering + ranks conversion round-trip", {
  skip_if_not(exists("ord_to_ranks"), "ord_to_ranks() not found in package")

  set.seed(1)
  n <- 6
  p <- 2
  m <- 3
  X <- matrix(round(rnorm(n * p), 3), n, p)
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

  n <- 5
  p <- 2
  m <- 2
  X <- matrix(seq_len(n * p), n, p)
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
  n <- 4
  p <- 2
  m <- 2
  X <- matrix(0, n, p)
  X0 <- matrix(0, m, p)

  bad_ord <- matrix(1L, n - 1L, m) # wrong n
  expect_error(.resolve_ranking(list(ord = bad_ord), X, X0),
    "`ord` must be n x m",
    fixed = TRUE
  )

  bad_ranks <- matrix(1L, n, m)
  bad_ranks[, 1] <- c(1L, 1L, 2L, 3L) # not a permutation
  expect_error(.resolve_ranking(list(ranks = bad_ranks), X, X0),
    "Each column of `ranks` must be a permutation",
    fixed = TRUE
  )

  bad_dx_shape <- matrix(0, n - 1L, m)
  expect_error(.resolve_ranking(list(dx = bad_dx_shape), X, X0), "`dx` must be n x m", fixed = TRUE)

  bad_dx_na <- matrix(0, n, m)
  bad_dx_na[1L, 1L] <- NA_real_
  expect_error(.resolve_ranking(list(dx = bad_dx_na), X, X0), "dx.*NA", ignore.case = TRUE)

  bad_dx_inf <- matrix(0, n, m)
  bad_dx_inf[1L, 1L] <- Inf
  expect_error(.resolve_ranking(list(dx = bad_dx_inf), X, X0), "dx.*non-finite", ignore.case = TRUE)

  bad_ord_frac <- cbind(seq_len(n), rev(seq_len(n)))
  storage.mode(bad_ord_frac) <- "double"
  bad_ord_frac[1L, 1L] <- 1.5
  expect_error(.resolve_ranking(list(ord = bad_ord_frac), X, X0), "ord.*integer", ignore.case = TRUE)

  bad_ranks_frac <- cbind(seq_len(n), rev(seq_len(n)))
  storage.mode(bad_ranks_frac) <- "double"
  bad_ranks_frac[1L, 1L] <- 1.5
  expect_error(.resolve_ranking(list(ranks = bad_ranks_frac), X, X0), "ranks.*integer", ignore.case = TRUE)

  bad_ranker <- function(X, X0, ctx) {
    out <- matrix(0, nrow(X), nrow(X0))
    out[1L, 1L] <- Inf
    list(dx = out)
  }
  expect_error(.resolve_ranking(bad_ranker, X, X0), "dx.*non-finite", ignore.case = TRUE)
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

  n <- 5
  p <- 2
  m <- 3
  set.seed(42)
  X <- matrix(rnorm(n * p), n, p)
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
    ranking = list(dx = DX), # bypass registry; use precomputed dx
    weight = rw,
    return = c("W", "ord", "ranks", "dx")
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
    weight = rw,
    return = "ord"
  )
  expect_equal(res2$ord[, 1], order(DX_tie[, 1], seq_along(DX_tie[, 1])))
})

test_that("weights_matrix passes ctx through to rank_weight evaluators", {
  n <- 5
  p <- 2
  m <- 2
  set.seed(24)
  X <- matrix(rnorm(n * p), n, p)
  X0 <- matrix(rnorm(m * p), m, p)

  DX <- vapply(seq_len(m), function(j) {
    rowSums((X - matrix(X0[j, ], n, p, byrow = TRUE))^2)
  }, numeric(n))

  rw_ctx <- new_rank_weight(
    name = "exp_ctx",
    eval_batch = function(R, ctx = NULL) {
      lambda <- if (is.null(ctx) || is.null(ctx$lambda)) 0.25 else ctx$lambda
      W <- exp(-lambda * (R - 1L))
      sweep(W, 2L, colSums(W), "/")
    }
  )

  res_lo <- weights_matrix(
    X = X, X0 = X0,
    ranking = list(dx = DX),
    weight = rw_ctx,
    ctx = list(lambda = 0.10),
    return = c("W", "ranks")
  )
  res_hi <- weights_matrix(
    X = X, X0 = X0,
    ranking = list(dx = DX),
    weight = rw_ctx,
    ctx = list(lambda = 1.00),
    return = c("W", "ranks")
  )

  expect_true(any(abs(res_lo$W - res_hi$W) > 1e-12))

  W_manual <- rank_weight_eval_batch(rw_ctx, res_lo$ranks, ctx = list(lambda = 0.10))
  expect_equal(res_lo$W, W_manual, tolerance = 1e-12)
})

test_that("weights_matrix handles x0 vs X0 and basic shape errors", {
  skip_if_not(exists("as_rank_weight"), "as_rank_weight() not found")

  n <- 4
  p <- 2
  m <- 1
  X <- matrix(1, n, p)
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
  n <- 3
  m <- 2
  p <- 1
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
  n <- 6
  p <- 2
  m <- 2
  X <- matrix(rnorm(n * p), n, p)
  X0 <- matrix(rnorm(m * p), m, p)

  res <- weights_dnn_exact(X, X0, distance = "l2", k = 2L, s = 3L, return = c("W", "ranks"))
  expect_equal(dim(res$W), c(n, m))
  expect_true(all(is.finite(res$W)))
  expect_equal(dim(res$ranks), c(n, m))
  expect_true(all(abs(colSums(res$W) - 1) < 1e-8))
})

# ------------------------------------------------------------------------------
# Extra robustness checks for weights API
# ------------------------------------------------------------------------------

test_that("end-to-end parity: dx vs ord vs ranks produce identical W", {
  skip_if_not(exists("as_rank_weight"), "as_rank_weight() not found")

  n <- 7
  p <- 3
  m <- 2
  set.seed(123)
  X <- matrix(rnorm(n * p), n, p)
  X0 <- matrix(rnorm(m * p), m, p)

  # Distances and induced ord/ranks
  DX <- vapply(seq_len(m), function(j) {
    rowSums((X - matrix(X0[j, ], n, p, byrow = TRUE))^2)
  }, numeric(n))
  ORD <- apply(DX, 2L, .order_with_ties)
  RKS <- apply(ORD, 2L, ord_to_ranks)

  # Simple rank weight
  rw <- as_rank_weight(function(r) {
    w <- 1 / as.numeric(r)
    w / sum(w)
  }, name = "inv_rank_parity")

  res_dx <- weights_matrix(X, X0, ranking = list(dx = DX), weight = rw, return = "W")
  res_ord <- weights_matrix(X, X0, ranking = list(ord = ORD), weight = rw, return = "W")
  res_rks <- weights_matrix(X, X0, ranking = list(ranks = RKS), weight = rw, return = "W")

  expect_equal(res_dx$W, res_ord$W, tolerance = 0)
  expect_equal(res_dx$W, res_rks$W, tolerance = 0)
})

test_that("return= filters extras: only requested components are present", {
  skip_if_not(exists("as_rank_weight"), "as_rank_weight() not found")

  n <- 6
  p <- 2
  m <- 2
  set.seed(1)
  X <- matrix(rnorm(n * p), n, p)
  X0 <- matrix(rnorm(m * p), m, p)

  DX <- vapply(seq_len(m), function(j) {
    rowSums((X - matrix(X0[j, ], n, p, byrow = TRUE))^2)
  }, numeric(n))

  rw <- as_rank_weight(function(r) {
    w <- 0.8^(as.numeric(r) - 1L)
    w / sum(w)
  }, name = "geom_0.8")

  res_onlyW <- weights_matrix(X, X0, ranking = list(dx = DX), weight = rw, return = "W")
  expect_true(is.matrix(res_onlyW$W))
  expect_false("ord" %in% names(res_onlyW))
  expect_false("ranks" %in% names(res_onlyW))
  expect_false("dx" %in% names(res_onlyW))

  res_some <- weights_matrix(X, X0, ranking = list(dx = DX), weight = rw, return = c("W", "ranks"))
  expect_true("W" %in% names(res_some))
  expect_true("ranks" %in% names(res_some))
  expect_false("ord" %in% names(res_some)) # not requested
  expect_false("dx" %in% names(res_some)) # not requested
})

test_that("ranking name path works via distance registry", {
  skip_if_not(
    exists("new_distance") && exists("register_distance") &&
      exists("get_distance") && exists("distance_eval_batch") &&
      exists("as_rank_weight"),
    "distance API pieces not found; skipping."
  )

  n <- 5
  p <- 2
  m <- 3
  set.seed(7)
  X <- matrix(rnorm(n * p), n, p)
  X0 <- matrix(rnorm(m * p), m, p)

  # Define a tiny cosine distance for batch only
  cos_spec <- new_distance(
    name = "cosine_weights_api",
    eval_batch = function(X, X0, ctx) {
      nx <- sqrt(rowSums(X^2))
      nx0 <- sqrt(rowSums(X0^2))
      S <- X %*% t(X0)
      D <- 1 - sweep(sweep(S, 1L, nx, `/`), 2L, nx0, `/`)
      pmax(D, 0)
    }
  )
  register_distance("cosine_weights_api", cos_spec)

  # Rank-weight
  rw <- as_rank_weight(function(r) {
    w <- 1 / as.numeric(r)
    w / sum(w)
  }, name = "inv_rank_cos")

  # By-name
  res_name <- weights_matrix(X, X0, ranking = "cosine_weights_api", weight = rw, return = c("W", "ranks", "dx"))
  # Explicit dx (should match)
  DX <- distance_eval_batch(get_distance("cosine_weights_api"), X, X0)
  res_dx <- weights_matrix(X, X0, ranking = list(dx = DX), weight = rw, return = c("W", "ranks", "dx"))

  expect_equal(res_name$dx, res_dx$dx, tolerance = 0)
  expect_equal(res_name$ranks, res_dx$ranks, tolerance = 0)
  expect_equal(res_name$W, res_dx$W, tolerance = 0)
})

test_that("x0 convenience path equals 1-column X0 path", {
  skip_if_not(exists("as_rank_weight"), "as_rank_weight() not found")

  n <- 8
  p <- 3
  set.seed(99)
  X <- matrix(rnorm(n * p), n, p)
  x0 <- rnorm(p)

  DX1 <- matrix(rowSums((X - matrix(x0, n, p, byrow = TRUE))^2), n, 1)

  rw <- as_rank_weight(function(r) {
    w <- 1 / as.numeric(r)
    w / sum(w)
  }, name = "inv_rank_x0")

  res_x0 <- weights_matrix(X = X, x0 = x0, ranking = list(dx = DX1), weight = rw, return = "W")
  res_X01 <- weights_matrix(X = X, X0 = matrix(x0, 1), ranking = list(dx = DX1), weight = rw, return = "W")

  expect_equal(res_x0$W, res_X01$W, tolerance = 0)
})

test_that("metadata types and contents are stable", {
  skip_if_not(exists("as_rank_weight"), "as_rank_weight() not found")

  n <- 5
  p <- 4
  m <- 3
  set.seed(2)
  X <- matrix(rnorm(n * p), n, p)
  X0 <- matrix(rnorm(m * p), m, p)

  DX <- vapply(seq_len(m), function(j) {
    rowSums((X - matrix(X0[j, ], n, p, byrow = TRUE))^2)
  }, numeric(n))

  rw <- as_rank_weight(function(r) {
    w <- 0.7^(as.numeric(r) - 1L)
    w / sum(w)
  }, name = "geom_0.7_meta")

  res <- weights_matrix(X, X0, ranking = list(dx = DX), weight = rw, return = "W")

  expect_true(is.list(res$meta))
  expect_true(is.integer(res$meta$n))
  expect_true(is.integer(res$meta$m))
  expect_true(is.integer(res$meta$p))
  expect_equal(res$meta$n, n)
  expect_equal(res$meta$m, m)
  expect_equal(res$meta$p, p)
  expect_equal(res$meta$weight_kind, "rank")
  expect_equal(res$meta$weight_name, "geom_0.7_meta")
})

test_that("weights_matrix warns on unknown return values", {
  skip_if_not(exists("as_rank_weight"), "as_rank_weight() not found")

  n <- 5L; p <- 2L; m <- 2L
  set.seed(42)
  X  <- matrix(rnorm(n * p), n, p)
  X0 <- matrix(rnorm(m * p), m, p)
  DX <- vapply(seq_len(m), function(j) rowSums((X - matrix(X0[j, ], n, p, byrow = TRUE))^2), numeric(n))
  rw <- as_rank_weight(function(r) { w <- 1 / as.numeric(r); w / sum(w) }, name = "inv_rank_ret")

  expect_warning(
    weights_matrix(X, X0, ranking = list(dx = DX), weight = rw, return = c("W", "typo")),
    "Unknown `return` value"
  )
  # valid values do not warn
  expect_no_warning(
    weights_matrix(X, X0, ranking = list(dx = DX), weight = rw, return = c("W", "ranks", "ord", "dx"))
  )
})

test_that("error paths: missing weight and unsupported ranking type", {
  n <- 4
  p <- 2
  m <- 2
  X <- matrix(0, n, p)
  X0 <- matrix(1, m, p)

  DX <- vapply(seq_len(m), function(j) {
    rowSums((X - matrix(X0[j, ], n, p, byrow = TRUE))^2)
  }, numeric(n))

  # missing weight
  expect_error(
    weights_matrix(X, X0, ranking = list(dx = DX), weight = NULL),
    "`weight` must be a rank_weight",
    fixed = TRUE
  )
  msg_missing_weight <- tryCatch(
    {
      weights_matrix(X, X0, ranking = list(dx = DX), weight = NULL)
      NA_character_
    },
    error = function(e) e$message
  )
  expect_match(msg_missing_weight, "Got:", fixed = TRUE)

  # unsupported ranking type
  rw <- as_rank_weight(function(r) {
    w <- 1 / as.numeric(r)
    w / sum(w)
  }, name = "inv_rank_err")
  expect_error(
    weights_matrix(X, X0, ranking = 42L, weight = rw),
    "Unsupported `ranking` type",
    fixed = TRUE
  )
  msg_bad_ranking <- tryCatch(
    {
      weights_matrix(X, X0, ranking = 42L, weight = rw)
      NA_character_
    },
    error = function(e) e$message
  )
  expect_match(msg_bad_ranking, "Got:", fixed = TRUE)
})

test_that("deterministic ties with massive plateaus still stable", {
  skip_if_not(exists("as_rank_weight"), "as_rank_weight() not found")

  n <- 10
  p <- 2
  m <- 1
  X <- cbind(seq_len(n), seq_len(n))
  x0 <- c(0, 0)

  DX <- matrix(5, n, m) # all equal -> order should be 1:n
  rw <- as_rank_weight(function(r) {
    w <- 1 / as.numeric(r)
    w / sum(w)
  }, name = "inv_rank_plateau")

  res <- weights_matrix(X, x0 = x0, ranking = list(dx = DX), weight = rw, return = c("ord", "ranks", "W"))
  expect_equal(as.integer(res$ord[, 1]), seq_len(n))
  expect_setequal(res$ranks[, 1], seq_len(n))
  expect_true(abs(sum(res$W[, 1]) - 1) < 1e-12)
})

test_that("as_rank_weight supports 1-arg and 2-arg functions", {
  n <- 5L
  r <- 1:n

  rw1 <- as_rank_weight(function(r) {
    w <- 1 / r
    w / sum(w)
  }, name = "one_arg")
  rw2 <- as_rank_weight(function(r, ctx = NULL) {
    w <- 0.8^(r - 1L)
    w / sum(w)
  }, name = "two_arg")

  # single target
  w1 <- rank_weight_eval(rw1, r)
  w2 <- rank_weight_eval(rw2, r)
  expect_equal(length(w1), n)
  expect_equal(length(w2), n)
  expect_true(all(is.finite(w1)) && abs(sum(w1) - 1) < 1e-12)
  expect_true(all(is.finite(w2)) && abs(sum(w2) - 1) < 1e-12)

  # batch path still works
  R <- cbind(r, rev(r))
  W1 <- rank_weight_eval_batch(rw1, R)
  W2 <- rank_weight_eval_batch(rw2, R)
  expect_equal(dim(W1), c(n, 2))
  expect_equal(dim(W2), c(n, 2))
  expect_true(all(abs(colSums(W1) - 1) < 1e-12))
  expect_true(all(abs(colSums(W2) - 1) < 1e-12))
})
