# tests/testthat/test-pipeline-e2e.R
# ------------------------------------------------------------------------------

test_that("end-to-end: custom distance + custom rank-weight + weights + score", {
  # --- deterministic toy data ---
  n <- 6
  p <- 3
  m <- 2
  X <- matrix(
    c(
      1, 2, 3,
      4, 5, 6,
      7, 8, 9,
      10, 11, 12,
      13, 14, 15,
      16, 17, 18
    ),
    nrow = n, ncol = p, byrow = TRUE
  )
  X0 <- rbind(
    c(0, 1, 0),
    c(1, 0, 1)
  )
  Y <- c(1, -1, 2, 0, -2, 3)

  # --- user distance: cosine with optional precompute via $setup(X) ---
  cosine_spec <- list(
    name = "cosine_test_e2e",
    setup = function(X) list(nx = sqrt(rowSums(X^2))),
    eval_batch = function(X, X0, ctx) {
      nx <- if (!is.null(ctx$nx)) ctx$nx else sqrt(rowSums(X^2))
      nx0 <- sqrt(rowSums(X0^2))
      S <- X %*% t(X0) # n x m
      D <- 1 - sweep(sweep(S, 1L, nx, `/`), 2L, nx0, `/`)
      pmax(D, 0)
    }
  )
  class(cosine_spec) <- c("dnn_distance", "list") # optional, helps registry

  # No distance_setup() generic in this API; when ctx=NULL the evaluators call spec$setup(X).
  ctx_dist <- cosine_spec$setup(X)

  # --- user rank-weight: exponential decay in rank ---
  lambda <- 0.25
  rw <- as_rank_weight(function(r, ctx = NULL) {
    w <- exp(-lambda * (r - 1L))
    w / sum(w)
  }, name = "exp_decay_e2e")

  # --- compute weights (+ ranks/ord/dx) via top-level API ---
  Wres <- weights_matrix(
    X, X0,
    ranking = cosine_spec, # pass the spec directly
    weight = rw,
    ctx = ctx_dist, # or NULL; both paths tested later
    return = c("W", "ranks", "ord", "dx")
  )

  # structure checks
  expect_equal(dim(Wres$W), c(n, m))
  expect_equal(dim(Wres$ranks), c(n, m))
  expect_equal(dim(Wres$ord), c(n, m))
  expect_equal(dim(Wres$dx), c(n, m))
  expect_true(all(is.finite(Wres$W)))
  expect_true(all(abs(colSums(Wres$W) - 1) < 1e-12))

  for (j in seq_len(m)) {
    expect_setequal(Wres$ord[, j], seq_len(n))
    expect_setequal(Wres$ranks[, j], seq_len(n))
  }

  # recompute W from the returned ranks via the rank-weight spec
  W_manual <- rank_weight_eval_batch(rw, Wres$ranks)
  expect_equal(W_manual, Wres$W, tolerance = 0)

  # pure-score: theta(x0_q) = sum_i W_{iq} * g_i
  g <- Y
  theta_hat <- as.numeric(t(Wres$W) %*% g)
  theta_manual <- colSums(Wres$W * matrix(g, nrow = n, ncol = m))
  expect_equal(theta_hat, theta_manual, tolerance = 0)

  # tiny order-1 score to sanity-check eval path
  my_score <- new_score(
    name  = "pure_score_e2e",
    order = 1L,
    setup = function(data) list(g = as.numeric(data$Y)),
    eval  = function(idxs, data, ctx, rng) sum(ctx$g[idxs])
  )
  df <- data.frame(Y = Y)
  ctx_score <- score_setup(my_score, df)
  expect_equal(score_eval(my_score, 3L, df, ctx_score), Y[3])
  expect_equal(score_eval(my_score, 6L, df, ctx_score), Y[6])
})

test_that("end-to-end robustness: ctx invariance, registration-by-name, x0 path, tie handling, scaling, lambda sensitivity, meta types", {
  n <- 6
  p <- 3
  m <- 2
  X <- matrix(
    c(
      1, 2, 3,
      4, 5, 6,
      7, 8, 9,
      10, 11, 12,
      13, 14, 15,
      16, 17, 18
    ),
    nrow = n, ncol = p, byrow = TRUE
  )
  X0 <- rbind(
    c(0, 1, 0),
    c(1, 0, 1)
  )

  cosine_spec <- list(
    name = "cosine_test_e2e",
    setup = function(X) list(nx = sqrt(rowSums(X^2))),
    eval_batch = function(X, X0, ctx) {
      nx <- if (!is.null(ctx$nx)) ctx$nx else sqrt(rowSums(X^2))
      nx0 <- sqrt(rowSums(X0^2))
      S <- X %*% t(X0)
      D <- 1 - sweep(sweep(S, 1L, nx, `/`), 2L, nx0, `/`)
      pmax(D, 0)
    }
  )
  class(cosine_spec) <- c("dnn_distance", "list")
  ctx_dist <- cosine_spec$setup(X)

  make_rw <- function(lambda) {
    as_rank_weight(function(r, ctx = NULL) {
      w <- exp(-lambda * (r - 1L))
      w / sum(w)
    }, name = paste0("exp_decay_", format(lambda)))
  }
  rw <- make_rw(0.25)
  rw2 <- make_rw(1.00)

  # (1) With ctx vs without ctx should match (setup auto-runs when ctx=NULL)
  W_with_ctx <- weights_matrix(
    X, X0,
    ranking = cosine_spec, weight = rw, ctx = ctx_dist,
    return = c("W", "ranks", "ord", "dx")
  )
  W_no_ctx <- weights_matrix(
    X, X0,
    ranking = cosine_spec, weight = rw, ctx = NULL,
    return = c("W", "ranks", "ord", "dx")
  )
  expect_equal(W_with_ctx$dx, W_no_ctx$dx, tolerance = 0)
  expect_equal(W_with_ctx$ranks, W_no_ctx$ranks, tolerance = 0)
  expect_equal(W_with_ctx$W, W_no_ctx$W, tolerance = 0)

  # (2) Registration-by-name path equals passing spec directly
  register_distance("cosine_test_e2e", cosine_spec)
  register_rank("exp_decay_0.25", rw)
  W_by_name <- weights_matrix(
    X, X0,
    ranking = "cosine_test_e2e",
    weight = get_rank("exp_decay_0.25"),
    return = c("W", "ranks", "ord", "dx")
  )
  expect_equal(W_by_name$dx, W_with_ctx$dx, tolerance = 0)
  expect_equal(W_by_name$ranks, W_with_ctx$ranks, tolerance = 0)
  expect_equal(W_by_name$W, W_with_ctx$W, tolerance = 0)

  # (3) x0 convenience path equals col 1 of multi-target
  x0 <- X0[1, ]
  W_x0 <- weights_matrix(X, x0 = x0, ranking = cosine_spec, weight = rw, ctx = ctx_dist, return = "W")
  expect_equal(W_x0$W[, 1], W_with_ctx$W[, 1], tolerance = 0)

  # (4) Deterministic tie handling
  DX <- W_with_ctx$dx
  DX_tie <- DX
  DX_tie[1, 1] <- DX_tie[2, 1] <- 7
  DX_tie[3, 1] <- DX_tie[4, 1] <- 7
  W_tie <- weights_matrix(
    X, X0,
    ranking = list(dx = DX_tie), weight = rw, return = c("ord", "ranks", "W")
  )
  expect_equal(W_tie$ord[, 1], order(DX_tie[, 1], seq_along(DX_tie[, 1])))
  expect_equal(rank_weight_eval_batch(rw, W_tie$ranks), W_tie$W, tolerance = 0)

  # (5) Cosine scaling invariance (cosine distance is scale-free)
  Xs <- 10 * X
  X0s <- 0.5 * X0
  W_scaled <- weights_matrix(
    Xs, X0s,
    ranking = cosine_spec, weight = rw, ctx = NULL, return = c("dx", "ranks", "W")
  )
  # note: due to floating arithmetic, allow tiny slack
  expect_equal(W_scaled$dx, W_with_ctx$dx, tolerance = 1e-12)
  expect_equal(W_scaled$ranks, W_with_ctx$ranks, tolerance = 0)
  expect_equal(W_scaled$W, W_with_ctx$W, tolerance = 1e-12)

  # (6) Lambda sensitivity (weights change if lambda changes)
  W_rw2 <- weights_matrix(X, X0, ranking = cosine_spec, weight = rw2, ctx = ctx_dist, return = "W")
  expect_true(any(abs(W_rw2$W - W_with_ctx$W) > 1e-12))

  # (7) Meta types and contents
  expect_true(is.integer(W_with_ctx$meta$n))
  expect_true(is.integer(W_with_ctx$meta$m))
  expect_true(is.integer(W_with_ctx$meta$p))
  expect_equal(W_with_ctx$meta$weight_kind, "rank")
  expect_equal(W_with_ctx$meta$weight_name, rw$name)
})

test_that("robustness++: single-target (m=1), degenerate duplicates, ranker as function, large ties, and order-s>1 scores", {
  # --- single-target path (m=1) ---
  n <- 5
  p <- 2
  X <- cbind(1:n, 11:15)
  x0 <- c(1, 11)
  X0 <- matrix(x0, nrow = 1L)

  cosine_spec <- list(
    name = "cosine_test_single",
    eval_batch = function(X, X0, ctx) {
      nx <- sqrt(rowSums(X^2))
      nx0 <- sqrt(rowSums(X0^2))
      S <- X %*% t(X0)
      D <- 1 - sweep(sweep(S, 1L, nx, `/`), 2L, nx0, `/`)
      pmax(D, 0)
    }
  )
  class(cosine_spec) <- c("dnn_distance", "list")

  rw <- as_rank_weight(function(r, ctx = NULL) {
    w <- 0.9^(r - 1L)
    w / sum(w)
  }, name = "geom_0.9")

  W1 <- weights_matrix(X, x0 = x0, ranking = cosine_spec, weight = rw, return = c("W", "ranks"))
  expect_equal(dim(W1$W), c(n, 1))
  expect_setequal(W1$ranks[, 1], seq_len(n))
  expect_true(abs(sum(W1$W[, 1]) - 1) < 1e-12)

  # --- degenerate duplicate rows in X (ties galore) ---
  X_dup <- rbind(X, matrix(X[1, ], 2, p, byrow = TRUE)) # add two exact duplicates
  X0_dup <- matrix(x0, nrow = 1L)

  Wdup <- weights_matrix(X_dup, x0 = x0, ranking = cosine_spec, weight = rw, return = c("ord", "ranks", "dx"))

  # stable by index within ties
  ord_dup <- Wdup$ord
  if (is.null(ord_dup) && !is.null(Wdup$ranks)) {
    # ranks -> order for comparison
    ord_dup <- apply(Wdup$ranks, 2L, order)
  }

  # recompute cosine distances to form the tie-broken order we expect (m=1)
  DX_dup_cos <- distance_eval_batch(cosine_spec, X_dup, X0_dup, ctx = NULL)
  ord_tie <- order(as.vector(DX_dup_cos[, 1]), seq_len(nrow(X_dup)))

  expect_equal(ord_dup[, 1], ord_tie)

  # --- ranker provided as a function returning dx (cosine again) ---
  ranker_fun <- function(X, X0, ctx) {
    D <- distance_eval_batch(cosine_spec, X, X0, ctx = NULL)
    list(dx = D)
  }
  Wfun <- weights_matrix(X, x0 = x0, ranking = ranker_fun, weight = rw, return = c("W", "ranks"))
  expect_equal(dim(Wfun$W), c(n, 1))
  expect_equal(Wfun$ranks, W1$ranks)

  # --- extreme large ties in dx column to stress tie-breaker ---
  DX_all_equal <- matrix(7, nrow = n, ncol = 1) # all equal distances
  Wall <- weights_matrix(X, x0 = x0, ranking = list(dx = DX_all_equal), weight = rw, return = "ord")
  expect_equal(as.integer(Wall$ord[, 1]), seq_len(n)) # ties broken by row index

  # --- order-s > 1 score: s=3 subset-sum kernel using precomputed g in ctx ---
  s <- 3L
  Y <- 1:n
  sc3 <- new_score(
    name  = "sum3",
    order = s,
    setup = function(data) list(g = as.numeric(data$Y)),
    eval  = function(idxs, data, ctx, rng) sum(ctx$g[idxs])
  )
  df <- data.frame(Y = Y)
  ctx3 <- score_setup(sc3, df)
  subs <- rbind(
    c(1L, 2L, 3L),
    c(2L, 4L, 5L),
    c(3L, 4L, 5L)
  )
  vals <- score_eval_batch(sc3, subs, df, ctx3)
  expect_equal(vals, c(sum(1:3), 2 + 4 + 5, 3 + 4 + 5))
})


test_that("robustness++: distance_eval fallback, rank_weight_eval fallback, and registry idempotence", {
  # distance with only eval_batch -> distance_eval() fallback path
  dist_batch_only <- list(
    name = "dist_batch_only_e2e",
    eval_batch = function(X, X0, ctx) {
      vapply(seq_len(nrow(X0)), function(j) {
        rowSums((X - matrix(X0[j, ], nrow(X), ncol(X), byrow = TRUE))^2)
      }, numeric(nrow(X)))
    }
  )
  class(dist_batch_only) <- c("dnn_distance", "list")

  X <- matrix(1:6, 3, 2)
  x0 <- c(1, 2)
  d1 <- distance_eval(dist_batch_only, X, x0)
  d2 <- as.vector(distance_eval_batch(dist_batch_only, X, matrix(x0, 1), NULL))
  expect_equal(d1, d2)

  # rank_weight with only eval_batch -> rank_weight_eval() fallback path
  rw_batch_only <- {
    eval_batch <- function(R, ctx = NULL) {
      n <- nrow(R)
      m <- ncol(R)
      W <- matrix(0, n, m)
      for (j in seq_len(m)) {
        w <- 0.8^(R[, j] - 1L)
        W[, j] <- w / sum(w)
      }
      W
    }
    structure(list(name = "rw_batch_only", eval_batch = eval_batch, .__kind = "rank_weight"),
      class = c("rank_weight", "list")
    )
  }
  r <- 1:5
  w1 <- rank_weight_eval(rw_batch_only, r)
  w2 <- rank_weight_eval_batch(rw_batch_only, matrix(r, ncol = 1))[, 1]
  expect_equal(w1, w2)

  # registries: distance and rank idempotence (last wins, names stable)
  before_d <- sort(distance_names())
  register_distance("dist_batch_only_e2e", dist_batch_only)
  register_distance("dist_batch_only_e2e", dist_batch_only) # again
  after_d <- sort(distance_names())
  expect_true(setequal(after_d, union(before_d, c("dist_batch_only_e2e"))))

  before_r <- sort(rank_names())
  register_rank("rw_batch_only", rw_batch_only)
  register_rank("rw_batch_only", rw_batch_only) # again
  after_r <- sort(rank_names())
  expect_true(setequal(after_r, union(before_r, c("rw_batch_only"))))
})
