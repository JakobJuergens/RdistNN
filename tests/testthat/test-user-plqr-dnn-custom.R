skip_if_not(
  all(vapply(
    c("plqr_dnn", "weights_dnn", "as_rank_weight"),
    exists, logical(1),
    mode = "function"
  )),
  "plqr_dnn() or required helpers not found; skipping file"
)

test_that("custom distance works in plqr_dnn; custom rank_weight works via weights_matrix", {
  skip_if_not(
    all(vapply(
      c(
        "register_distance", "get_distance",
        "weights_dnn", "weights_matrix",
        "as_rank_weight", "rank_weight_eval_batch",
        "distance_eval_batch", "ord_to_ranks"
      ),
      exists, logical(1),
      mode = "function"
    )),
    "required APIs not available"
  )

  set.seed(124)
  d <- .make_plqr_data(n = 40L, p = 3L, p_d = 2L, m = 4L, tau = 0.5)

  nm <- paste0("linf_plqr_", as.integer(runif(1, 1e6, 2e6)))
  register_distance(nm, function(X, x0, ctx) {
    apply(abs(X - matrix(x0, nrow(X), ncol(X), TRUE)), 1L, max)
  })
  spec <- get_distance(nm)

  wt <- weights_dnn(
    X = d$X, X0 = d$X0, distance = nm,
    k = 3L, s = 6L, method = "exact",
    return = "W"
  )
  g_manual <- (d$tau - as.numeric(d$Y <= as.vector(d$D %*% d$beta) + d$q)) * (d$D - d$mD)
  moment_manual <- crossprod(g_manual, wt$W)

  out1 <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau, X0 = d$X0,
    distance = nm, k = 3L, s = 6L, method = "exact",
    se = FALSE, return = "W"
  )
  expect_equal(unname(out1$moment), unname(moment_manual), tolerance = 1e-12)

  DX <- distance_eval_batch(spec, d$X, d$X0)
  ORD <- apply(DX, 2L, function(col) order(col, seq_along(col)))
  RKS <- apply(ORD, 2L, ord_to_ranks)

  rw_power <- function(alpha) {
    as_rank_weight(function(r, ctx = NULL) {
      w <- 1 / (as.numeric(r)^alpha)
      w / sum(w)
    }, name = paste0("inv_rank^", alpha))
  }
  rw1 <- rw_power(1)
  rw2 <- rw_power(2)

  W1 <- rank_weight_eval_batch(rw1, RKS)
  W2 <- rank_weight_eval_batch(rw2, RKS)
  moment_rw1 <- crossprod(g_manual, W1)
  moment_rw2 <- crossprod(g_manual, W2)

  expect_true(all(abs(colSums(W1) - 1) < 1e-12))
  expect_true(all(abs(colSums(W2) - 1) < 1e-12))
  expect_false(isTRUE(all.equal(unname(moment_rw1), unname(moment_rw2))))
})
