test_that("covariate-aware ordering uses exact raw lexicographic tie keys", {
  .with_duplicate_policy_warnings({
  X <- rbind(c(1, 0), c(-1, 0), c(0, 1), c(0, -1), c(1, 0), c(2, 0))
  expected <- c(2L, 4L, 3L, 1L, 5L, 6L)
  for (distance in c("l1", "l2")) {
    out <- weights_dnn(X, x0 = c(0, 0), k = 1L, s = 2L, distance = distance,
      return = c("W", "ord", "ranks"))
    expect_identical(as.integer(out$ord), expected)
    expect_identical(as.integer(out$ranks), c(4L, 1L, 3L, 2L, 5L, 6L))
  }
  # A represented nonzero distance difference always precedes covariate keys.
  out <- weights_matrix(matrix(c(1, -1)), x0 = 0,
    ranking = list(dx = matrix(c(1, 1 + 2^-52))),
    weight = rank_weight_dnn_exact(k = 1L, s = 2L), return = "ord")
  expect_identical(as.integer(out$ord), 1:2)
  })
})

test_that("custom distances use raw covariates and supplied ranks remain explicit", {
  d <- .tie_fixture()
  zero <- new_distance("canonical_zero", eval_batch = function(X, X0, ctx) {
    matrix(0, nrow(X), nrow(X0))
  })
  expected <- c(2L, 4L, 3L, 6L, 1L, 5L)
  rw <- rank_weight_dnn_exact(k = 1L, s = 2L)
  for (ranking in list(zero, list(dx = matrix(0, 6L, 3L)),
      function(X, X0, ctx) list(dx = matrix(0, nrow(X), nrow(X0))))) {
    out <- weights_matrix(d$X, X0 = d$X0, ranking = ranking, weight = rw,
      return = c("ord", "ranks"))
    expect_identical(out$ord, matrix(expected, 6L, 3L))
  }
  explicit <- matrix(6:1, 6L, 3L)
  for (field in c("ord", "ranks")) {
    out <- weights_matrix(d$X, X0 = d$X0,
      ranking = stats::setNames(list(explicit), field), weight = rw,
      return = field)
    expect_identical(out[[field]], explicit)
  }
})

test_that("native and R covariate order helpers share canonical full and subset ties", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  d <- .tie_fixture()
  cases <- list(
    list(order_by_l1, list()), list(order_by_l2, list()),
    list(order_by_l1_scaled, list(rule = "custom", scale_vec = c(1, 1))),
    list(order_by_l2_scaled, list(rule = "custom", scale_vec = c(1, 1))))
  symbols <- c("rustat_order_by_l1", "rustat_order_by_l2",
    "rustat_order_by_l1_scaled", "rustat_order_by_l2_scaled")
  for (native in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = native)
    for (i in seq_along(cases)) {
      if (native) expect_true(is.function(.rdistnn_compiled_fun(symbols[[i]])))
      for (ids in list(1:6, d$perm, c(6L, 2L, 5L, 1L))) {
        local <- d$X[ids, , drop = FALSE]
        result <- do.call(cases[[i]][[1L]], c(list(local, c(0, 0)), cases[[i]][[2L]]))
        expected <- .tie_reference_order(result$dx, local, ids)
        expect_identical(ids[result$ord], ids[expected])
      }
    }
  }
})

test_that("estimator fits and saved predictions preserve distinct-row tie permutations", {
  d <- .tie_fixture()
  W <- .tie_reference_weights(d$X, d$X0, 1L, 3L)
  expected <- as.numeric(crossprod(d$y, W))
  for (kind in c("np", "tdnn", "cate", "clate")) {
    for (jk in c("fast", "loo")) {
      args <- if (kind == "tdnn") list() else list(s = 3L)
      fit <- do.call(.tie_fit, c(list(kind, d, jackknife = jk), args))
      moved <- do.call(.tie_fit, c(list(kind, d, at = d$perm, jackknife = jk), args))
      # Rank/support identities are exact; only signal accumulation gets slack.
      expect_identical(moved$ranks[order(d$perm), ], fit$ranks)
      expect_identical(moved$W[order(d$perm), ], fit$W)
      expect_identical(colSums(moved$W != 0), colSums(fit$W != 0))
      expect_equal(moved$ess, fit$ess, tolerance = 1e-14)
      expect_equal(coef(moved), coef(fit), tolerance = 1e-12)
      expect_equal(moved$se, fit$se, tolerance = 1e-12)
      expect_equal(vcov(moved), vcov(fit), tolerance = 1e-12)
      if (kind != "tdnn") expect_equal(unname(coef(fit)), expected, tolerance = 1e-12)
      restored <- unserialize(serialize(fit, NULL))
      targets <- d$X0[c(2L, 1L, 2L), , drop = FALSE]
      expect_equal(predict(restored, X0 = targets, se.fit = TRUE),
        predict(moved, X0 = targets, se.fit = TRUE), tolerance = 1e-12)
    }
  }
})

test_that("candidate views retain canonical full-sample order under permutations", {
  d <- .tie_fixture()
  tuning <- new_tuning("canonical_tie_grid", eval = function(...) list(s = 3L),
    eval_grid = function(...) data.frame(s = c(2L, 3L, 4L)))
  for (kind in c("np", "cate", "clate")) {
    fit <- .tie_fit(kind, d, tuning = tuning, tuning_mode = "grid")
    moved <- .tie_fit(kind, d, at = d$perm, tuning = tuning, tuning_mode = "grid")
    for (i in 1:3) {
      view <- candidate_fit(fit, i)
      other <- candidate_fit(moved, i)
      expect_identical(view$W, other$W[order(d$perm), ])
      expect_equal(coef(view), coef(other), tolerance = 1e-12)
      expect_equal(view$se, other$se, tolerance = 1e-12)
      expect_equal(predict(view, X0 = d$X0), predict(other, X0 = d$X0), tolerance = 1e-12)
    }
  }
})

test_that("PLQR moments and saved inference use covariate-only tie order", {
  d <- .tie_fixture()
  D <- seq_len(nrow(d$X))
  make <- function(ids, density = NULL) {
    plqr_dnn(X = d$X[ids, , drop = FALSE], X0 = d$X0,
      D = matrix(D[ids]), mD = matrix(0, length(ids)),
      g = rep(0, length(ids)), Y = d$y[ids], theta_eval = .3, tau = .5,
      s = 3L, density = density, return = c("W", "ranks", "ord"))
  }
  expected_W <- .tie_reference_weights(d$X, d$X0, 1L, 3L)
  score <- (.5 - as.numeric(d$y <= .3 * D)) * D
  for (density in list(NULL, rep(.5, nrow(d$X)))) {
    fit <- make(seq_len(nrow(d$X)), density)
    other <- make(d$perm, density)
    expect_equal(fit$W, expected_W, tolerance = 1e-14)
    expect_identical(fit$ranks, other$ranks[order(d$perm), ])
    expect_equal(as.numeric(fit$moment), as.numeric(crossprod(score, expected_W)),
      tolerance = 1e-12)
    expect_equal(fit$moment, other$moment, tolerance = 1e-12)
    expect_equal(fit$moment_se, other$moment_se, tolerance = 1e-12)
    restored <- unserialize(serialize(fit, NULL))
    target <- d$X0[c(2L, 1L, 2L), , drop = FALSE]
    expect_equal(predict(restored, X0 = target, se.fit = !is.null(density)),
      predict(other, X0 = target, se.fit = !is.null(density)), tolerance = 1e-12)
  }
})

test_that("identical covariate rows retain ranks but share weight mass", {
  .with_duplicate_policy_warnings({
  X <- matrix(c(1, -1, 1, -1))
  out <- weights_dnn(X, x0 = 0, k = 1L, s = 2L, return = "ord")
  expect_identical(as.integer(out$ord), c(2L, 4L, 1L, 3L))
  # Swapping signals on identical rows preserves their averaged block mass.
  first <- cate_dnn(X = X, x0 = 0, s = 2L, pseudo_outcome = c(1, 2, 3, 4), se = FALSE)
  second <- cate_dnn(X = X, x0 = 0, s = 2L, pseudo_outcome = c(1, 4, 3, 2), se = FALSE)
  expect_equal(first$tau, second$tau)
  })
})

test_that("weight plotting respects retained canonical orders", {
  d <- .tie_fixture()
  for (kind in c("np", "tdnn", "cate", "clate")) {
    args <- if (kind == "tdnn") list() else list(s = 3L)
    fit <- do.call(.tie_fit, c(list(kind, d), args))
    view <- plot_data_weights(fit)
    expect_identical(view$obs_id[view$target_id == 1L], as.integer(fit$ord[, 1L]))
    fit$ord <- fit$ranks <- fit$.model$ranks <- NULL
    if (!is.null(fit$.model$geometry)) fit$.model$geometry$ranks <- NULL
    expect_identical(plot_data_weights(fit)$obs_id, view$obs_id)
  }
})

test_that("raw covariates decide ties even in null directions of fitted geometry", {
  X <- cbind(c(2, -2, 1, -1), rep(2^600, 4L))
  zero_projection <- list(version = 2L, scale = 1, origin = c(0, 0),
    factor = matrix(c(0, 0, 0, 1), 2L), rank = 1L)
  fit <- weights_dnn(X, x0 = c(0, 2^600), k = 1L, s = 2L,
    distance = distance_mahalanobis(), ctx = zero_projection, return = "ord")
  expect_identical(as.integer(fit$ord), c(2L, 4L, 3L, 1L))
})

test_that("scaled ties retain raw keys when transformed coordinates underflow", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  X <- matrix(c(2, -2, 1, -1) * .Machine$double.xmin)
  expected <- c(2L, 4L, 3L, 1L)
  for (native in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = native)
    for (ordering in list(order_by_l1_scaled, order_by_l2_scaled)) {
      result <- ordering(X, 0, rule = "custom", scale_vec = .Machine$double.xmax)
      expect_identical(result$dx, rep(0, 4L))
      expect_identical(result$ord, expected)
    }
    for (distance in list(distance_l1_scaled("custom", scale_vec = .Machine$double.xmax),
        distance_l2_scaled("custom", scale_vec = .Machine$double.xmax))) {
      fit <- weights_dnn(X, x0 = 0, k = 1L, s = 2L, distance = distance, return = "ord")
      expect_identical(as.integer(fit$ord), expected)
    }
  }
})
