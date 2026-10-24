test_that("duplicate raw covariates share block mass and permutation-invariant inference", {
  X <- matrix(c(0, 0, 0, 1, 2, 2, 3, 4), ncol = 1L)
  y <- c(1, 7, -3, 4, 10, -2, 5, 11)
  targets <- matrix(c(0, 1.5), ncol = 1L)
  perm <- c(3, 1, 8, 6, 4, 2, 7, 5)
  fit <- suppressWarnings(np_dnn(X = X, Y = y, X0 = targets, k = 1, s = 5,
                                 return = c("W", "ranks")))
  reordered <- suppressWarnings(np_dnn(X = X[perm, , drop = FALSE], Y = y[perm],
                                       X0 = targets, k = 1, s = 5,
                                       return = "W"))
  expect_equal(fit$W[1, ], fit$W[2, ])
  expect_equal(fit$W[2, ], fit$W[3, ])
  expect_equal(fit$mu, reordered$mu, tolerance = 1e-13)
  expect_equal(fit$se, reordered$se, tolerance = 1e-13)
  expect_equal(suppressWarnings(vcov(fit)), suppressWarnings(vcov(reordered)),
               tolerance = 1e-13)
})

test_that("uniform-weight advice distinguishes duplicate blocks from a small order", {
  messages <- character()
  fit <- withCallingHandlers(np_dnn(X = matrix(0, 8L, 1L), x0 = 0,
    Y = seq_len(8L), s = 5L), warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_length(messages, 2L)
  advice <- messages[grepl("near-", messages, fixed = TRUE)]
  expect_match(messages[1L], "averaged within distance blocks", fixed = TRUE)
  expect_match(advice, "A larger `s` cannot distinguish observations at equal computed distance", fixed = TRUE)
  expect_equal(fit$W, matrix(1 / 8, 8L, 1L), tolerance = 1e-15)
})

test_that("block means preserve represented subnormals and recover finite overflow means", {
  eta <- .Machine$double.xmin * .Machine$double.eps
  X <- matrix(c(0, 1, 1, 1), ncol = 1L)
  ranks <- matrix(1:4, ncol = 1L)
  jk <- RdistNN:::.dnn_loo_centered(c(0, eta, eta, eta), ranks, rep(1, 3L),
    rep(1, 4L), groups = c(1L, 2L, 2L, 2L))
  expect_equal(as.numeric(jk$loo) / eta, c(3, 2, 2, 2), tolerance = 0)
  rw <- as_rank_weight(function(r) rep(1e308, length(r)))
  W <- suppressWarnings(weights_matrix(matrix(0, 3L, 1L), x0 = 0,
                                       ranking = "l2", weight = rw))$W
  expect_true(all(is.finite(W)))
  expect_equal(as.numeric(W) / 1e308, rep(1, 3L), tolerance = 1e-15)
})

test_that("scalar paired-ratio tuning reconstruction derives groups from retained raw X", {
  X <- matrix(c(0, 0, 0, 1, 2, 2, 3, 4), ncol = 1L)
  y <- c(1, 7, -3, 4, 10, -2, 5, 11)
  tuning <- new_tuning("duplicate_scalar_sensitivity", eval = function(...) list(s = 3L),
    eval_grid = function(...) data.frame(s = c(2L, 3L, 5L)))
  f <- suppressWarnings(affine_ratio_dnn(X = X, x0 = 0, numerator = y,
    denominator = rep(1, 8L), s = 3L))
  payload <- RdistNN:::.plot_tuning_payload(f, tuning, TRUE)
  grid <- RdistNN:::.affine_ratio_tuning_grid(f, f$.model$signal, f$.model$ranks,
    f$.model, TRUE, "fast", f$.model$denominator_tol, payload, reuse_default = FALSE)
  for (i in seq_len(nrow(grid$candidates))) {
    s <- grid$candidates$s[i]
    W <- .block_oracle_weights(X, matrix(0), .block_exact_table(8L, 1L, s))
    rows <- grid$estimates$candidate_id == grid$candidates$candidate_id[i]
    expect_equal(grid$estimates$estimate[rows], as.numeric(crossprod(y, W)), tolerance = 1e-13)
  }
})

test_that("duplicate-block estimates and uncertainty agree across all fitted metrics and engines", {
  X <- cbind(c(0, 0, 0, 1, 2, 2, 3, 4), c(0, 0, 0, 1, 0, 0, 2, 3))
  y <- c(1, 7, -3, 4, 10, -2, 5, 11)
  X0 <- rbind(c(0, 0), c(1.5, 1), c(4, 3))
  perm <- c(3L, 1L, 8L, 6L, 4L, 2L, 7L, 5L)
  metrics <- list("l1", "l2", distance_l1_scaled("custom", scale_vec = c(2, 3)),
                  distance_l2_scaled("custom", scale_vec = c(2, 3)), distance_mahalanobis())
  original <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(original), add = TRUE)
  for (metric in metrics) for (method in c("exact", "asymp")) {
    native <- NULL
    for (compiled in c(TRUE, FALSE)) {
      options(RdistNN.use_compiled = compiled)
      fit <- suppressWarnings(np_dnn(X = X, Y = y, X0 = X0, distance = metric,
        k = 2L, s = 5L, method = method, normalize = FALSE, return = c("W", "ranks")))
      moved <- suppressWarnings(np_dnn(X = X[perm, , drop = FALSE], Y = y[perm],
        X0 = X0[3:1, , drop = FALSE], distance = metric,
        k = 2L, s = 5L, method = method, normalize = FALSE, return = "W"))
      expect_equal(fit$W, moved$W[order(perm), 3:1], tolerance = 2e-13)
      expect_equal(fit$mu, moved$mu[3:1], tolerance = 2e-12)
      expect_equal(fit$se, moved$se[3:1], tolerance = 2e-12)
      if (compiled) native <- fit
      else {
        expect_equal(fit$W, native$W, tolerance = 2e-13)
        expect_equal(fit$mu, native$mu, tolerance = 2e-12)
        expect_equal(fit$se, native$se, tolerance = 2e-12)
      }
    }
  }
})

test_that("supplied PLQR coefficients use independently averaged deletion meat and density bread", {
  X <- matrix(c(0, 0, 0, 1, 2, 2, 3, 4), ncol = 1L)
  X0 <- matrix(c(0, 2), ncol = 1L)
  y <- c(1, 7, -3, 4, 10, -2, 5, 11)
  D <- cbind(seq_len(8L), c(1, 0, -2, 1, 3, 2, -1, 1))
  mD <- matrix(.1, 8L, 2L)
  g <- rep(.2, 8L)
  theta <- c(.3, -.1)
  density <- seq(.3, 1, length.out = 8L)
  score <- (D - mD) * (.5 - as.numeric(y <= as.numeric(D %*% theta) + g))
  W <- .block_oracle_weights(X, X0, .block_exact_table(8L, 1L, 5L))
  loo <- .block_oracle_loo(X, X0, score, .block_exact_table(7L, 1L, 5L))
  for (meat in c("analytic", "jackknife")) {
    fit <- suppressWarnings(plqr_dnn(X = X, X0 = X0, Y = y, D = D,
      mD = mD, g = g, theta_eval = theta, tau = .5, s = 5L,
      density = density, theta_meat = meat, return = "W"))
    expect_equal(unname(fit$moment), crossprod(score, W), tolerance = 1e-13)
    for (j in 1:2) {
      centered <- sweep(loo[, , j], 2L, colMeans(loo[, , j]), "-")
      moment_cov <- crossprod(centered) * 7 / 8
      bread <- crossprod(D - mD, D * (density * W[, j]))
      inverse <- solve(bread)
      actual_meat <- if (meat == "analytic") {
        crossprod(D - mD, (D - mD) * W[, j]) * .25 * sum(W[, j]^2)
      } else moment_cov
      covariance <- inverse %*% actual_meat %*% t(inverse)
      expect_equal(unname(fit$theta_vcov[[j]]), covariance, tolerance = 1e-12)
      expect_equal(unname(fit$jacobian[[j]]), bread, tolerance = 1e-13)
      if (meat == "jackknife") {
        expect_equal(unname(fit$moment_vcov[[j]]), moment_cov, tolerance = 1e-12)
        expect_equal(unname(fit$moment_se[, j]^2), diag(moment_cov), tolerance = 1e-12)
      }
    }
  }
})

test_that("complete support blocks are centered before averaging and unused extremes stay unused", {
  X <- matrix(c(rep(0, 5L), 1, 2, 3), ncol = 1L)
  changes <- c(0, 4, 8, -4, 12)
  signal <- c(2^54 + changes, .Machine$double.xmax, -.Machine$double.xmax, 0)
  ranks <- matrix(1:8, ncol = 1L)
  groups <- c(rep(1L, 5L), 2L, 3L, 4L)
  for (compiled in c(TRUE, FALSE)) {
    local_options <- options(RdistNN.use_compiled = compiled)
    jk <- RdistNN:::.dnn_loo_centered(signal, ranks, .block_exact_table(7L, 1L, 5L),
      .block_exact_table(8L, 1L, 5L), groups = groups)
    reference <- c(vapply(1:5, function(i) mean(changes[-i]), numeric(1L)),
                   rep(mean(changes), 3L))
    expect_equal(as.numeric(jk$loo), reference, tolerance = 1e-13)
    options(local_options)
  }
})

test_that("bare supplied rank APIs stay strict and explicit distance metadata defines averaging", {
  X <- matrix(c(-1, 1, -1, 2, 2, 3), ncol = 1L)
  X0 <- matrix(0)
  ranks <- matrix(c(6L, 2L, 1L, 4L, 3L, 5L), ncol = 1L)
  rw <- rank_weight_dnn_exact(1L, 3L)
  strict <- weights_matrix(X, X0, ranking = ranks, weight = rw)$W
  expect_equal(as.numeric(strict), .block_exact_table(6L, 1L, 3L)[ranks[, 1L]])
  dx <- matrix(abs(X[, 1L]), ncol = 1L)
  W <- suppressWarnings(weights_matrix(X, X0, ranking = list(ranks = ranks, dx = dx), weight = rw))$W
  expect_equal(W, .block_oracle_weights(X, X0, .block_exact_table(6L, 1L, 3L), ranks), tolerance = 1e-14)
  signal <- c(1, 7, -3, 4, 10, -2)
  jk <- RdistNN:::.dnn_loo_centered(signal, ranks, .block_exact_table(5L, 1L, 3L),
    .block_exact_table(6L, 1L, 3L), groups = RdistNN:::.dnn_distance_groups(dx))
  ref <- .block_oracle_loo(X, X0, signal, .block_exact_table(5L, 1L, 3L), ranks)
  expect_equal(as.numeric(RdistNN:::.dnn_restore_loo(jk$loo, jk$offset)), as.numeric(ref), tolerance = 1e-13)
  d <- suppressWarnings(weights_dnn(matrix(c(-1, 1, -1)), x0 = 0, k = 1L, s = 2L))
  expect_equal(as.numeric(d$W), rep(1 / 3, 3L), tolerance = 1e-15)
})

test_that("GMM localization and scalar pure-score subset aggregation share the X-aware rule", {
  X <- matrix(c(0, 0, 0, 1, 2, 2, 3, 4), ncol = 1L)
  y <- c(1, 7, -3, 4, 10, -2, 5, 11)
  perm <- c(3L, 1L, 8L, 6L, 4L, 2L, 7L, 5L)
  data <- data.frame(X = X[, 1L], Y1 = y, Y2 = y^2)
  model <- .test_mean2_gmm_model()
  f <- suppressWarnings(gmm_dnn(data = data, x_cols = "X", x0 = 0,
    model = model, k = 1L, s = 5L, return = "W"))
  other <- suppressWarnings(gmm_dnn(data = data[perm, ], x_cols = "X", x0 = 0,
    model = model, k = 1L, s = 5L, return = "W"))
  expect_equal(f$W, .block_oracle_weights(X, matrix(0), .block_exact_table(8L, 1L, 5L)))
  expect_equal(f$theta, other$theta, tolerance = 1e-9)
  expect_equal(f$se, other$se, tolerance = 1e-9)
  expect_equal(suppressWarnings(vcov(f)), suppressWarnings(vcov(other)), tolerance = 1e-9)
  sc <- score_np_y(x0 = 0, x_cols = "X", y_col = "Y1")
  ctx <- score_setup(sc, data)
  scores <- score_eval_batch(sc, matrix(1:8, ncol = 1L), data, ctx)
  expect_equal(as.numeric(crossprod(scores, f$W)),
               as.numeric(crossprod(y, f$W)))
})

test_that("all estimators, saved inference and candidate reconstructions use averaged blocks", {
  X <- cbind(c(0, 0, 0, 1, 2, 2, 3, 4), c(0, 0, 0, 1, 0, 0, 2, 3))
  X0 <- rbind(c(0, 0), c(2, 0))
  y <- c(1, 7, -3, 4, 10, -2, 5, 11)
  den <- c(1, 3, .5, 2, 4, 2, 1, 3)
  D <- cbind(seq_len(8L), c(1, 0, -2, 1, 3, 2, -1, 1))
  perm <- c(3L, 1L, 8L, 6L, 4L, 2L, 7L, 5L)
  tuning <- new_tuning("duplicate_reconstruction", eval = function(...) list(s = 3L),
    eval_grid = function(...) data.frame(s = c(2L, 3L, 5L)))
  make <- function(kind, ids, grid = FALSE, engine = "fast", meat = "analytic") {
    args <- list(X = X[ids, , drop = FALSE], X0 = X0,
                 return = c("W", "ranks", "dx", "ess"), jackknife = engine)
    if (grid) args <- c(args, list(tuning = tuning, tuning_mode = "grid"))
    else if (kind != "tdnn") args$s <- 5L
    f <- switch(kind,
      np = np_dnn, tdnn = np_tdnn, cate = cate_dnn,
      ratio = affine_ratio_dnn, clate = clate_dnn, plqr = plqr_dnn)
    args <- c(args, switch(kind,
      np = list(Y = y[ids]),
      tdnn = list(Y = y[ids], s1 = 2L, s2 = 5L, dim = 2L),
      cate = list(pseudo_outcome = y[ids]),
      ratio = list(numerator = y[ids], denominator = den[ids]),
      clate = list(Y = y[ids], D = rep(1, length(ids)), A = rep(1, length(ids)),
        e = rep(.5, length(ids)), m_y0 = rep(0, length(ids)), m_y1 = y[ids],
        m_d0 = rep(0, length(ids)), m_d1 = rep(1, length(ids))),
      plqr = list(Y = y[ids], D = D[ids, , drop = FALSE],
        g = rep(0, length(ids)), mD = matrix(0, length(ids), 2L),
        theta_eval = c(.3, -.1), tau = .5,
        density = if (grid) NULL else rep(.5, length(ids)), theta_meat = meat)))
    suppressWarnings(do.call(f, args))
  }
  for (kind in c("np", "tdnn", "cate", "ratio", "clate", "plqr")) {
    for (jk in c("fast", "loo")) {
      f <- make(kind, 1:8, engine = jk)
      other <- make(kind, perm, engine = jk)
      expect_equal(f$W, other$W[order(perm), ], tolerance = 1e-14)
      estimate <- if (kind == "plqr") f$moment else suppressWarnings(coef(f))
      reordered <- if (kind == "plqr") other$moment else suppressWarnings(coef(other))
      expect_equal(estimate, reordered, tolerance = 1e-12)
      expect_equal(f$se, other$se, tolerance = 1e-12)
      if (kind != "plqr") {
        expect_equal(suppressWarnings(vcov(f)), suppressWarnings(vcov(other)), tolerance = 1e-12)
      } else {
        expect_equal(f$moment_se, other$moment_se, tolerance = 1e-12)
        expect_equal(f$theta_vcov, other$theta_vcov, tolerance = 1e-12)
      }
      restored <- unserialize(serialize(f, NULL))
      expect_equal(suppressWarnings(predict(restored, X0 = X0[2:1, , drop = FALSE], se.fit = TRUE)),
        suppressWarnings(predict(other, X0 = X0[2:1, , drop = FALSE], se.fit = TRUE)), tolerance = 1e-12)
      if (kind != "plqr") {
        actual <- suppressWarnings(diagnostics(f))
        expect_equal(actual$n_eff, 1 / colSums(f$W^2), tolerance = 1e-13)
      }
    }
    if (kind == "tdnn") next
    fit <- make(kind, 1:8, grid = TRUE)
    for (j in 1:3) {
      view <- suppressWarnings(candidate_fit(fit, j))
      s <- view$meta$s
      expected_W <- .block_oracle_weights(X, X0, .block_exact_table(8, 1L, s))
      expect_equal(view$W, expected_W, tolerance = 1e-14)
      expect_equal(view$ess, 1 / colSums(expected_W^2), tolerance = 1e-13)
      if (kind != "plqr") expect_equal(suppressWarnings(diagnostics(view))$n_eff,
        view$ess, tolerance = 1e-13)
      profile <- suppressWarnings(plot_data_weights(view))
      expect_equal(profile$weight[profile$target_id == 1L], expected_W[order(view$ranks[, 1L]), 1L])
    }
    profiles <- RdistNN:::.plot_candidate_weight_frame(fit, "all", "test")$data
    expect_true("target_id" %in% names(profiles))
    for (j in 1:3) for (target in 1:2) {
      row <- profiles$candidate_id == fit$tuning_grid$candidates$candidate_id[j] &
        profiles$target_id == target
      view <- suppressWarnings(candidate_fit(fit, j))
      expect_equal(profiles$weight[row], view$W[order(view$ranks[, target]), target])
    }
    for (j in 1:3) {
      view <- suppressWarnings(candidate_fit(fit, j))
      actual <- fit$tuning_grid$diagnostics[j, ]
      expect_equal(actual$ess, min(1 / colSums(view$W^2)), tolerance = 1e-13)
      expect_equal(actual$max_weight, max(abs(view$W)), tolerance = 1e-13)
    }
  }
  plqr <- make("plqr", 1:8, meat = "jackknife")
  shifted <- make("plqr", perm, meat = "jackknife")
  expect_equal(plqr$moment_vcov, shifted$moment_vcov, tolerance = 1e-12)
  expect_equal(plqr$theta_vcov, shifted$theta_vcov, tolerance = 1e-12)
})

test_that("both actual original-signal block engines average each deletion", {
  X <- cbind(c(0, 0, 0, 1, 2, 2, 3, 4), 0)
  X0 <- rbind(c(0, 0), c(1.5, 0), c(4, 0))
  signal <- cbind(c(1, 7, -3, 4, 10, -2, 5, 11), c(-2, 3, 9, 1, 5, -7, 4, 8))
  groups <- RdistNN:::.dnn_distance_groups(distance_eval_batch(get_distance("l2"), X, X0))
  ranks <- suppressWarnings(weights_dnn(X, X0, k = 1, s = 5, return = "ranks"))$ranks
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  native_calls <- fallback_calls <- 0L
  native <- RdistNN:::.rdistnn_compiled_fun("rustat_block_rank_loo")
  expect_true(is.function(native))
  fallback <- RdistNN:::.dnn_block_loo_R
  testthat::local_mocked_bindings(
    rustat_block_rank_loo = function(...) {
      native_calls <<- native_calls + 1L
      native(...)
    },
    .dnn_block_loo_R = function(...) {
      fallback_calls <<- fallback_calls + 1L
      fallback(...)
    }, .package = "RdistNN")
  rules <- list(
    function(n) .block_exact_table(n, 1L, 5L),
    function(n) .block_exact_table(n, 2L, 5L),
    function(n) rank_weight_eval(rank_weight_dnn_asymp(2L, 5L, normalize = FALSE), seq_len(n)),
    function(n) tdnn_rank_table_exact(n, 2L, 5L, 2L),
    function(n) tdnn_rank_table_asymp(n, 2L, 5L, 2L, normalize = FALSE))
  for (rule in rules) {
    ref <- .block_oracle_loo(X, X0, signal, rule(nrow(X) - 1L))
    for (compiled in c(TRUE, FALSE)) {
      options(RdistNN.use_compiled = compiled)
      counts <- c(native_calls, fallback_calls)
      jk <- RdistNN:::.dnn_loo_centered(signal, ranks, rule(nrow(X) - 1L),
        tab_full = rule(nrow(X)), groups = groups)
      actual <- RdistNN:::.dnn_restore_loo(jk$loo, jk$offset)
      expect_equal(unname(actual), unname(ref), tolerance = 2e-12)
      expect_equal(c(native_calls, fallback_calls) - counts,
                   if (compiled) c(nrow(X0), 0L) else c(0L, nrow(X0)))
      covariance <- RdistNN:::.dnn_jackknife_vcov_from_loo(jk$loo, 2L, nrow(X), nrow(X0))$vcov
      for (j in seq_len(nrow(X0))) {
        centered <- sweep(ref[, , j], 2L, colMeans(ref[, , j]), "-")
        expect_equal(unname(covariance[[j]]), crossprod(centered) * 7 / 8,
                     tolerance = 2e-12)
      }
    }
  }
})

test_that("all-identical raw rows give sample and usual deletion means", {
  X <- matrix(0, 8L, 2L)
  y <- c(1, 8, -3, 4, 7, 2, -5, 6)
  targets <- rbind(c(0, 0), c(1, 1))
  for (k in c(1L, 2L)) for (jk in c("fast", "loo")) {
    f <- suppressWarnings(np_dnn(X = X, Y = y, X0 = targets, k = k, s = 5,
                                  jackknife = jk))
    expect_equal(f$W, matrix(1 / 8, 8L, 2L), tolerance = 1e-15)
    expect_equal(f$mu, rep(mean(y), 2L), tolerance = 1e-13)
    expect_equal(f$se, rep(sd(y) / sqrt(8), 2L), tolerance = 1e-13)
    expect_equal(unname(suppressWarnings(vcov(f))),
                 matrix(var(y) / 8, 2L, 2L), tolerance = 1e-13)
    survivor <- suppressWarnings(weights_dnn(X[-1L, , drop = FALSE], X0 = targets,
                                              k = k, s = 5L))$W
    expect_equal(survivor, matrix(1 / 7, 7L, 2L), tolerance = 1e-15)
    expect_gt(survivor[1L, 1L], f$W[2L, 1L])
  }
})

test_that("warning scopes avoid deletion and candidate floods and reject old tied fits", {
  X <- matrix(c(0, 0, 0, 1, 2, 2, 3, 4), ncol = 1L)
  y <- c(1, 7, -3, 4, 10, -2, 5, 11)
  tuning <- new_tuning("duplicate_grid", eval = function(...) list(s = 3L),
    eval_grid = function(...) data.frame(s = c(2L, 3L, 5L)))
  seen <- character()
  fit <- withCallingHandlers(np_dnn(X = X, Y = y, x0 = 0,
    jackknife = "loo", tuning = tuning, tuning_mode = "grid"), warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_length(seen, 1L)
  expect_match(seen, "continuous distance distribution", fixed = TRUE)
  seen <- character()
  views <- withCallingHandlers(candidate_fits(fit), warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_length(seen, 1L)
  expect_length(views, 3L)
  old <- fit
  old$.model$tie_policy <- NULL
  for (consume in list(function(x) predict(x, x0 = 0),
      function(x) candidate_fit(x, 1L), vcov, diagnostics, plot_data_weights)) {
    err <- tryCatch(consume(old), error = identity)
    expect_s3_class(err, "rdistnn_state_error")
    expect_identical(err$code, "state_invalid")
    expect_identical(err$argument, "object$.model$tie_policy")
    expect_match(conditionMessage(err), "Refit", fixed = TRUE)
  }
})
