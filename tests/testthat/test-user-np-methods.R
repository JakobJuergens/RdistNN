# Independent small-sample exact DNN weights by enumerating all subsamples.
.np_ref_weights <- function(X, X0, s, k = 1L) {
  n <- nrow(X)
  subsets <- combn(n, s)
  W <- matrix(0, n, nrow(X0))
  for (j in seq_len(nrow(X0))) {
    d <- rowSums(sweep(X, 2, X0[j, ], "-")^2)
    for (a in seq_len(ncol(subsets))) {
      ids <- subsets[, a]
      take <- ids[.tie_reference_order(d[ids], X[ids, , drop = FALSE], ids)][seq_len(k)]
      W[take, j] <- W[take, j] + 1 / (k * ncol(subsets))
    }
  }
  .block_oracle_average_weights(W, X)
}

.np_ref_cov <- function(X, X0, y, s, k = 1L, s2 = NULL, dim = ncol(X)) {
  n <- nrow(X)
  loo <- t(vapply(seq_len(n), function(i) {
    W <- .np_ref_weights(X[-i, , drop = FALSE], X0, s, k)
    if (!is.null(s2)) {
      b <- (s2 / s)^(2 / dim)
      W <- (b * .np_ref_weights(X[-i, , drop = FALSE], X0, s2) - W) / (b - 1)
    }
    as.numeric(crossprod(y[-i], W))
  }, numeric(nrow(X0))))
  centered <- sweep(loo, 2, colMeans(loo), "-")
  (n - 1) / n * crossprod(centered)
}

.np_method_data <- function() list(
  X = cbind(a = c(-2, -1, 0, 0, 1, 2, 4), b = c(3, 1, 0, 0, 3, 2, -1)),
  X0 = rbind(c(-.2, .1), c(1.7, 1), c(0, 0)),
  y = c(2, -1, 3, 7, -2, 4, 1)
)

test_that("NP methods expose means and joint common-deletion covariance", {
  .with_duplicate_policy_warnings({
  d <- .np_method_data()
  for (jk in c("fast", "loo")) {
    f <- np_dnn(X = d$X, X0 = d$X0, Y = d$y, k = 2, s = 3, jackknife = jk, return = character())
    ref <- as.numeric(crossprod(d$y, .np_ref_weights(d$X, d$X0, 3, 2)))
    V <- .np_ref_cov(d$X, d$X0, d$y, 3, 2)
    expect_equal(unname(coef(f)), ref, tolerance = 1e-12)
    expect_equal(unname(vcov(f)), V, tolerance = 1e-12)
    expect_gt(abs(vcov(f)[1, 2]), .01)
    expect_equal(unname(diag(vcov(f))), f$se^2, tolerance = 1e-12)
    expect_named(coef(f), paste0("target_", 1:3))
    sm <- summary(f, level = .9)
    expect_s3_class(sm, "summary.np_result")
    expect_equal(sm$coefficients$estimate, ref)
    expect_equal(sm$coefficients$statistic, ref / sqrt(diag(V)))
    expect_equal(unname(confint(f, parm = 2, level = .9)),
      matrix(ref[2] + c(-1, 1) * qnorm(.95) * sqrt(V[2, 2]), 1), tolerance = 1e-12)
    expect_output(print(sm), "summary.np_result")
    expect_equal(predict(f, X0 = d$X0), ref)
    pr <- predict(f, newdata = d$X0, interval = "confidence", se.fit = TRUE, level = .9)
    expect_equal(pr$se.fit, sqrt(diag(V)))
    expect_equal(unname(pr$fit[, 2:3]), unname(confint(f, level = .9)))
  }
  })
})

test_that("NP prediction retains tabular inputs and fitted geometry through serialization", {
  .with_duplicate_policy_warnings({
  d <- .np_method_data()
  dat <- data.frame(a = d$X[, 1], b = d$X[, 2], y = d$y)
  metric <- distance_l2_scaled("sd")
  scale <- apply(d$X, 2, sd)
  refX <- sweep(d$X, 2, scale, "/")
  ref0 <- sweep(d$X0, 2, scale, "/")
  f <- np_dnn(data = dat, x_cols = c("a", "b"), y_col = "y", X0 = d$X0,
    distance = metric, s = 3, return = character())
  restored <- unserialize(serialize(f, NULL))
  nd <- data.frame(extra = 1:3, b = d$X0[, 2], a = d$X0[, 1])
  expect_equal(predict(restored, newdata = nd), f$mu)
  expect_equal(unname(vcov(restored)), .np_ref_cov(refX, ref0, d$y, 3), tolerance = 1e-12)
  expect_equal(predict(restored, x0 = d$X0[1, ]), f$mu[1])
  expect_error(predict(restored, newdata = nd["a"]), "b.*x_cols")
  expect_null(f$W)
  })
})

test_that("NP zero variance and absent uncertainty remain distinguishable", {
  .with_duplicate_policy_warnings({
  d <- .np_method_data()
  f <- np_dnn(X = d$X, X0 = d$X0, Y = rep(5, 7), s = 3)
  expect_equal(unname(vcov(f)), matrix(0, 3, 3))
  sm <- summary(f)$coefficients
  expect_equal(sm$std.error, rep(0, 3))
  expect_equal(sm$statistic, rep(Inf, 3))
  expect_equal(sm$p.value, rep(0, 3))
  expect_equal(unname(confint(f)), matrix(5, 3, 2))
  expect_equal(unname(predict(f, X0 = d$X0, interval = "confidence")[, 2:3]), matrix(5, 3, 2))
  zero <- np_dnn(X = d$X, X0 = d$X0, Y = rep(0, 7), s = 3)
  expect_true(all(is.na(summary(zero)$coefficients$statistic)))
  expect_true(all(is.na(summary(zero)$coefficients$p.value)))
  expect_equal(unname(confint(zero)), matrix(0, 3, 2))
  f0 <- np_dnn(X = d$X, X0 = d$X0, Y = d$y, s = 3, se = FALSE, return = character())
  expect_true(all(is.na(summary(f0)$coefficients$std.error)))
  expect_error(vcov(f0), "se = TRUE")
  expect_error(confint(f0), "se = TRUE")
  expect_equal(predict(f0, X0 = d$X0, se.fit = TRUE)$se.fit,
    sqrt(diag(.np_ref_cov(d$X, d$X0, d$y, 3))))
  })
})

test_that("NP candidate views update covariance and prediction controls", {
  .with_duplicate_policy_warnings({
  d <- .np_method_data()
  tun <- new_tuning("np_methods", eval = function(...) list(s = 3L),
    eval_grid = function(...) data.frame(s = c(2L, 3L, 4L)))
  f <- np_dnn(X = d$X, X0 = d$X0, Y = d$y, tuning = tun, tuning_mode = "grid", return = character())
  for (i in 1:3) {
    v <- candidate_fit(f, i)
    s <- f$tuning_grid$candidates$s[i]
    expect_equal(v$.model$s, s)
    expect_equal(predict(v, X0 = d$X0), v$mu)
    expect_equal(unname(vcov(v)), .np_ref_cov(d$X, d$X0, d$y, s), tolerance = 1e-12)
  }
  expect_equal(f$.model$s, 3L)
  })
})

test_that("NP methods validate public targets and inference controls", {
  .with_duplicate_policy_warnings({
  d <- .np_method_data()
  f <- np_dnn(X = d$X, X0 = d$X0, Y = d$y, s = 3)
  expect_error(summary(f, level = 0), "level")
  expect_error(confint(f, parm = 1.5), "integer target indices")
  expect_error(predict(f, X0 = matrix(1, 2, 1)), "ncol")
  expect_error(predict(f, X0 = d$X0, x0 = d$X0[1, ]), "either")
  expect_error(predict(f, X0 = d$X0, se.fit = NA), "se.fit")
  expect_error(predict(f, x0 = c(NA, 0)), "NA")
  f$.model <- NULL
  expect_error(predict(f, X0 = d$X0), "fitted model state")
  })
})

test_that("TDNN inherits NP methods with actual signed covariance", {
  .with_duplicate_policy_warnings({
  d <- .np_method_data()
  for (jk in c("fast", "loo")) {
    f <- np_tdnn(X = d$X, X0 = d$X0, Y = d$y, s1 = 2, s2 = 4, dim = 2,
      jackknife = jk, return = character())
    expect_equal(unname(vcov(f)), .np_ref_cov(d$X, d$X0, d$y, 2, s2 = 4, dim = 2), tolerance = 1e-12)
    expect_equal(predict(unserialize(serialize(f, NULL)), X0 = d$X0), f$mu)
    expect_equal(predict(f, X0 = d$X0, se.fit = TRUE)$se.fit, f$se)
    expect_equal(unname(coef(f)), f$mu)
  }
  })
})

test_that("loo covariance preserves deletion-specific fitted geometry", {
  .with_duplicate_policy_warnings({
  d <- .np_method_data()
  loo <- t(vapply(seq_len(nrow(d$X)), function(i) {
    Xi <- d$X[-i, , drop = FALSE]
    scales <- apply(Xi, 2, sd)
    W <- .np_ref_weights(sweep(Xi, 2, scales, "/"), sweep(d$X0, 2, scales, "/"), 3)
    as.numeric(crossprod(d$y[-i], W))
  }, numeric(nrow(d$X0))))
  C <- sweep(loo, 2, colMeans(loo), "-")
  V <- 6 / 7 * crossprod(C)
  f <- np_dnn(X = d$X, X0 = d$X0, Y = d$y, s = 3,
    distance = distance_l2_scaled("sd"), jackknife = "loo", return = character())
  expect_equal(unname(vcov(f)), V, tolerance = 1e-12)
  expect_equal(f$se^2, diag(V), tolerance = 1e-12)
  scales <- apply(d$X, 2, sd)
  frozen <- .np_ref_cov(sweep(d$X, 2, scales, "/"), sweep(d$X0, 2, scales, "/"), d$y, 3)
  expect_equal(predict(f, X0 = d$X0, se.fit = TRUE)$se.fit^2, diag(frozen), tolerance = 1e-12)
  })
})

test_that("NP resolves registered custom geometry once for saved predictions", {
  .with_duplicate_policy_warnings({
  d <- .np_method_data()
  metric <- new_distance("np_saved_custom",
    setup = function(X) list(scale = apply(X, 2, sd)),
    eval = function(X, x0, ctx) {
      rowSums(sweep(sweep(X, 2, x0, "-"), 2, ctx$scale, "/")^2)
    })
  old <- get_distance("l2")
  on.exit(register_distance("l2", old))
  register_distance("l2", metric)
  f <- np_dnn(X = d$X, X0 = d$X0, Y = d$y, s = 3, return = character())
  register_distance("l2", distance_l1())
  f <- unserialize(serialize(f, NULL))
  expect_equal(predict(f, X0 = d$X0), f$mu)
  expect_equal(predict(f, X0 = d$X0, se.fit = TRUE)$se.fit, f$se)
  expect_equal(unname(diag(vcov(f))), f$se^2, tolerance = 1e-12)
  # An explicit context is canonical training state and bypasses setup.
  metric$setup <- function(X) stop("setup must not run")
  fc <- np_dnn(X = d$X, X0 = d$X0, Y = d$y, s = 3,
    distance = metric, ctx = list(scale = c(1, 5)), return = character())
  expect_equal(predict(fc, X0 = d$X0), fc$mu)
  })
})

test_that("asymptotic NP and TDNN covariance matches explicit beta-grid deletions", {
  .with_duplicate_policy_warnings({
  d <- .np_method_data()
  weights <- function(X, X0, s, k) {
    n <- nrow(X)
    table <- rowMeans(vapply(seq_len(k), function(a) {
      diff(pbeta((0:n) / n, a, s - a + 1))
    }, numeric(n)))
    W <- vapply(seq_len(nrow(X0)), function(j) {
      distance <- rowSums(sweep(X, 2, X0[j, ], "-")^2)
      table[order(.tie_reference_order(distance, X))]
    }, numeric(n))
    .block_oracle_average_weights(W, X)
  }
  for (tdnn in c(FALSE, TRUE)) for (normalize in c(FALSE, TRUE)) {
    n <- nrow(d$X)
    loo <- t(vapply(seq_len(n), function(i) {
      X <- d$X[-i, , drop = FALSE]
      W <- if (tdnn) {
        2 * weights(X, d$X0, 4, 1) - weights(X, d$X0, 2, 1)
      } else weights(X, d$X0, 3, 2)
      as.numeric(crossprod(d$y[-i], W))
    }, numeric(nrow(d$X0))))
    C <- sweep(loo, 2, colMeans(loo), "-")
    V <- (n - 1) / n * crossprod(C)
    f <- if (tdnn) {
      np_tdnn(X = d$X, X0 = d$X0, Y = d$y, s1 = 2, s2 = 4, dim = 2,
        method = "asymp", normalize = normalize, return = character())
    } else {
      np_dnn(X = d$X, X0 = d$X0, Y = d$y, k = 2, s = 3,
        method = "asymp", normalize = normalize, return = character())
    }
    expect_equal(unname(vcov(f)), V, tolerance = 1e-12)
    expect_equal(predict(f, X0 = d$X0, se.fit = TRUE)$se.fit^2, diag(V), tolerance = 1e-12)
  }
  })
})

test_that("NP covariance keeps support-local variation with extreme unused outcomes", {
  X <- matrix(0:7, ncol = 1)
  y <- c(1e15 + c(0, 1, 3), rep(-1e300, 5))
  f <- np_dnn(X = X, X0 = matrix(c(0, .1), ncol = 1), Y = y, s = 7, return = character())
  ref <- .np_ref_cov(X, matrix(c(0, .1), ncol = 1), c(0, 1, 3, rep(0, 5)), 7)
  expect_equal(unname(vcov(f)), ref, tolerance = 1e-12)
  expect_gt(vcov(f)[1, 1], 0)
})

test_that("NP singular geometry and scalar/batch prediction are coherent", {
  X <- cbind(a = c(-2, -1, 0, 1, 2, 3), b = c(-4, -2, 0, 2, 4, 6))
  X0 <- rbind(c(0, 0), c(1, 2))
  f <- np_dnn(X = X, X0 = X0, Y = c(1, 3, -1, 4, 2, 5),
    distance = distance_mahalanobis(), s = 3, return = character())
  expect_equal(predict(f, X0 = X0), f$mu)
  for (j in 1:2) {
    expect_equal(predict(f, x0 = X0[j, ], se.fit = TRUE)$se.fit, f$se[j])
    expect_equal(unname(vcov(f)[j, j]), f$se[j]^2)
  }
  no_deletions <- np_dnn(X = X, X0 = X0, Y = 1:6, s = 6, se = FALSE)
  expect_equal(predict(no_deletions, X0 = X0), no_deletions$mu)
  expect_error(predict(no_deletions, X0 = X0, se.fit = TRUE), "s < n")
})

test_that("TDNN tabular saved state works without optional outputs or fitted SE", {
  .with_duplicate_policy_warnings({
  d <- .np_method_data()
  dat <- data.frame(a = d$X[, 1], b = d$X[, 2], y = d$y)
  f <- np_tdnn(data = dat, x_cols = 1:2, y_col = "y", X0 = d$X0,
    distance = distance_l1_scaled("custom", scale_vec = c(2, 5)),
    s1 = 2, s2 = 4, se = FALSE, return = character())
  f <- unserialize(serialize(f, NULL))
  nd <- data.frame(b = d$X0[, 2], y = 0, a = d$X0[, 1])
  expect_equal(predict(f, newdata = nd), f$mu)
  pr <- predict(f, newdata = nd, interval = "confidence", se.fit = TRUE)
  expect_true(all(is.finite(pr$se.fit)))
  expect_equal(unname(pr$fit[, 1]), f$mu)
  expect_null(f$W)
  expect_null(f$y)
  expect_error(confint(f), "se = TRUE")
  })
})

test_that("NP rejects corrupted saved model state before prediction", {
  .with_duplicate_policy_warnings({
  d <- .np_method_data()
  f <- np_dnn(X = d$X, X0 = d$X0, Y = d$y, s = 3, return = character())
  bad <- f
  bad$.model$y_vec <- bad$.model$y_vec[-1]
  expect_error(predict(bad, X0 = d$X0), "y_vec.*length")
  bad <- f
  bad$.model$x_mat[1, 1] <- Inf
  expect_error(predict(bad, X0 = d$X0), "x_mat.*non-finite")
  bad <- f
  bad$.model$k <- 4
  expect_error(predict(bad, X0 = d$X0), "k <= s")
  bad <- f
  bad$.model$method <- "other"
  expect_error(predict(bad, X0 = d$X0), "method.*exact")
  bad <- f
  bad$.model$normalize <- NA
  expect_error(predict(bad, X0 = d$X0), "normalize.*TRUE or FALSE")
  bad <- f
  bad$.model$kind <- "other"
  expect_error(predict(bad, X0 = d$X0), "kind")
  td <- np_tdnn(X = d$X, X0 = d$X0, Y = d$y, s1 = 2, s2 = 4)
  td$.model$s2 <- td$.model$s1
  expect_error(predict(td, X0 = d$X0), "s1 < s2")
  })
})

test_that("print.np_result runs without error", {
  .test_rng_scope({
    d <- .make_np_data(n = 16L, p = 2L, m = 2L)
    out <- np_dnn(
      X = d$X, X0 = d$X0, Y = d$y,
      k = 1L, s = 2L, se = FALSE
    )
    txt <- capture.output(print(out))
    expect_true(any(grepl("<np_result>", txt, fixed = TRUE)))
  }, seed = 123L)
})

test_that("print.tdnn_np_result runs and shows key fields", {
  .test_rng_scope({
    d   <- .make_tdnn_data(n = 20L, m = 2L)
    out <- np_tdnn(X = d$X, X0 = d$X0, Y = d$y, s1 = 1L, s2 = 2L, se = FALSE)
    txt <- capture.output(print(out))
    expect_true(any(grepl("<tdnn_np_result>", txt, fixed = TRUE)))
    expect_true(any(grepl("TDNN", txt)))
    expect_true(any(grepl("s1=", txt)))
    expect_true(any(grepl("s2=", txt)))
  }, seed = 123L)
})
