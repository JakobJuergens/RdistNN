# These regressions distinguish subtract-then-scale from pre-scaling inputs,
# and sequential double accumulation from rowSums()'s extended accumulator.
scaled_arithmetic_reference <- function(X, X0, inv_s, metric, squared = FALSE) {
  out <- matrix(0, nrow(X), nrow(X0))
  for (q in seq_len(nrow(X0))) {
    for (i in seq_len(nrow(X))) {
      acc <- 0
      for (j in seq_len(ncol(X))) {
        d <- (X[i, j] - X0[q, j]) * inv_s[j]
        acc <- acc + if (metric == "l1") abs(d) else d * d
      }
      out[i, q] <- if (metric == "l2" && !squared) sqrt(acc) else acc
    }
  }
  out
}

check_scaled_arithmetic <- function(X, X0, metric, rule, scale_vec = NULL,
                                    inv_s, squared = FALSE, hand = NULL) {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  prefix <- paste0("rustat_dist_", metric, "_scaled")
  cpp <- .rdistnn_compiled_fun(prefix)
  cpp_batch <- .rdistnn_compiled_fun(paste0(prefix, "_batch"))
  # Missing native code must fail this suite instead of silently testing R twice.
  expect_true(is.function(cpp))
  expect_true(is.function(cpp_batch))
  if (!is.function(cpp) || !is.function(cpp_batch)) {
    stop("Scaled-distance regression requires the compiled package.")
  }
  flags <- if (metric == "l2") list(squared = squared) else list()
  spec <- do.call(get(paste0("distance_", metric, "_scaled")),
                  c(list(rule = rule, scale_vec = scale_vec, eps = 1e-8), flags))
  ctx <- spec$setup(X)
  expect_equal(ctx$inv_s, inv_s, tolerance = 0)
  literal <- scaled_arithmetic_reference(X, X0, inv_s, metric, squared)
  native <- do.call(cpp_batch, c(list(X, X0, inv_s), flags))
  expect_equal(native, literal, tolerance = 1e-15)
  if (!is.null(hand)) expect_equal(native[, 1L], hand, tolerance = 1e-15)
  rank_columns <- function(D) apply(D, 2L, rank, ties.method = "first")

  for (compiled in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = compiled)
    if (!compiled) {
      expect_null(.rdistnn_compiled_fun(prefix))
      expect_null(.rdistnn_compiled_fun(paste0(prefix, "_batch")))
    }
    for (context in list(ctx, list(inv_s = inv_s), NULL)) {
      batch <- spec$eval_batch(X, X0, context)
      expect_equal(batch, native, tolerance = 1e-14)
      expect_identical(rank_columns(batch), rank_columns(native))
      scalar <- vapply(seq_len(nrow(X0)), function(q) {
        spec$eval(X, X0[q, ], context)
      }, numeric(nrow(X)))
      expect_equal(scalar, native, tolerance = 1e-14)
      expect_identical(rank_columns(scalar), rank_columns(native))
    }
    for (q in seq_len(nrow(X0))) {
      dx <- spec$eval(X, X0[q, ], list(inv_s = inv_s))
      ordered <- .distance_ranking(X, X0[q, ], spec)
      expect_equal(dx, native[, q], tolerance = 1e-14)
      expect_identical(.tie_reference_order(dx), ordered$ord)
      expect_equal(ordered$dx, native[, q], tolerance = 1e-14)
    }
    # Reuse fitted scales on a new training subset; the cached scaled X cannot
    # replace the caller's X, and the subset must not refit the scales.
    keep <- c(2L, 4L)
    expect_equal(spec$eval_batch(X[keep, , drop = FALSE], X0, ctx),
                 native[keep, , drop = FALSE], tolerance = 1e-14)
  }
}

test_that("scaled distances preserve small differences at large offsets", {
  # At 1e12 these increments are exactly representable multiples of 2^-13.
  X <- 1e12 + rbind(c(8, 16), c(16, 8), c(25, 0), c(0, 25), c(0, 0)) / 8192
  X0 <- 1e12 + rbind(c(0, 0), c(8, 8), c(0, 0)) / 8192
  hand_l1 <- c(8 / 3 + 16 / 7, 16 / 3 + 8 / 7, 25 / 3, 25 / 7, 0) / 8192
  hand_l2_sq <- c((8 / 3)^2 + (16 / 7)^2, (16 / 3)^2 + (8 / 7)^2,
                  (25 / 3)^2, (25 / 7)^2, 0) / 8192^2
  check_scaled_arithmetic(X, X0, "l1", "custom", c(3, 7), c(1 / 3, 1 / 7),
                          hand = hand_l1)
  for (squared in c(FALSE, TRUE)) {
    check_scaled_arithmetic(X, X0, "l2", "custom", c(3, 7), c(1 / 3, 1 / 7),
                            squared, if (squared) hand_l2_sq else sqrt(hand_l2_sq))
  }
})

test_that("scaled ranks retain stable ties with floored scales and signed targets", {
  X <- cbind(c(-3, -1, 0, 0, 2, 8), c(1, 4, 2, 2, 9, 15), rep(5, 6))
  X0 <- rbind(c(0, 2, 5), c(1, 3, 6), c(-1, -3, -6), c(8, 15, 5))
  for (metric in c("l1", "l2")) {
    for (squared in if (metric == "l2") c(FALSE, TRUE) else FALSE) {
      check_scaled_arithmetic(X, X0, metric, "custom", c(-2, 0, 4),
                              c(0.5, 1e8, 0.25), squared)
      check_scaled_arithmetic(X, X0, metric, "range", inv_s = c(1 / 11, 1 / 14, 1e8),
                              squared = squared)
    }
  }
})

test_that("unscaled ranks use the same coordinate accumulation as native code", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  X0 <- matrix(0, 2L, 3L)
  for (metric in c("l1", "l2")) {
    X <- rbind(c(if (metric == "l1") 1e16 else 1e8, 1, 1),
               c(if (metric == "l1") 1e16 else 1e8, 0, 0))
    for (squared in if (metric == "l2") c(FALSE, TRUE) else FALSE) {
      options(RdistNN.use_compiled = TRUE)
      symbols <- paste0("rustat_", c(paste0("dist_", metric),
                                      paste0("dist_", metric, "_batch")))
      native_funs <- lapply(symbols, .rdistnn_compiled_fun)
      expect_true(all(vapply(native_funs, is.function, logical(1))))
      if (!all(vapply(native_funs, is.function, logical(1)))) {
        stop("Distance regression requires the compiled package.")
      }
      flags <- if (metric == "l2") list(squared = squared) else list()
      native <- do.call(native_funs[[2L]], c(list(X, X0), flags))
      # Each small contribution rounds away when added to the first coordinate.
      hand <- if (metric == "l1" || squared) 1e16 else 1e8
      expect_equal(native, matrix(hand, 2L, 2L), tolerance = 0)
      spec <- do.call(get(paste0("distance_", metric)), flags)
      for (compiled in c(TRUE, FALSE)) {
        options(RdistNN.use_compiled = compiled)
        if (!compiled) {
          for (symbol in symbols) expect_null(.rdistnn_compiled_fun(symbol))
        }
        expect_equal(spec$eval_batch(X, X0), native, tolerance = 1e-15)
        dx <- spec$eval(X, X0[1L, ])
        expect_equal(dx, native[, 1L], tolerance = 1e-15)
        expect_identical(order(dx, seq_along(dx)), 1:2)
        ordered <- .distance_ranking(X, X0[1L, ], spec)
        expect_identical(ordered$ord, 1:2)
      }
    }
  }
})

test_that("L2 rounds each square before addition even when FMA is available", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  a <- 0.5790338377234229
  b <- 1.4758492921057114
  unit_X <- rbind(c(a, b), c(b, a))
  X0 <- matrix(0, 2L, 2L)
  # The two rounded products are identical up to permutation. Fusing the
  # second square with addition instead gives distinct native distances.
  a_sq <- a * a
  b_sq <- b * b
  hand_sq <- a_sq + b_sq
  for (scaled in c(FALSE, TRUE)) {
    X <- if (scaled) sweep(unit_X, 2L, c(2, 4), "*") else unit_X
    inv_s <- c(0.5, 0.25)
    suffix <- if (scaled) "_scaled" else ""
    for (squared in c(FALSE, TRUE)) {
      options(RdistNN.use_compiled = TRUE)
      symbols <- paste0("rustat_", c(paste0("dist_l2", suffix),
                                      paste0("dist_l2", suffix, "_batch")))
      native_funs <- lapply(symbols, .rdistnn_compiled_fun)
      expect_true(all(vapply(native_funs, is.function, logical(1))))
      if (!all(vapply(native_funs, is.function, logical(1)))) {
        stop("FMA regression requires the compiled L2 package functions.")
      }
      spec <- if (scaled) {
        distance_l2_scaled("custom", scale_vec = c(2, 4), squared = squared)
      } else {
        distance_l2(squared)
      }
      ctx <- spec$setup(X)
      args <- list(X, X0[1L, ])
      batch_args <- list(X, X0)
      if (scaled) {
        args <- c(args, list(inv_s))
        batch_args <- c(batch_args, list(inv_s))
      }
      native <- do.call(native_funs[[1L]], c(args, list(squared)))
      native_batch <- do.call(native_funs[[2L]], c(batch_args, list(squared)))
      hand <- if (squared) hand_sq else sqrt(hand_sq)
      expect_equal(native, rep(hand, 2L), tolerance = 0)
      expect_equal(native_batch, matrix(hand, 2L, 2L), tolerance = 0)
      for (compiled in c(TRUE, FALSE)) {
        options(RdistNN.use_compiled = compiled)
        if (!compiled) {
          for (symbol in symbols) expect_null(.rdistnn_compiled_fun(symbol))
        }
        dx <- spec$eval(X, X0[1L, ], ctx)
        batch <- spec$eval_batch(X, X0, ctx)
        expect_equal(dx, rep(hand, 2L), tolerance = 0)
        expect_equal(batch, matrix(hand, 2L, 2L), tolerance = 0)
        expect_identical(order(dx, seq_along(dx)), 1:2)
        public_order <- if (scaled) {
          .distance_ranking(X, X0[1L, ], distance_l2_scaled(rule = "custom", scale_vec = c(2, 4), squared = squared))
        } else {
          .distance_ranking(X, X0[1L, ], distance_l2(squared = squared))
        }
        expect_identical(public_order$ord, 1:2)
      }
    }
  }
})
