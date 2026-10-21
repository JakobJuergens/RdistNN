# Pin large decimal fixtures to binary64 bits: R builds can parse decimals
# with different intermediate precision. Exact score assertions stay exact.
.plqr_comparison_case <- function(kind) {
  x <- c(1, 2, 3, 5, 7, 11) * 0x1.249ad2594c37dp+332
  switch(kind,
    positive = list(D = matrix(x), theta = 0x1.c5416bb92e3e6p+730, below = TRUE),
    negative = list(D = matrix(-x), theta = 0x1.c5416bb92e3e6p+730, below = FALSE),
    mixed_positive = list(D = cbind(x, x), theta = c(0x1.33234de7ad7e3p+797, -0x1.33234de7ad7e3p+796), below = TRUE),
    mixed_negative = list(D = cbind(x, x), theta = c(0x1.33234de7ad7e3p+796, -0x1.33234de7ad7e3p+797), below = FALSE),
    partial_overflow = list(D = matrix(0x1.249ad2594c37dp+332, 6L, 5L),
      theta = c(0x1.f25c186a6f04cp+690, 0x1.f25c186a6f04cp+690, -0x1.f25c186a6f04cp+690, -0x1.f25c186a6f04cp+690, -0x1.f25c186a6f04cp+690), below = FALSE))
}

.plqr_comparison_data <- function(D, g = 0) {
  p <- ncol(D)
  d <- data.frame(x = seq_len(nrow(D)), Y = 0, g = g)
  d[paste0("D", seq_len(p))] <- as.data.frame(D)
  d[paste0("mD", seq_len(p))] <- rep(list(rep(0, nrow(D))), p)
  d
}

.plqr_comparison_score <- function(d, theta) {
  p <- length(theta)
  plqr_score_matrix(d, "Y", paste0("D", seq_len(p)), "g",
    paste0("mD", seq_len(p)), theta, .5)
}

.plqr_comparison_fit <- function(d, theta) {
  p <- length(theta)
  plqr_dnn(data = d, x_cols = "x", D_cols = paste0("D", seq_len(p)),
    g_col = "g", mD_cols = paste0("mD", seq_len(p)), y_col = "Y",
    theta_eval = theta, tau = .5, x0 = 0, s = 2L, return = c("W", "score"))
}

test_that("every PLQR constructor certifies decisive overflowing comparisons", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  for (native in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = native)
    if (native) expect_true(is.function(.rdistnn_compiled_fun("rustat_plqr_g_all")))
    for (kind in c("positive", "negative", "mixed_positive", "mixed_negative", "partial_overflow")) {
      a <- .plqr_comparison_case(kind)
      d <- .plqr_comparison_data(a$D)
      score <- (.5 - a$below) * a$D
      expected_w <- (6 - seq_len(6L)) / choose(6, 2)
      # Independent complete-subset deletion means, with ascending x selecting
      # the smaller retained observation index in each pair. Variance is
      # translation invariant; center before averaging huge constant scores.
      reference_score <- sweep(score, 2L, score[1L, ], "-")
      loo <- t(matrix(vapply(seq_len(6L), function(i) {
        sets <- utils::combn(setdiff(seq_len(6L), i), 2L, simplify = FALSE)
        rowMeans(matrix(vapply(sets, function(ids) reference_score[min(ids), ], numeric(ncol(score))),
          nrow = ncol(score)))
      }, numeric(ncol(score))), nrow = ncol(score), ncol = 6L))
      expected_se <- sqrt(5 / 6 * colSums(sweep(loo, 2L, colMeans(loo), "-")^2))
      expect_equal(unname(.plqr_comparison_score(d, a$theta)), unname(score), tolerance = 1e-14)
      fit <- .plqr_comparison_fit(d, a$theta)
      expect_equal(as.numeric(fit$W), expected_w, tolerance = 1e-14)
      expect_equal(as.numeric(fit$moment), as.numeric(crossprod(expected_w, score)), tolerance = 1e-14)
      expect_equal(as.numeric(fit$moment_se), expected_se, tolerance = 1e-14)
      expect_identical(fit$.model$below_threshold, rep(as.numeric(a$below), 6L))
      for (j in seq_len(ncol(score))) {
        sc <- score_plqr_component(0, "x", paste0("D", seq_len(ncol(score))), "g",
          paste0("mD", seq_len(ncol(score))), "Y", a$theta, j = j, tau = .5, s = 2L)
        ctx <- score_setup(sc, d)
        if (native) expect_identical(typeof(ctx$ptr_ctx), "externalptr") else expect_null(ctx$ptr_ctx)
        expect_equal(ctx$score, score[, j], tolerance = 0)
        expect_equal(score_eval(sc, c(1L, 3L), d, ctx), sum(score[c(1L, 3L), j]), tolerance = 1e-14)
      }
      restored <- unserialize(serialize(fit, NULL))
      expect_equal(restored$.model$below_threshold, fit$.model$below_threshold)
    }
  }
})

test_that("PLQR comparison risk decisions agree across ordinary matrix engines", {
  old <- options(matprod = "default", RdistNN.use_compiled = FALSE)
  on.exit(options(old), add = TRUE)
  a <- .plqr_comparison_case("partial_overflow")
  d <- .plqr_comparison_data(a$D)
  expected <- rep(0x1.249ad2594c37dp+331, nrow(a$D))
  for (product in c("default", "blas", "internal")) {
    options(matprod = product)
    expect_equal(.plqr_comparison_score(d, a$theta)[, 1L], expected, tolerance = 0)
    expect_identical(.plqr_comparison_fit(d, a$theta)$.model$below_threshold, rep(0, 6L))
  }
})

test_that("PLQR comparisons retain exact residual signs after huge cancellation", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  for (native in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = native)
    for (g in c(-1, -2^-1074, 0, 2^-1074)) {
      d <- .plqr_comparison_data(matrix(0x1.249ad2594c37dp+332, 6L, 2L), g = g)
      theta <- c(0x1.33234de7ad7e3p+797, -0x1.33234de7ad7e3p+797)
      expected <- rep(if (g >= 0) -0x1.249ad2594c37dp+331 else 0x1.249ad2594c37dp+331, 6L)
      expect_equal(.plqr_comparison_score(d, theta)[, 1L], expected, tolerance = 0)
      expect_identical(.plqr_comparison_fit(d, theta)$.model$below_threshold,
        rep(as.numeric(g >= 0), 6L))
      sc <- score_plqr_component(0, "x", c("D1", "D2"), "g", c("mD1", "mD2"),
        "Y", theta, j = 1L, tau = .5, s = 2L)
      expect_equal(score_setup(sc, d)$score, expected, tolerance = 0)
    }
  }
})

test_that("the common PLQR risk boundary preserves safe ordinary equalities", {
  safe <- .plqr_comparison_data(matrix(2^1020, 6L, 2L))
  boundary <- .plqr_comparison_data(matrix(2^1021, 6L, 2L))
  expect_equal(.plqr_comparison_score(safe, c(1, -1)),
    as.matrix(safe[c("D1", "D2")]) * -.5, tolerance = 0)
  expect_equal(.plqr_comparison_score(boundary, c(1, -1)),
    as.matrix(boundary[c("D1", "D2")]) * -.5, tolerance = 0)
})

test_that("ordinary PLQR cancellation decisions agree across matrix products", {
  old <- options(matprod = "default", RdistNN.use_compiled = FALSE)
  on.exit(options(old), add = TRUE)
  D <- matrix(rep(c(2^54, 1, -2^54), each = 6L), 6L)
  d <- .plqr_comparison_data(D)
  d$Y <- .5
  for (mode in c("default", "blas", "internal")) {
    options(matprod = mode)
    expect_equal(unname(.plqr_comparison_score(d, rep(1, 3L))), -.5 * D, tolerance = 0)
    expect_identical(.plqr_comparison_fit(d, rep(1, 3L))$.model$below_threshold, rep(1, 6L))
  }
})

test_that("PLQR comparisons preserve subnormal product signs and exact equality", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  below <- getFromNamespace(".plqr_below_threshold", "RdistNN")
  huge <- .Machine$double.xmax
  eta <- 2^-1074
  cases <- list(
    list(D = c(huge, huge, 1), theta = c(huge, -huge, 0), g = eta, expected = TRUE),
    list(D = c(huge, huge, 1), theta = c(huge, -huge, 0), g = -eta, expected = FALSE),
    list(D = c(huge, huge, 1), theta = c(huge, -huge, 0), g = 0, expected = TRUE),
    list(D = c(eta, 1), theta = c(eta, 0), g = 0, expected = TRUE),
    list(D = c(eta, 1), theta = c(-eta, 0), g = 0, expected = FALSE))
  for (native in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = native)
    for (a in cases) {
      D <- matrix(rep(a$D, each = 6L), 6L)
      d <- .plqr_comparison_data(D, g = a$g)
      expect_identical(below(D, a$theta, d$g, d$Y), rep(a$expected, 6L))
      j <- ncol(D)
      expected <- rep(.5 - a$expected, 6L)
      expect_identical(unname(.plqr_comparison_score(d, a$theta)[, j]), expected)
      sc <- score_plqr_component(0, "x", paste0("D", seq_len(j)), "g",
        paste0("mD", seq_len(j)), "Y", a$theta, j = j, tau = .5, s = 2L)
      ctx <- score_setup(sc, d)
      expect_identical(unname(ctx$score), expected)
    }
  }
})

test_that("the public PLQR score matrix validates computed scores", {
  d <- .plqr_comparison_data(matrix(.Machine$double.xmax, 6L, 1L))
  d$mD1 <- -.Machine$double.xmax
  expect_error(.plqr_comparison_score(d, 1), "Computed PLQR score.*[Rr]escale")
})

test_that("saved PLQR plot reconstruction uses the certified retained indicator", {
  skip_if_not_installed("ggplot2")
  a <- .plqr_comparison_case("partial_overflow")
  fit <- .plqr_comparison_fit(.plqr_comparison_data(a$D), a$theta)
  fit$score <- NULL
  fit <- unserialize(serialize(fit, NULL))
  tuning <- new_tuning("certified_comparison_grid", eval = function(...) list(s = 2L),
    eval_grid = function(...) data.frame(s = c(2L, 3L)))
  p <- plot_tuning(fit, tuning = tuning, se = FALSE)
  expect_equal(p$data$estimate, rep(0x1.249ad2594c37dp+331, nrow(p$data)), tolerance = 0)
  calibration <- plot_calibration(fit)
  expect_identical(calibration$data$share, 0)
})
