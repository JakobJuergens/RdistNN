test_that("NP marginal prediction uncertainty agrees with explicit deletions", {
  .with_localization_policy_warnings({
  X <- cbind(c(-2, -1, 0, 0, 1, 2, 4), c(3, 1, 0, 0, 3, 2, -1))
  X0 <- rbind(c(1.7, 1), c(-.2, .1), c(0, 0), c(1.7, 1))
  y <- c(2, -1, 3, 7, -2, 4, 1)
  scale <- c(2, 5)
  scaled_X <- sweep(X, 2, scale, "/")
  scaled_X0 <- sweep(X0, 2, scale, "/")

  # Exact weights enumerate subsamples; asymptotic weights use direct Beta
  # CDF differences. Neither reference uses the rank-deletion implementation.
  reference_weights <- function(x_mat, s, k, method) {
    n <- nrow(x_mat)
    W <- matrix(0, n, nrow(X0))
    if (method == "exact") {
      W <- .tie_reference_weights(x_mat, scaled_X0, k, s)
    } else {
      table <- rowMeans(vapply(seq_len(k), function(a) {
        diff(pbeta((0:n) / n, a, s - a + 1))
      }, numeric(n)))
      for (j in seq_len(nrow(X0))) {
        d <- rowSums(sweep(x_mat, 2, scaled_X0[j, ], "-")^2)
        W[, j] <- table[order(.tie_reference_order(d))]
      }
    }
    .block_oracle_average_weights(W, vapply(seq_len(nrow(scaled_X0)), function(j) sqrt(rowSums(sweep(x_mat, 2L, scaled_X0[j, ], "-")^2)), numeric(nrow(x_mat))))
  }

  for (tdnn in c(FALSE, TRUE)) for (method in c("exact", "asymp")) {
    reference_mean <- function(x_mat, y_vec) {
      W <- if (tdnn) {
        2 * reference_weights(x_mat, 4, 1, method) -
          reference_weights(x_mat, 2, 1, method)
      } else reference_weights(x_mat, 3, 2, method)
      as.numeric(crossprod(y_vec, W))
    }
    n <- nrow(X)
    loo <- t(vapply(seq_len(n), function(i) {
      reference_mean(scaled_X[-i, , drop = FALSE], y[-i])
    }, numeric(nrow(X0))))
    centered <- sweep(loo, 2, colMeans(loo), "-")
    V <- (n - 1) / n * crossprod(centered)
    expected_se <- sqrt(diag(V))
    expected_mu <- reference_mean(scaled_X, y)
    expected_interval <- cbind(fit = expected_mu,
      lwr = expected_mu - qnorm(.95) * expected_se,
      upr = expected_mu + qnorm(.95) * expected_se)
    rownames(expected_interval) <- paste0("target_", seq_len(nrow(X0)))

    for (normalize in c(FALSE, TRUE)) {
      args <- list(X = X, X0 = X0, Y = y, method = method,
        normalize = normalize, distance = distance_l2_scaled("custom", scale_vec = scale),
        return = character())
      fit <- if (tdnn) {
        do.call(np_tdnn, c(args, list(s1 = 2, s2 = 4, dim = 2)))
      } else do.call(np_dnn, c(args, list(s = 3, k = 2)))
      prediction <- predict(fit, X0 = X0, se.fit = TRUE,
        interval = "confidence", level = .9)
      expect_equal(prediction$fit, expected_interval, tolerance = 1e-12)
      expect_equal(prediction$se.fit, expected_se, tolerance = 1e-12)
      expect_equal(unname(vcov(fit)), V, tolerance = 1e-12)
      expect_equal(prediction$se.fit^2, unname(diag(vcov(fit))), tolerance = 1e-12)
      expect_equal(predict(fit, X0 = X0, interval = "confidence", level = .9),
        prediction$fit, tolerance = 1e-12)
      expect_equal(predict(fit, x0 = X0[2, ], se.fit = TRUE)$se.fit,
        expected_se[2], tolerance = 1e-12)
    }
  }
  })
})

.np_allocation_threshold <- function(n, m) 8 * n * m / 2

test_that("the allocation trace retains buffers forbidden by the prediction guard", {
  skip_if_not(isTRUE(capabilities("profmem")), "R memory profiling is unavailable")
  n <- 8L
  m <- 256L
  limit <- 8 * n * m * 8
  profile <- tempfile("np-allocation-control-", fileext = ".log")
  on.exit(unlink(profile), add = TRUE)
  on.exit(Rprofmem(NULL), add = TRUE)
  Rprofmem(profile, threshold = .np_allocation_threshold(n, m))
  # A permitted n-by-m work buffer and a forbidden m-by-m covariance buffer.
  allowed <- matrix(0, n, m)
  forbidden <- matrix(0, m, m)
  Rprofmem(NULL)
  lines <- readLines(profile, warn = FALSE)
  allocations <- as.numeric(sub(" .*", "", lines[grepl("^[0-9]+", lines)]))
  expect_true(any(allocations >= 8 * n * m & allocations < limit))
  expect_true(any(allocations >= limit))
})

test_that("NP prediction uncertainty avoids quadratic target allocations", {
  skip_if_not(isTRUE(capabilities("profmem")), "R memory profiling is unavailable")
  X <- matrix(1:20, ncol = 1)
  X0 <- matrix(seq(1, 20, length.out = 2000), ncol = 1)
  y <- sin(1:20)
  fits <- list(
    dnn = np_dnn(X = X, x0 = 0, Y = y, s = 5, se = FALSE, return = character()),
    tdnn = np_tdnn(X = X, x0 = 0, Y = y, s1 = 3, s2 = 7,
      se = FALSE, return = character())
  )
  profile <- tempfile("np-prediction-memory-", fileext = ".log")
  on.exit(unlink(profile), add = TRUE)
  on.exit(Rprofmem(NULL), add = TRUE)

  for (fit in fits) for (interval in c("none", "confidence")) {
    # Load lazy code before measuring the target-dependent allocations.
    predict(fit, X0 = X0[1:2, , drop = FALSE], se.fit = TRUE, interval = interval)
    # Retain large work buffers and every allocation at the rejection limit,
    # without recording the many smaller allocations and their call stacks.
    Rprofmem(profile, threshold = .np_allocation_threshold(nrow(X), nrow(X0)))
    prediction <- tryCatch(
      predict(fit, X0 = X0, se.fit = interval == "none", interval = interval),
      finally = Rprofmem(NULL)
    )
    lines <- readLines(profile, warn = FALSE)
    allocations <- as.numeric(sub(" .*", "", lines[grepl("^[0-9]+", lines)]))
    expect_true(length(allocations) > 0L)
    # Permit several n-by-m work buffers, but reject the 32 MB m-by-m
    # covariance allocation (over 100 times the 320 KB weight matrix).
    expect_lt(max(allocations), 8 * nrow(X) * nrow(X0) * 8)
    if (interval == "none") {
      expect_length(prediction$fit, nrow(X0))
      expect_true(all(is.finite(prediction$se.fit)))
    } else {
      expect_equal(dim(prediction), c(nrow(X0), 3L))
      expect_true(all(is.finite(prediction)))
    }
  }
})
