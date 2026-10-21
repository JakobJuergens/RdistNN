# Pin large decimal fixtures to binary64 bits: R builds can parse decimals
# with different intermediate precision. Exact score assertions stay exact.
test_that("native PLQR comparisons recover exact overflow signs", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  create <- .rdistnn_compiled_fun("rustat_plqr_ctx_create")
  all <- .rdistnn_compiled_fun("rustat_plqr_g_all")
  indexed <- .rdistnn_compiled_fun("rustat_plqr_g_indexed")
  expect_true(is.function(create))
  expect_true(is.function(all))
  expect_true(is.function(indexed))
  # Exact Fraction.from_float references, recorded independently in the audit,
  # give respectively negative, positive, negative and positive thresholds.
  cases <- list(
    list(D = rep(0x1.249ad2594c37dp+332, 5), theta = c(1, 1, -1, -1, -1) * 0x1.f25c186a6f04cp+690, below = 0),
    list(D = rep(0x1.249ad2594c37dp+332, 5), theta = c(-1, -1, 1, 1, 1) * 0x1.f25c186a6f04cp+690, below = 1),
    list(D = rep(0x1.1ccf385ebc8a0p+1023, 3), theta = c(2, -2, -1), below = 0),
    list(D = rep(0x1.1ccf385ebc8a0p+1023, 3), theta = c(2, -2, 1), below = 1),
    list(D = 0x1.249ad2594c37dp+332, theta = 0x1.f25c186a6f04cp+690, below = 1),
    list(D = -0x1.249ad2594c37dp+332, theta = 0x1.f25c186a6f04cp+690, below = 0)
  )
  for (case in cases) {
    D <- matrix(case$D, 1L)
    expected <- (.5 - case$below) * D[1L, 1L]
    ctx <- create(0, D, 0, matrix(0, 1L, ncol(D)), case$theta, .5, 1L)
    expect_identical(all(ctx), expected)
    expect_identical(indexed(c(1L, 1L), ctx), rep(expected, 2L))
  }
})

test_that("native PLQR comparisons retain exact cancellation across the binary range", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  create <- .rdistnn_compiled_fun("rustat_plqr_ctx_create")
  all <- .rdistnn_compiled_fun("rustat_plqr_g_all")
  indexed <- .rdistnn_compiled_fun("rustat_plqr_g_indexed")
  tiny <- .Machine$double.xmin * .Machine$double.eps
  cases <- list(
    list(D = c(0x1.249ad2594c37dp+332, 0x1.249ad2594c37dp+332), theta = c(0x1.f25c186a6f04cp+690, -0x1.f25c186a6f04cp+690), g = 1, y = 0, below = 1),
    list(D = c(0x1.1ccf385ebc8a0p+1023, 0x1.1ccf385ebc8a0p+1023), theta = c(2, -2), g = 0, y = 0, below = 1),
    list(D = 1, theta = 1, g = 2^1022, y = 2^1022, below = 1),
    list(D = 1, theta = 0, g = 2^1021, y = 2^1021, below = 1),
    list(D = c(0x1.249ad2594c37dp+332, 0x1.249ad2594c37dp+332), theta = c(0x1.33234de7ad7e3p+797, -0x1.33234de7ad7e3p+797), g = tiny, y = 0, below = 1),
    list(D = c(0x1.249ad2594c37dp+332, 0x1.249ad2594c37dp+332), theta = c(0x1.33234de7ad7e3p+797, -0x1.33234de7ad7e3p+797), g = -tiny, y = 0, below = 0),
    # The last exact product is smaller than the smallest representable double.
    list(D = c(2^1022, 2^1022, tiny), theta = c(2^1022, -2^1022, tiny), g = 0, y = 0, below = 1),
    list(D = c(2^1022, 2^1022, tiny), theta = c(2^1022, -2^1022, -tiny), g = 0, y = 0, below = 0),
    list(D = 1, theta = 0, g = 0, y = 0, below = 1)
  )
  for (case in cases) {
    D <- matrix(case$D, 1L)
    expected <- (.5 - case$below) * D[1L, 1L]
    ctx <- create(case$y, D, case$g, matrix(0, 1L, ncol(D)), case$theta, .5, 1L)
    expect_identical(all(ctx), expected)
    expect_identical(indexed(c(1L, 1L), ctx), rep(expected, 2L))
  }
})

test_that("native ordinary cancellation uses an exact adaptive comparison", {
  old <- options(RdistNN.use_compiled = TRUE, RdistNN.num_threads = 2L)
  on.exit(options(old), add = TRUE)
  permutations <- rbind(c(1L, 2L, 3L), c(1L, 3L, 2L), c(2L, 1L, 3L),
    c(2L, 3L, 1L), c(3L, 1L, 2L), c(3L, 2L, 1L))
  for (n in c(1L, 32L, 1024L)) {
    D <- matrix(rep(c(2^54, 1, -2^54), each = n), n, 3L)
    for (i in seq_len(nrow(permutations))) {
      at <- permutations[i, ]
      ctx <- .rdistnn_compiled_fun("rustat_plqr_ctx_create")(rep(.5, n),
        D[, at, drop = FALSE], rep(0, n), matrix(0, n, 3L), rep(1, 3L), .5,
        as.integer(match(2L, at)))
      # Exact threshold1 exceeds Y=.5 in every algebraically equivalent order.
      expect_identical(.rdistnn_compiled_fun("rustat_plqr_g_all")(ctx), rep(-.5, n))
      expect_identical(.rdistnn_compiled_fun("rustat_plqr_g_indexed")(c(n, 1L), ctx), rep(-.5, 2L))
    }
  }
})

test_that("native exceptional PLQR comparisons are safe across parallel row shapes", {
  old <- options(RdistNN.use_compiled = TRUE, RdistNN.num_threads = 2L)
  on.exit(options(old), add = TRUE)
  D <- matrix(0x1.249ad2594c37dp+332, 1024L, 5L)
  theta <- c(1, 1, -1, -1, -1) * 0x1.f25c186a6f04cp+690
  create <- .rdistnn_compiled_fun("rustat_plqr_ctx_create")
  ctx <- create(rep(0, nrow(D)), D, rep(0, nrow(D)), matrix(0, nrow(D), ncol(D)),
    theta, .5, 1L)
  expect_identical(.rdistnn_compiled_fun("rustat_plqr_g_all")(ctx), rep(0x1.249ad2594c37dp+331, 1024L))
  ids <- rep(c(1024L, 1L), 512L)
  expect_identical(.rdistnn_compiled_fun("rustat_plqr_g_indexed")(ids, ctx), rep(0x1.249ad2594c37dp+331, 1024L))
})
