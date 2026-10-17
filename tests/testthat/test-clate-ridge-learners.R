# Dense independent normal equations, deliberately separate from the block engine.
.clate_reference_dictionary <- function(d, M) {
  J <- M^4
  out <- matrix(0, nrow(d), 3L * J)
  for (i in seq_len(nrow(d))) {
    z <- as.numeric(d[i, paste0("Z", 1:4)])
    cell <- 1L + sum(pmin(floor(z * M), M - 1L) * M^(0:3))
    q <- c(1, sqrt(12) * (d$X1[i] - .5), sqrt(12) * (d$X2[i] - .5))
    out[i, 3L * (cell - 1L) + 1:3] <- sqrt(J) * q
  }
  out
}

test_that("ridge-series adapter solves the stated ambient-normalized objectives", {
  d <- clate_demo_data(n = 18L, seed = 331L)$data
  d[1L, paste0("Z", 1:4)] <- 1
  d[2L, paste0("Z", 1:4)] <- 0
  M <- 2L
  lambda <- .1
  p <- .clate_reference_dictionary(d, M)
  for (deleted in list(integer(), 3L, c(2L, 7L))) {
    fit <- .clate_ridge_predictions(d, M, lambda, nrow(d), delete = deleted)
    retain <- setdiff(seq_len(nrow(d)), deleted)
    weights <- cbind(e = rep(1, nrow(d)), m_y0 = 1 - d$A, m_y1 = d$A,
                     m_d0 = 1 - d$A, m_d1 = d$A)
    responses <- cbind(e = d$A, m_y0 = d$Y, m_y1 = d$Y, m_d0 = d$D, m_d1 = d$D)
    for (j in seq_len(5L)) {
      L <- weights[retain, j]
      P <- p[retain, , drop = FALSE]
      gram <- crossprod(P, L * P) / nrow(d) + lambda * diag(ncol(p))
      rhs <- crossprod(P, L * responses[retain, j]) / nrow(d)
      beta <- as.numeric(solve(gram, rhs))
      expect_equal(fit$coefficients[, j], beta, ignore_attr = TRUE, tolerance = 1e-11)
      raw <- as.numeric(p %*% beta)
      used <- switch(j, pmin(.75, pmax(.25, raw)),
                     pmin(1, pmax(-1, raw)), pmin(1, pmax(-1, raw)),
                     pmin(1, pmax(0, raw)), pmin(1, pmax(0, raw)))
      expect_equal(fit$predictions[[j]], used, tolerance = 1e-11)
    }
    expect_equal(fit$training$normalizer_n, nrow(d))
    expect_equal(fit$training$delete, deleted)
    expect_identical(fit$row_ids, d$row_id)
  }
})

test_that("all one- and two-deleted ridge fits obey the exact rank-one identity", {
  d <- clate_demo_data(n = 7L, seed = 77L)$data
  p <- .clate_reference_dictionary(d, 1L)
  lambda <- .3
  fits <- lapply(c(list(integer()), as.list(seq_len(7L))), function(I) {
    .clate_ridge_predictions(d, 1L, lambda, 7L, delete = I)
  })
  for (I in c(list(integer()), as.list(seq_len(7L)))) {
    full <- fits[[if (!length(I)) 1L else I + 1L]]
    for (j in setdiff(seq_len(7L), I)) {
      retained <- .clate_ridge_predictions(d, 1L, lambda, 7L, delete = c(I, j))
      keep <- setdiff(seq_len(7L), I)
      L <- d$A
      G <- crossprod(p[keep, , drop = FALSE], L[keep] * p[keep, , drop = FALSE]) / 7 + lambda * diag(3)
      rhs <- L[j] * p[j, ] / 7 * (d$Y[j] - sum(p[j, ] * retained$coefficients[, "m_y1"]))
      expect_equal(
        full$coefficients[, "m_y1"] - retained$coefficients[, "m_y1"],
        as.numeric(solve(G, rhs)), ignore_attr = TRUE, tolerance = 1e-11
      )
    }
  }
})

test_that("empty arms cells and deleted samples still have unique bounded ridge predictions", {
  d <- clate_demo_data(n = 2L, seed = 71L)$data
  d$A <- 0
  f <- .clate_ridge_predictions(d, 2L, .1, 2L)
  expect_equal(f$predictions$m_y1, c(0, 0))
  expect_equal(f$predictions$m_d1, c(0, 0))
  empty <- .clate_ridge_predictions(d, 2L, .1, 2L, delete = 1:2)
  expect_equal(empty$coefficients, matrix(0, 48L, 5L), ignore_attr = TRUE)
  expect_equal(empty$predictions$e, c(.25, .25))
  for (bad in list(0, -1, NA_real_, Inf)) {
    expect_error(.clate_ridge_predictions(d, 1L, bad, 2L), "lambda")
  }
  expect_error(.clate_ridge_predictions(d, 0, .1, 2L), "M")
  expect_error(.clate_ridge_predictions(d, 1L, .1, 0), "normalizer_n")
  expect_error(.clate_ridge_predictions(d, 1L, .1, 2L, delete = c(1, 1)), "delete")
  expect_error(.clate_ridge_predictions(d, 1L, .1, 2L, delete = 3L), "delete")
  d$Y[1L] <- 2
  expect_error(.clate_ridge_predictions(d, 1L, .1, 2L), "Y")
})

test_that("bounded ridge predictions feed the public full-sample workflow reproducibly", {
  for (heterogeneous in c(FALSE, TRUE)) {
    z <- clate_demo_data(n = 90L, seed = 231L, heterogeneous = heterogeneous)
    d <- z$data
    learner <- .clate_ridge_predictions(d, M = 1L, lambda = 1 / sqrt(nrow(d)))
    at <- rbind(c(.25, .25), c(.5, .5), c(.75, .75))
    arguments <- c(list(
      X = as.matrix(d[c("X1", "X2")]), Y = d$Y, D = d$D, A = d$A,
      X0 = at, s = 12L, training = "full_sample", row_ids = learner$row_ids
    ), learner$predictions)
    fit <- do.call(clate_dnn, arguments)
    expect_identical(fit$status$ratio_status, rep("ok", 3L))
    expect_true(all(is.finite(fit$theta)))
    expect_true(all(is.finite(fit$se)))
    expect_equal(do.call(clate_dnn, arguments), fit)
    expect_identical(fit$meta$training, "full_sample")
    expect_identical(fit$row_ids, d$row_id)
  }
})

test_that("ridge predictions remain aligned by observation identity after row permutation", {
  d <- clate_demo_data(n = 18L, seed = 331L)$data
  d$row_id <- paste0("observation_", c(83, 14, 65, 21, 98, 37, 42, 56, 101:110))
  permutation <- c(18:10, 1:9)
  reordered <- d[permutation, ]
  original_delete <- c(2L, 7L)
  reordered_delete <- match(d$row_id[original_delete], reordered$row_id)
  original <- .clate_ridge_predictions(d, 2L, .1, nrow(d), delete = original_delete)
  shuffled <- .clate_ridge_predictions(reordered, 2L, .1, nrow(d), delete = reordered_delete)
  expect_identical(shuffled$row_ids, reordered$row_id)
  join <- match(original$row_ids, shuffled$row_ids)
  for (nm in names(original$predictions)) {
    expect_equal(shuffled$predictions[[nm]][join], original$predictions[[nm]], tolerance = 1e-12)
    expect_equal(shuffled$raw_predictions[[nm]][join], original$raw_predictions[[nm]], tolerance = 1e-12)
  }
  aligned <- lapply(shuffled$predictions, function(x) x[join])
  arguments <- list(Y = d$Y, D = d$D, A = d$A, row_ids = d$row_id, training = "full_sample")
  signals <- do.call(iv_augmented_signals, c(arguments, aligned))
  direct <- do.call(iv_augmented_signals, c(arguments, original$predictions))
  expect_equal(signals$phi_y, direct$phi_y)
  expect_equal(signals$phi_d, direct$phi_d)
  expect_identical(signals$row_ids, d$row_id)
})
