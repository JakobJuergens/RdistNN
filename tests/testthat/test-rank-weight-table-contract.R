skip_if_not(
  all(vapply(
    c(
      "rank_weight_eval",
      "rank_weight_eval_batch",
      "rank_weight_dnn_exact",
      "rank_weight_dnn_asymp",
      "rank_weight_tdnn_exact",
      "rank_weight_tdnn_asymp"
    ),
    exists, logical(1),
    mode = "function"
  )),
  "rank-weight helpers not found; skipping file"
)

.layer1_rank_matrix <- function(n, m, seed = 101L) {
  set.seed(seed)
  R <- replicate(m, sample.int(n))
  storage.mode(R) <- "integer"
  R
}

.layer1_reference_batch <- function(R, ref_one) {
  n <- nrow(R)
  vapply(seq_len(ncol(R)), function(j) ref_one(R[, j]), numeric(n))
}

.layer1_expect_table_contract <- function(spec, n, R, ref_table, ref_one) {
  tab <- RdistNN:::.rank_weight_table(spec, n)

  expect_type(tab, "double")
  expect_length(tab, n)
  expect_true(all(is.finite(tab)))
  expect_equal(tab, ref_table, tolerance = 1e-12)

  W_table <- matrix(tab[R], nrow = n, ncol = ncol(R))
  expect_equal(rank_weight_eval_batch(spec, R), W_table, tolerance = 1e-12)
  expect_equal(rank_weight_eval_batch(spec, R), .layer1_reference_batch(R, ref_one), tolerance = 1e-12)
  expect_equal(rank_weight_eval(spec, R[, 1L]), ref_one(R[, 1L]), tolerance = 1e-12)
}

test_that("DNN rank-weight specs expose validated tables matching existing evaluators", {
  n <- 12L
  R <- .layer1_rank_matrix(n, m = 4L, seed = 202L)

  .layer1_expect_table_contract(
    spec = rank_weight_dnn_exact(k = 3L, s = 3L),
    n = n,
    R = R,
    ref_table = dnn_weights_exact_full_from_ranks(seq_len(n), k = 3L, s = 3L)$w_rank,
    ref_one = function(r) dnn_weights_exact_full_from_ranks(r, k = 3L, s = 3L)$w
  )

  .layer1_expect_table_contract(
    spec = rank_weight_dnn_exact(k = 2L, s = n),
    n = n,
    R = R,
    ref_table = dnn_weights_exact_full_from_ranks(seq_len(n), k = 2L, s = n)$w_rank,
    ref_one = function(r) dnn_weights_exact_full_from_ranks(r, k = 2L, s = n)$w
  )

  .layer1_expect_table_contract(
    spec = rank_weight_dnn_asymp(k = 2L, s = 5L, normalize = TRUE),
    n = n,
    R = R,
    ref_table = rank_weight_dnn_asymp(k = 2L, s = 5L, normalize = TRUE)$eval(seq_len(n)),
    ref_one = function(r) rank_weight_dnn_asymp(k = 2L, s = 5L, normalize = TRUE)$eval(r)
  )

  .layer1_expect_table_contract(
    spec = rank_weight_dnn_asymp(k = 2L, s = 5L, normalize = FALSE),
    n = n,
    R = R,
    ref_table = rank_weight_dnn_asymp(k = 2L, s = 5L, normalize = FALSE)$eval(seq_len(n)),
    ref_one = function(r) rank_weight_dnn_asymp(k = 2L, s = 5L, normalize = FALSE)$eval(r)
  )
})

test_that("TDNN rank-weight specs expose validated tables matching existing evaluators", {
  n <- 16L
  R <- .layer1_rank_matrix(n, m = 3L, seed = 303L)

  .layer1_expect_table_contract(
    spec = rank_weight_tdnn_exact(s1 = 3L, s2 = 7L, dim = 2L),
    n = n,
    R = R,
    ref_table = tdnn_rank_table_exact(n, s1 = 3L, s2 = 7L, dim = 2L),
    ref_one = function(r) dnn_weights_tdnn_exact_full_from_ranks(r, s1 = 3L, s2 = 7L, dim = 2L)$w
  )

  .layer1_expect_table_contract(
    spec = rank_weight_tdnn_asymp(s1 = 3L, s2 = 7L, dim = 2L, normalize = TRUE),
    n = n,
    R = R,
    ref_table = tdnn_rank_table_asymp(n, s1 = 3L, s2 = 7L, dim = 2L, normalize = TRUE),
    ref_one = function(r) {
      dnn_weights_tdnn_asymp_full_from_ranks(r, s1 = 3L, s2 = 7L, dim = 2L, normalize = TRUE)$w
    }
  )

  .layer1_expect_table_contract(
    spec = rank_weight_tdnn_asymp(s1 = 3L, s2 = 7L, dim = 2L, normalize = FALSE),
    n = n,
    R = R,
    ref_table = tdnn_rank_table_asymp(n, s1 = 3L, s2 = 7L, dim = 2L, normalize = FALSE),
    ref_one = function(r) {
      dnn_weights_tdnn_asymp_full_from_ranks(r, s1 = 3L, s2 = 7L, dim = 2L, normalize = FALSE)$w
    }
  )
})

test_that("rank_weight_eval_batch prefers table path when present", {
  calls <- new.env(parent = emptyenv())
  calls$table <- 0L
  calls$eval <- 0L
  calls$batch <- 0L

  spec <- list(
    name = "table_preferred",
    table = function(n) {
      calls$table <- calls$table + 1L
      w <- 1 / seq_len(n)
      w / sum(w)
    },
    eval = function(r, ctx = NULL) {
      calls$eval <- calls$eval + 1L
      rep(NA_real_, length(r))
    },
    eval_batch = function(R, ctx = NULL) {
      calls$batch <- calls$batch + 1L
      matrix(NA_real_, nrow = nrow(R), ncol = ncol(R))
    },
    .__kind = "rank_weight"
  )
  class(spec) <- c("rank_weight", class(spec))

  R <- .layer1_rank_matrix(n = 7L, m = 2L, seed = 404L)
  tab <- spec$table(nrow(R))
  calls$table <- 0L

  expect_equal(rank_weight_eval_batch(spec, R), matrix(tab[R], nrow = nrow(R)), tolerance = 1e-12)
  expect_equal(calls$table, 1L)
  expect_equal(calls$eval, 0L)
  expect_equal(calls$batch, 0L)

  expect_equal(rank_weight_eval(spec, R[, 1L]), tab[R[, 1L]], tolerance = 1e-12)
  expect_equal(calls$table, 2L)
  expect_equal(calls$eval, 0L)
  expect_equal(calls$batch, 0L)
})

test_that("rank-weight validation accepts table-only specs through the registry", {
  nm <- paste0("table_only_", as.integer(runif(1, 1e6, 2e6)))
  spec <- list(
    name = nm,
    table = function(n) {
      w <- seq_len(n)
      w / sum(w)
    }
  )

  register_rank(nm, spec)
  got <- get_rank(nm)
  expect_true(is_rank_weight(got))

  R <- .layer1_rank_matrix(n = 6L, m = 2L, seed = 454L)
  tab <- got$table(nrow(R))
  expect_equal(rank_weight_eval_batch(got, R), matrix(tab[R], nrow = nrow(R)), tolerance = 1e-12)
  expect_equal(rank_weight_eval(got, R[, 1L]), tab[R[, 1L]], tolerance = 1e-12)
})

test_that("custom rank weights without tables keep the old per-column fallback", {
  calls <- new.env(parent = emptyenv())
  calls$eval <- 0L

  spec <- as_rank_weight(function(r) {
    calls$eval <- calls$eval + 1L
    w <- 1 / as.numeric(r)
    w / sum(w)
  }, name = "no_table_fallback")

  R <- .layer1_rank_matrix(n = 8L, m = 3L, seed = 505L)
  expect_null(RdistNN:::.rank_weight_table(spec, nrow(R)))

  W <- rank_weight_eval_batch(spec, R)
  expect_equal(calls$eval, ncol(R))
  expect_equal(W, .layer1_reference_batch(R, function(r) spec$eval(r)), tolerance = 1e-12)
})

test_that("rank-weight table accessor rejects malformed tables", {
  malformed <- function(table_fun) {
    spec <- list(name = "bad_table", table = table_fun, .__kind = "rank_weight")
    class(spec) <- c("rank_weight", class(spec))
    spec
  }

  expect_error(
    RdistNN:::.rank_weight_table(malformed(function(n) rep(1 / n, n - 1L)), 5L),
    "length 5"
  )
  expect_error(
    RdistNN:::.rank_weight_table(malformed(function(n) rep(Inf, n)), 5L),
    "finite"
  )
  expect_error(
    RdistNN:::.rank_weight_table(malformed(function(n) rep(0, n)), 5L),
    "positive finite sum"
  )
})
