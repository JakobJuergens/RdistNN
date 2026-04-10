skip_if_not(
  all(vapply(
    c("np_dnn", "diagnostics"),
    exists, logical(1),
    mode = "function"
  )),
  "np_dnn() or diagnostics() not found; skipping file"
)

set.seed(777)

.diag_top_share <- function(w, k) {
  k_use <- max(1L, min(length(w), as.integer(k)))
  sum(sort(as.numeric(w), decreasing = TRUE)[seq_len(k_use)])
}

.diag_weighted_median <- function(x, w) {
  ord <- order(x)
  x_ord <- x[ord]
  w_ord <- w[ord]
  cw <- cumsum(w_ord) / sum(w_ord)
  x_ord[which(cw >= 0.5)[1L]]
}

test_that("diagnostics(np_result) returns expected schema and identities", {
  n <- 45L
  p <- 2L
  m <- 3L
  x_mat <- matrix(rnorm(n * p), n, p)
  x0_mat <- matrix(rnorm(m * p), m, p)
  y_vec <- x_mat[, 1L] - 0.4 * x_mat[, 2L] + rnorm(n, sd = 0.2)

  fit <- np_dnn(
    X = x_mat, X0 = x0_mat, Y = y_vec,
    distance = "l2", method = "exact",
    k = 2L, s = 4L, se = FALSE,
    return = c("W", "dx")
  )
  d <- diagnostics(fit)

  expect_s3_class(d, "dnn_diagnostics")
  expect_equal(nrow(d), m)
  expect_true(all(c(
    "target_id", "n_eff", "max_weight",
    "top1_weight_share", "top5_weight_share", "top10_weight_share",
    "radius_rank_k", "radius_rank_s",
    "weighted_mean_distance", "weighted_median_distance"
  ) %in% names(d)))

  expect_equal(d$n_eff, 1 / colSums(fit$W^2), tolerance = 1e-12)
  expect_equal(d$max_weight, apply(fit$W, 2L, max), tolerance = 1e-12)
  expect_equal(
    d$top1_weight_share,
    apply(fit$W, 2L, .diag_top_share, k = 1L),
    tolerance = 1e-12
  )
  expect_equal(
    d$top5_weight_share,
    apply(fit$W, 2L, .diag_top_share, k = 5L),
    tolerance = 1e-12
  )
  expect_equal(
    d$top10_weight_share,
    apply(fit$W, 2L, .diag_top_share, k = 10L),
    tolerance = 1e-12
  )

  k_use <- fit$meta$k
  s_use <- fit$meta$s
  expect_equal(
    d$radius_rank_k,
    apply(fit$dx, 2L, function(col) sort(col, partial = k_use)[k_use]),
    tolerance = 1e-12
  )
  expect_equal(
    d$radius_rank_s,
    apply(fit$dx, 2L, function(col) sort(col, partial = s_use)[s_use]),
    tolerance = 1e-12
  )
  expect_equal(
    d$weighted_mean_distance,
    colSums(fit$dx * fit$W) / colSums(fit$W),
    tolerance = 1e-12
  )
  expect_equal(
    d$weighted_median_distance,
    vapply(seq_len(ncol(fit$W)), function(j) .diag_weighted_median(fit$dx[, j], fit$W[, j]), numeric(1L)),
    tolerance = 1e-12
  )
})

test_that("diagnostics(np_result) reports unavailable distance fields when dx is absent", {
  n <- 40L
  x_mat <- matrix(rnorm(n * 2L), n, 2L)
  y_vec <- x_mat[, 1L] + rnorm(n, sd = 0.2)

  fit <- np_dnn(
    X = x_mat, x0 = c(0, 0), Y = y_vec,
    method = "exact", k = 2L, s = 3L, se = FALSE,
    return = "W"
  )
  d <- diagnostics(fit)

  expect_true(all(!d$distance_diagnostics_available))
  expect_true(all(is.na(d$radius_rank_k)))
  expect_true(all(is.na(d$radius_rank_s)))
  expect_true(all(is.na(d$weighted_mean_distance)))
  expect_true(all(is.na(d$weighted_median_distance)))
})

test_that("diagnostics(np_result) errors clearly when W is not retained", {
  n <- 35L
  x_mat <- matrix(rnorm(n * 2L), n, 2L)
  y_vec <- x_mat[, 1L] + rnorm(n, sd = 0.2)

  fit <- np_dnn(
    X = x_mat, x0 = c(0, 0), Y = y_vec,
    method = "exact", k = 2L, s = 3L, se = FALSE,
    return = "ess"
  )

  expect_error(
    diagnostics(fit),
    "refit with `return` including \"W\""
  )
})

test_that("diagnostics(np_result) includes tuning metadata summary", {
  demo <- np_demo_data(n = 80L, p = 2L, seed = 19L)
  x_mat <- as.matrix(demo$data[, c("X1", "X2")])
  y_vec <- demo$data$Y
  tun <- tuning_dnn_paper(kappa = 1.25)

  fit <- np_dnn(
    X = x_mat, x0 = demo$x0, Y = y_vec,
    method = "exact", k = 2L, tuning = tun, se = FALSE,
    return = "W"
  )
  d <- diagnostics(fit)

  expect_true(all(d$tuning_source == "tuning"))
  expect_true(all(d$tuning_name == "paper_undersmooth"))
  expect_true(all(d$tuning_grid_available))
  expect_true(all(d$tuning_default_s == fit$meta$s))
})

test_that("print.dnn_diagnostics runs", {
  n <- 30L
  x_mat <- matrix(rnorm(n * 2L), n, 2L)
  y_vec <- x_mat[, 1L] + rnorm(n, sd = 0.2)
  fit <- np_dnn(
    X = x_mat, x0 = c(0, 0), Y = y_vec,
    method = "exact", k = 2L, s = 3L, se = FALSE, return = "W"
  )
  d <- diagnostics(fit)
  txt <- capture.output(print(d))
  expect_true(any(grepl("<dnn_diagnostics>", txt, fixed = TRUE)))
})
