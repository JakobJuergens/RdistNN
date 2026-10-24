.expect_candidate_schema <- function(candidates) {
  expect_s3_class(candidates, "data.frame")
  expect_true(all(c(
    "candidate_id", "computed_key", "is_default", "feasible",
    "reason", "status", "lambda", "s"
  ) %in% names(candidates)))
  expect_type(candidates$candidate_id, "character")
  expect_type(candidates$computed_key, "character")
  expect_type(candidates$is_default, "logical")
  expect_type(candidates$status, "character")
  expect_equal(sum(candidates$is_default), 1L)
  expect_equal(candidates$status, ifelse(candidates$feasible, "ok", "infeasible"))
  expect_equal(anyDuplicated(candidates$candidate_id), 0L)
}

test_that("manual scalar s gets the common one-row candidate schema", {
  res <- RdistNN:::.dnn_resolve_s(
    s = 5L,
    tuning = NULL,
    n = 40L,
    p = 2L,
    k = 1L,
    se = FALSE,
    jackknife = "none"
  )

  candidates <- res$tuning$grid$candidates
  .expect_candidate_schema(candidates)
  expect_equal(res$s, 5L)
  expect_false(res$tuning_used)
  expect_equal(candidates$candidate_id, "s_5")
  expect_equal(candidates$computed_key, "s_5")
  expect_true(candidates$is_default[1L])
  expect_true(candidates$feasible[1L])
  expect_equal(candidates$reason, "ok")
  expect_equal(candidates$status, "ok")
  expect_true(is.na(candidates$lambda[1L]))
})

test_that("paper tuning formats candidate ids and identifies the default row", {
  tun <- tuning_dnn_paper(lambda = 1, lambda_grid = c(0.75, 1, 1.25), rounding = "floor")
  res <- RdistNN:::.dnn_resolve_s(
    s = NULL,
    tuning = tun,
    n = 1000L,
    p = 2L,
    k = 5L,
    se = TRUE,
    jackknife = "fast"
  )

  candidates <- res$tuning$grid$candidates
  .expect_candidate_schema(candidates)
  expect_equal(res$s, 177L)
  expect_equal(candidates$candidate_id, c("lambda_0p75", "lambda_1", "lambda_1p25"))
  expect_equal(candidates$computed_key, c("s_133", "s_177", "s_222"))
  expect_equal(candidates$is_default, c(FALSE, TRUE, FALSE))
  expect_true(all(candidates$feasible))
  expect_equal(candidates$status, rep("ok", 3L))
})

test_that("duplicate candidate ids are suffixed and duplicate s shares computed keys", {
  tun <- new_tuning(
    name = "duplicate_probe",
    eval = function(...) list(s = 4L, meta = list(lambda = 1)),
    eval_grid = function(...) {
      data.frame(
        lambda = c(1, 1, 1.25, NA, NA),
        s = c(4L, 4L, 4L, 7L, 7L),
        stringsAsFactors = FALSE
      )
    }
  )
  res <- RdistNN:::.dnn_resolve_s(
    s = NULL,
    tuning = tun,
    n = 50L,
    p = 2L,
    k = 2L,
    se = FALSE,
    jackknife = "none"
  )

  candidates <- res$tuning$grid$candidates
  .expect_candidate_schema(candidates)
  expect_equal(
    candidates$candidate_id,
    c("lambda_1", "lambda_1__2", "lambda_1p25", "s_7", "s_7__2")
  )
  expect_equal(candidates$computed_key, c("s_4", "s_4", "s_4", "s_7", "s_7"))
  expect_equal(candidates$is_default, c(TRUE, FALSE, FALSE, FALSE, FALSE))
})

test_that("infeasible non-default candidates are recorded in status", {
  tun <- new_tuning(
    name = "infeasible_probe",
    eval = function(...) list(s = 3L, meta = list(lambda = 1)),
    eval_grid = function(...) {
      data.frame(
        lambda = c(0.5, 1, 2),
        s = c(1L, 3L, 40L),
        stringsAsFactors = FALSE
      )
    }
  )
  res <- RdistNN:::.dnn_resolve_s(
    s = NULL,
    tuning = tun,
    n = 40L,
    p = 2L,
    k = 2L,
    se = TRUE,
    jackknife = "fast"
  )

  candidates <- res$tuning$grid$candidates
  .expect_candidate_schema(candidates)
  expect_equal(candidates$reason, c("k > s", "ok", "jackknife requires s < n"))
  expect_equal(candidates$feasible, c(FALSE, TRUE, FALSE))
  expect_equal(candidates$status, c("infeasible", "ok", "infeasible"))
  expect_equal(candidates$is_default, c(FALSE, TRUE, FALSE))
})

test_that("custom grids get a synthetic default row when the default s is absent", {
  tun <- new_tuning(
    name = "synthetic_default_probe",
    eval = function(...) list(s = 6L, meta = list(lambda = 2)),
    eval_grid = function(...) {
      data.frame(
        lambda = c(0.5, 1),
        s = c(3L, 4L),
        label = c("low", "mid"),
        stringsAsFactors = FALSE
      )
    }
  )
  res <- RdistNN:::.dnn_resolve_s(
    s = NULL,
    tuning = tun,
    n = 50L,
    p = 2L,
    k = 2L,
    se = FALSE,
    jackknife = "none"
  )

  candidates <- res$tuning$grid$candidates
  .expect_candidate_schema(candidates)
  expect_equal(nrow(candidates), 3L)
  expect_equal(candidates$s, c(6L, 3L, 4L))
  expect_equal(candidates$candidate_id, c("s_6", "lambda_0p5", "lambda_1"))
  expect_equal(candidates$is_default, c(TRUE, FALSE, FALSE))
  expect_true(is.na(candidates$lambda[1L]))
  expect_true(is.na(candidates$label[1L]))
})

test_that("wrappers still evaluate only the default candidate", {
  withr::local_preserve_seed()
  n_calls <- new.env(parent = emptyenv())
  n_calls$count <- 0L
  dist <- new_distance(
    name = "counted_l2",
    eval_batch = function(X, X0, ctx = NULL) {
      n_calls$count <- n_calls$count + 1L
      vapply(seq_len(nrow(X0)), function(j) {
        rowSums((X - matrix(X0[j, ], nrow(X), ncol(X), byrow = TRUE))^2)
      }, numeric(nrow(X)))
    }
  )
  tun <- new_tuning(
    name = "multi_candidate_no_eval",
    eval = function(...) list(s = 4L, meta = list(lambda = 1)),
    eval_grid = function(...) {
      data.frame(lambda = c(0.75, 1, 1.25), s = c(3L, 4L, 5L), stringsAsFactors = FALSE)
    }
  )
  set.seed(3107)
  X <- matrix(rnorm(48), nrow = 24L, ncol = 2L)
  Y <- rnorm(24L)

  fit <- np_dnn(
    X = X,
    Y = Y,
    x0 = c(0, 0),
    distance = dist,
    method = "exact",
    k = 2L,
    tuning = tun,
    se = FALSE,
    return = "tuning"
  )

  expect_equal(n_calls$count, 1L)
  expect_equal(fit$meta$s, 4L)
  expect_equal(nrow(fit$tuning$grid$candidates), 3L)
  expect_equal(fit$tuning$grid$candidates$is_default, c(FALSE, TRUE, FALSE))
})
