.diag_grid_tuning <- function(default_s = 4L, s_grid = c(2L, 4L, 6L)) {
  new_tuning(
    name = "diag_grid",
    eval = function(...) list(s = default_s, meta = list(lambda = 1)),
    eval_grid = function(...) {
      data.frame(
        lambda = c(0.5, 1, 1.5),
        s = s_grid,
        stringsAsFactors = FALSE
      )
    }
  )
}

test_that("diagnostics(np_result) includes tuning metadata summary", {
  .test_rng_scope({
    demo <- np_demo_data(n = 80L, p = 2L, seed = 19L)
    x_mat <- as.matrix(demo$data[, c("X1", "X2")])
    y_vec <- demo$data$Y
    tun <- tuning_dnn_paper(lambda = 1.25)

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
    expect_true(all(d$tuning_grid_n_candidates == nrow(attr(d, "tuning_grid"))))
    expect_true(all(d$tuning_lambda_default == 1.25))
  }, seed = 777L)
})

test_that("candidate grid diagnostics are available without retained default W", {
  .test_rng_scope({
    n <- 36L
    p <- 2L
    m <- 3L
    x_mat <- matrix(rnorm(n * p), n, p)
    x0_mat <- matrix(rnorm(m * p), m, p)
    y_vec <- x_mat[, 1L] - 0.25 * x_mat[, 2L] + rnorm(n, sd = 0.2)

    fit <- suppressWarnings(np_dnn(
      X = x_mat, X0 = x0_mat, Y = y_vec,
      method = "exact", k = 2L,
      tuning = .diag_grid_tuning(default_s = 4L, s_grid = c(1L, 4L, 6L)),
      tuning_mode = "grid",
      se = TRUE, jackknife = "fast",
      return = character(0)
    ))

    expect_false("W" %in% names(fit))
    expect_error(
      diagnostics(fit),
      "refit with `return` including \"W\""
    )

    d_all <- diagnostics(fit, candidates = "all")
    expect_s3_class(d_all, "dnn_diagnostics")
    expect_equal(d_all$candidate_id, fit$tuning_grid$candidates$candidate_id)
    expect_equal(d_all$status, fit$tuning_grid$candidates$status)
    expect_true(all(c("candidate_id", "n_eff", "max_weight", "top5_weight_share") %in% names(d_all)))
    expect_true(all(is.na(d_all$n_eff[d_all$status != "ok"])))
    expect_true(all(is.finite(d_all$n_eff[d_all$status == "ok"])))
  }, seed = 777L)
})

test_that("candidate grid diagnostics validate selections and stale grid rows", {
  .test_rng_scope({
    n <- 36L
    p <- 2L
    m <- 3L
    x_mat <- matrix(rnorm(n * p), n, p)
    x0_mat <- matrix(rnorm(m * p), m, p)
    y_vec <- x_mat[, 1L] - 0.25 * x_mat[, 2L] + rnorm(n, sd = 0.2)

    fit <- suppressWarnings(np_dnn(
      X = x_mat, X0 = x0_mat, Y = y_vec,
      method = "exact", k = 2L,
      tuning = .diag_grid_tuning(default_s = 4L, s_grid = c(1L, 4L, 6L)),
      tuning_mode = "grid",
      se = TRUE, jackknife = "fast",
      return = character(0)
    ))
    ok_id <- fit$tuning_grid$candidates$candidate_id[
      fit$tuning_grid$candidates$status == "ok"
    ][[1L]]

    expect_error(
      diagnostics(fit, candidates = character(0)),
      "must be \"default\", \"all\", or a character vector"
    )
    expect_error(
      diagnostics(fit, candidates = 1L),
      "must be \"default\", \"all\", or a character vector"
    )
    expect_error(
      diagnostics(fit, candidates = "not_a_candidate"),
      "Unknown candidate id"
    )

    missing_row <- fit
    missing_row$tuning_grid$diagnostics <- missing_row$tuning_grid$diagnostics[
      missing_row$tuning_grid$diagnostics$candidate_id != ok_id,
      ,
      drop = FALSE
    ]
    expect_error(
      diagnostics(missing_row, candidates = ok_id),
      "missing row\\(s\\) for selected candidate"
    )

    duplicate_row <- fit
    duplicate_row$tuning_grid$diagnostics <- rbind(
      duplicate_row$tuning_grid$diagnostics,
      duplicate_row$tuning_grid$diagnostics[
        duplicate_row$tuning_grid$diagnostics$candidate_id == ok_id,
        ,
        drop = FALSE
      ]
    )
    expect_error(
      diagnostics(duplicate_row, candidates = ok_id),
      "duplicate candidate row"
    )

    missing_top_share <- fit
    missing_top_share$tuning_grid$diagnostics$top5_weight_share <- NULL
    expect_error(
      diagnostics(missing_top_share, candidates = ok_id),
      "top5_weight_share"
    )
  }, seed = 777L)
})
