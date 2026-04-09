test_that("built-in paper tuning is registered and classed", {
  expect_true("paper_undersmooth" %in% tuning_names())
  expect_true(is_tuning(get_tuning("paper_undersmooth")))
})

test_that("new_tuning()/as_tuning()/is_tuning() support eval and eval_grid specs", {
  tun_eval <- new_tuning(
    name = "fixed_7",
    eval = function(n, kx, k, se, jackknife, ctx = NULL) list(s = 7L)
  )
  expect_true(is_tuning(tun_eval))

  tun_grid_only <- new_tuning(
    name = "grid_only",
    eval_grid = function(n, kx, k, se, jackknife, ctx = NULL) {
      data.frame(kappa = c(0.8, 1.1), s = c(5L, 7L), stringsAsFactors = FALSE)
    }
  )
  expect_true(is_tuning(tun_grid_only))

  tun_as <- as_tuning(function(...) 9L, name = "fixed_9")
  expect_true(is_tuning(tun_as))

  out_eval <- tuning_eval(tun_eval, n = 200L, kx = 2L, k = 5L, mode = "default")
  expect_equal(out_eval$s, 7L)

  out_grid_default <- tuning_eval(tun_grid_only, n = 200L, kx = 2L, k = 5L, mode = "default")
  expect_equal(out_grid_default$s, 5L)

  out_grid <- tuning_eval(tun_grid_only, n = 200L, kx = 2L, k = 5L, mode = "grid")
  expect_equal(out_grid$candidates$s, c(5L, 7L))
})

test_that("register_tuning()/get_tuning()/tuning_names() support function and list specs", {
  nm_fun <- paste0("tmp_fixed_", sample.int(1e6, 1L))
  register_tuning(nm_fun, function(...) 11L)
  expect_true(nm_fun %in% tuning_names())
  expect_equal(tuning_eval(get_tuning(nm_fun), n = 120L, kx = 2L, k = 5L)$s, 11L)

  nm_list <- paste0("tmp_list_", sample.int(1e6, 1L))
  register_tuning(nm_list, list(eval = function(...) list(s = 13L)))
  got <- get_tuning(nm_list)
  expect_true(is_tuning(got))
  expect_equal(got$name, nm_list)
  expect_equal(tuning_eval(got, n = 100L, kx = 2L, k = 5L)$s, 13L)

  nm_a <- paste0("a_", sample.int(1e6, 1L))
  nm_z <- paste0("z_", sample.int(1e6, 1L))
  register_tuning(nm_z, function(...) 4L)
  register_tuning(nm_a, function(...) 4L)
  tnames <- tuning_names()
  expect_true(all(sort(tnames) == tnames))
})

test_that("registry and constructor validation failures are explicit", {
  expect_error(register_tuning("", function(...) 1L), "non-empty character scalar")
  expect_error(new_tuning(name = "bad_new"), "at least one of `eval` or `eval_grid`")
  expect_error(as_tuning(1), "`fun` must be a function")

  expect_error(get_tuning(1), "character scalar")
  expect_error(get_tuning("__definitely_missing__"), "Unknown tuning")
  msg <- tryCatch(
    get_tuning("__definitely_missing__"),
    error = function(e) conditionMessage(e)
  )
  expect_match(msg, "Expected:", fixed = TRUE)
  expect_match(msg, "Got:", fixed = TRUE)
})

test_that("tuning_eval validates argument types and accepts bare function tuning", {
  expect_error(
    tuning_eval(tuning_dnn_paper(), n = 10.5, kx = 2L, k = 2L),
    "integer scalar"
  )
  expect_error(
    tuning_eval(tuning_dnn_paper(), n = 100L, kx = 2L, k = 2L, se = NA),
    "TRUE or FALSE"
  )
  expect_error(
    tuning_eval(tuning_dnn_paper(), n = 100L, kx = 2L, k = 2L, jackknife = 1),
    "character scalar"
  )
  expect_error(
    tuning_eval(tuning_dnn_paper(), n = 100L, kx = 2L, k = 2L, jackknife = "bad"),
    "should be one of"
  )
  expect_error(
    tuning_eval(list(), n = 100L, kx = 2L, k = 2L),
    "must be a tuning spec"
  )

  out_fun <- tuning_eval(function(n, kx) floor(n^((kx + 1) / (kx + 2))), n = 100L, kx = 2L, k = 3L)
  expect_equal(out_fun$s, 31L)
})

test_that("default and grid output validators reject malformed user outputs", {
  tun_missing_s <- new_tuning(name = "bad_missing_s", eval = function(...) list(meta = list(a = 1)))
  expect_error(
    tuning_eval(tun_missing_s, n = 50L, kx = 2L, k = 2L, mode = "default"),
    "must return a list with field `s` or a numeric scalar"
  )

  tun_bad_s <- new_tuning(name = "bad_s", eval = function(...) list(s = 2.5))
  expect_error(
    tuning_eval(tun_bad_s, n = 50L, kx = 2L, k = 2L, mode = "default"),
    "integer scalar"
  )

  tun_grid_not_df <- new_tuning(name = "bad_grid_not_df", eval_grid = function(...) list(s = 2L))
  expect_error(
    tuning_eval(tun_grid_not_df, n = 50L, kx = 2L, k = 2L, mode = "grid"),
    "must return a data.frame"
  )

  tun_grid_missing_s <- new_tuning(
    name = "bad_grid_missing_s",
    eval_grid = function(...) data.frame(kappa = 1, stringsAsFactors = FALSE)
  )
  expect_error(
    tuning_eval(tun_grid_missing_s, n = 50L, kx = 2L, k = 2L, mode = "grid"),
    "must return a data.frame with column `s`"
  )

  tun_grid_bad_s <- new_tuning(
    name = "bad_grid_s",
    eval_grid = function(...) data.frame(kappa = 1, s = 2.5, stringsAsFactors = FALSE)
  )
  expect_error(
    tuning_eval(tun_grid_bad_s, n = 50L, kx = 2L, k = 2L, mode = "grid"),
    "column `s` must be finite integer-like numeric values"
  )
})

test_that("feasibility flags and grid reasons are deterministic and correctly prioritized", {
  tun <- new_tuning(
    name = "feasibility_probe",
    eval = function(...) list(s = 2L),
    eval_grid = function(...) {
      data.frame(kappa = c(0.5, 1.0, 1.5), s = c(1L, 5L, 3L), stringsAsFactors = FALSE)
    }
  )

  out_default <- tuning_eval(
    tun,
    n = 5L, kx = 2L, k = 3L,
    se = TRUE, jackknife = "fast", mode = "default"
  )
  expect_false(out_default$feasibility$k_le_s)
  expect_true(out_default$feasibility$s_lt_n)
  expect_equal(out_default$meta$name, "feasibility_probe")
  expect_equal(out_default$meta$n, 5L)
  expect_equal(out_default$meta$kx, 2L)

  out_grid <- tuning_eval(
    tun,
    n = 5L, kx = 2L, k = 3L,
    se = TRUE, jackknife = "fast", mode = "grid"
  )
  expect_equal(out_grid$candidates$reason, c("k > s", "jackknife requires s < n", "ok"))
  expect_equal(out_grid$candidates$feasible, c(FALSE, FALSE, TRUE))
})

test_that("eval-only and grid-only tunings use documented fallback behavior", {
  tun_eval_only <- new_tuning(name = "eval_only", eval = function(...) list(s = 8L))
  out_grid <- tuning_eval(tun_eval_only, n = 100L, kx = 2L, k = 3L, mode = "grid")
  expect_equal(nrow(out_grid$candidates), 1L)
  expect_true(is.na(out_grid$candidates$kappa[1L]))
  expect_equal(out_grid$candidates$s, 8L)

  tun_grid_empty <- new_tuning(
    name = "grid_empty",
    eval_grid = function(...) data.frame(s = integer(0L), stringsAsFactors = FALSE)
  )
  expect_error(
    tuning_eval(tun_grid_empty, n = 100L, kx = 2L, k = 3L, mode = "default"),
    "at least one candidate row"
  )
})

test_that("paper_undersmooth default and grid are deterministic", {
  out_default <- tuning_eval("paper_undersmooth", n = 1000L, kx = 2L, k = 5L, mode = "default")
  expect_equal(out_default$s, 177L)
  expect_equal(out_default$meta$kappa, 1)
  expect_equal(out_default$meta$rounding, "floor")
  expect_equal(out_default$meta$formula, "kappa * n^((kx + 1)/(kx + 2))")

  out_grid <- tuning_eval("paper_undersmooth", n = 1000L, kx = 2L, k = 5L, mode = "grid")
  expect_s3_class(out_grid$candidates, "data.frame")
  expect_equal(out_grid$candidates$kappa, c(0.75, 1, 1.25, 1.5))
  expect_equal(out_grid$candidates$s, c(133L, 177L, 222L, 266L))
  expect_true(all(out_grid$candidates$feasible))
})

test_that("paper_undersmooth supports rounding choices and kx_override", {
  expect_equal(
    tuning_eval(tuning_dnn_paper(rounding = "floor"), n = 50L, kx = 1L, k = 2L)$s,
    13L
  )
  expect_equal(
    tuning_eval(tuning_dnn_paper(rounding = "ceiling"), n = 50L, kx = 1L, k = 2L)$s,
    14L
  )
  expect_equal(
    tuning_eval(tuning_dnn_paper(rounding = "round"), n = 50L, kx = 1L, k = 2L)$s,
    14L
  )

  tun_kx_override <- tuning_dnn_paper(kx_override = 5L)
  out_override <- tuning_eval(tun_kx_override, n = 100L, kx = 2L, k = 2L)
  expect_equal(out_override$s, 51L)
})

test_that("paper_undersmooth input validation is explicit", {
  expect_error(tuning_dnn_paper(kappa = 0), "positive finite numeric scalar")
  expect_error(tuning_dnn_paper(kappa_grid = c(0.5, 0)), "positive finite numeric vector")
  expect_error(tuning_dnn_paper(kx_override = 2.5), "integer scalar")
})
