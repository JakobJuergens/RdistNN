test_that("internal validators split NA vs non-finite diagnostics", {
  v_vec_n <- getFromNamespace(".validate_vector_nrow", "RdistNN")
  v_mat <- getFromNamespace(".validate_matrix", "RdistNN")
  v_x0 <- getFromNamespace(".validate_target_vector", "RdistNN")
  v_x0_inv_s <- getFromNamespace(".validate_x0_inv_s", "RdistNN")
  v_bin <- getFromNamespace(".validate_binary01", "RdistNN")
  v_subs <- getFromNamespace(".validate_subset_indices", "RdistNN")

  expect_error(
    v_vec_n(c(1, NA_real_), n = 2L, arg = "y", finite = TRUE),
    "y.*NA",
    ignore.case = TRUE
  )
  expect_error(
    v_vec_n(c(1, Inf), n = 2L, arg = "y", finite = TRUE),
    "y.*non-finite",
    ignore.case = TRUE
  )

  expect_error(
    v_mat(matrix(c(1, NA_real_), ncol = 1L), arg = "X", finite = TRUE),
    "X.*NA",
    ignore.case = TRUE
  )
  expect_error(
    v_mat(matrix(c(1, NaN), ncol = 1L), arg = "X", finite = TRUE),
    "X.*non-finite",
    ignore.case = TRUE
  )

  expect_error(
    v_x0(c(0, NA_real_), p = 2L, x0_arg = "x0", finite = TRUE),
    "x0.*NA",
    ignore.case = TRUE
  )
  expect_error(
    v_x0(c(0, -Inf), p = 2L, x0_arg = "x0", finite = TRUE),
    "x0.*non-finite",
    ignore.case = TRUE
  )

  expect_error(
    v_x0_inv_s(c(0, NA_real_), c(1, 1), p = 2L, finite = TRUE),
    "x0/inv_s.*NA",
    ignore.case = TRUE
  )
  expect_error(
    v_x0_inv_s(c(0, 1), c(1, Inf), p = 2L, finite = TRUE),
    "x0/inv_s.*non-finite",
    ignore.case = TRUE
  )

  expect_error(v_bin(c(0, 1, NA_real_), arg = "t_col"), "t_col.*NA", ignore.case = TRUE)
  expect_error(v_bin(c(0, 1, Inf), arg = "t_col"), "t_col.*non-finite", ignore.case = TRUE)
  expect_error(v_bin(c(0, 1, 2), arg = "t_col"), "binary", ignore.case = TRUE)

  expect_error(v_subs(c(1, NA_real_), n = 5L, s = 2L), "idxs.*NA", ignore.case = TRUE)
  expect_error(v_subs(c(1, Inf), n = 5L, s = 2L), "idxs.*non-finite", ignore.case = TRUE)
  expect_error(v_subs(c(1, 2.5), n = 5L, s = 2L), "integer indices", ignore.case = TRUE)

  expect_error(
    v_subs(matrix(c(1, NA_real_), nrow = 1L), n = 5L, s = 2L),
    "subs.*NA",
    ignore.case = TRUE
  )
  expect_error(
    v_subs(matrix(c(1, NaN), nrow = 1L), n = 5L, s = 2L),
    "subs.*non-finite",
    ignore.case = TRUE
  )
})
