# -------------------------------------------------------------------
# Distance registry basics
# Tags:
#   [registry]Register/list/get
# -------------------------------------------------------------------

test_that("[registry] list_distances includes built-ins and custom registration works", {
  skip_if_not(exists("list_distances", mode = "function"))
  ds <- list_distances()
  expect_true(all(c("l2", "l1") %in% ds))

  # Register a custom distance under a temp name and retrieve it
  tmp <- new_distance(name = "tmp_d", eval = function(X, x0, ctx = NULL) rep(0, nrow(X)))
  expect_silent(register_distance(tmp))
  expect_true("tmp_d" %in% list_distances())
  got <- get_distance("tmp_d")
  expect_s3_class(got, "dnn_distance")
})

# tests/testthat/test-zzz-autoreg.R
testthat::local_edition(3)

test_that("[onLoad] scaled distances are auto-registered when available", {
  skip_if_not(exists("list_distances", mode = "function"))
  ds <- list_distances()

  if (exists("distance_l1_scaled", mode = "function")) {
    expect_true(any(grepl("^l1_scaled\\[(sd|mad)\\]$", ds)))
  }
  if (exists("distance_l2_scaled", mode = "function")) {
    expect_true(any(grepl("^l2_scaled\\[(sd|mad)\\]_sq$", ds)))
  }
})
