# -------------------------------------------------------------------
# Distance registry basics
# Tags:
#   [registry]Register/list/get
# -------------------------------------------------------------------

test_that("[registry] list_distances includes built-ins and custom registration works", {
  skip_if_not(exists("list_distances", mode = "function"))
  skip_if_not(exists("distance_l1", mode = "function"))
  skip_if_not(exists("distance_l2", mode = "function"))
  skip_if_not(exists("new_distance", mode = "function"))

  if (exists("reset_distance_registry", mode = "function")) reset_distance_registry()

  # Ensure core distances are present (don’t assume .onLoad did it)
  ensure_registered(list(
    l1 = distance_l1(),
    l2 = distance_l2(FALSE) # canonical unsquared "l2"
  ))

  ds_names <- registry_names(list_distances())
  expect_true(all(c("l2", "l1") %in% ds_names))

  # Register a custom distance and retrieve it
  tmp <- new_distance(
    name = "tmp_d",
    eval = function(X, x0, ctx = NULL) rep(0, nrow(X))
  )
  ensure_registered(list(tmp_d = tmp))

  ds_names <- registry_names(list_distances())
  expect_true("tmp_d" %in% ds_names)

  got <- get_distance("tmp_d")
  expect_s3_class(got, "dnn_distance")
})
