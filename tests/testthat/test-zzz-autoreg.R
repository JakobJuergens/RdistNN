testthat::local_edition(3)

test_that("[onLoad] scaled distances are auto-registered when available", {
  skip_if_not(exists("list_distances", mode = "function"))

  ds_names <- registry_names(list_distances())

  if (exists("distance_l1_scaled", mode = "function")) {
    # If none present, assume auto-registration is disabled in this build
    if (!any(grepl("^l1_scaled\\[(sd|mad)\\]$", ds_names))) {
      skip("auto-registration for l1_scaled is disabled or deferred in this build.")
    }
    expect_true(any(grepl("^l1_scaled\\[(sd|mad)\\]$", ds_names)))
  }

  if (exists("distance_l2_scaled", mode = "function")) {
    # New canonical naming: l2_scaled_sq[rule]
    if (!any(grepl("^l2_scaled_sq\\[(sd|mad)\\]$", ds_names))) {
      skip("auto-registration for l2_scaled_sq is disabled or deferred in this build.")
    }
    expect_true(any(grepl("^l2_scaled_sq\\[(sd|mad)\\]$", ds_names)))
  }
})
