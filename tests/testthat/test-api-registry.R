.registry_name_contract_cases <- function() {
  registries <- list(
    distance = list(register = register_distance, get = get_distance, names = distance_names,
      env = RdistNN:::.dist_registry, spec = distance_l2()),
    rank = list(register = register_rank_weight, get = get_rank_weight, names = rank_weight_names,
      env = RdistNN:::.rank_registry,
      spec = new_rank_weight("registry_rank", eval = function(ranks, ctx) rep(1, length(ranks)))),
    score = list(register = register_score, get = get_score, names = list_scores,
      env = RdistNN:::.score_registry,
      spec = as_score(function(idxs, data) sum(data[idxs]), s = 1L, name = "registry_score")),
    tuning = list(register = register_tuning, get = get_tuning, names = tuning_names,
      env = RdistNN:::.tuning_registry,
      spec = new_tuning("registry_tuning", eval = function(...) list(s = 2L)))
  )
  snapshots <- lapply(registries, function(registry) as.list.environment(registry$env, all.names = TRUE))
  on.exit(for (i in seq_along(registries)) {
    env <- registries[[i]]$env
    rm(list = ls(env, all.names = TRUE), envir = env)
    list2env(snapshots[[i]], env)
  }, add = TRUE)
  for (registry in registries) {
    for (bad in list(NA_character_, "", character(), NULL, c("one", "two"), 1, TRUE)) {
      expect_error(registry$register(bad, registry$spec), "name.*non-empty character scalar")
      expect_error(registry$get(bad), "name.*non-empty character scalar")
    }
    for (name in c(".audit registry / name", " ")) {
      registry$register(name, registry$spec)
      expect_identical(registry$get(name)$name, registry$spec$name)
      expect_true(exists(name, registry$env, inherits = FALSE))
      expect_true(name %in% registry$names())
      expect_identical(anyDuplicated(class(registry$get(name))), 0L)
    }
    missing <- tryCatch(registry$get("__audit_unregistered_name__"), error = identity)
    expect_s3_class(missing, "rdistnn_input_error")
    expect_identical(missing$code, "argument_invalid")
    expect_identical(missing$argument, "name")
    expect_match(missing$actual, "__audit_unregistered_name__", fixed = TRUE)
  }
}

test_that("registry boundaries reject missing and empty names without restricting valid keys", {
  .registry_name_contract_cases()
})
