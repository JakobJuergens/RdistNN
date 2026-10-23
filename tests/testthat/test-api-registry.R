.registry_name_contract_cases <- function() {
  registries <- list(
    distance = list(register = register_distance, get = get_distance,
      env = RdistNN:::.dist_registry, spec = distance_l2()),
    rank = list(register = register_rank_weight, get = get_rank_weight,
      env = RdistNN:::.rank_registry,
      spec = new_rank_weight("registry_rank", eval = function(ranks, ctx) rep(1, length(ranks)))),
    score = list(register = register_score, get = get_score,
      env = RdistNN:::.score_registry,
      spec = as_score(function(idxs, data) sum(data[idxs]), s = 1L, name = "registry_score")),
    tuning = list(register = register_tuning, get = get_tuning,
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
    }
    expect_error(registry$get("__audit_unregistered_name__"), "Unknown|No score registered")
  }
}

test_that("registry boundaries reject missing and empty names without restricting valid keys", {
  .registry_name_contract_cases()
})
