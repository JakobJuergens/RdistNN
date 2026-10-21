test_that("loading fixture helpers preserves the caller RNG state", {
  .test_rng_scope({
    for (file in c("helper-distances.R", "helper-jackknife.R", "helper-user-cate-dnn.R", "helper-user-plqr-dnn.R", "helper-user-np-dnn.R")) {
      set.seed(271L)
      before <- .Random.seed
      sys.source(test_path(file), envir = new.env(parent = environment()))
      expect_identical(.Random.seed, before, info = file)
      rm(".Random.seed", envir = .GlobalEnv)
      sys.source(test_path(file), envir = new.env(parent = environment()))
      expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE), info = file)
    }
  })
})

test_that("shared random fixtures are deterministic without consuming caller draws", {
  .test_rng_scope({
    fixtures <- list(
      query = function() .make_X0(3L, 2L),
      distance_training = function() gen_X(),
      distance_targets = function() gen_X0(),
      distance_target = function() gen_x0(),
      np_interface = function() .make_np_iface_data(),
      tdnn_interface = function() .make_tdnn_iface_data(),
      np = function() .make_np_data(),
      tdnn = function() .make_tdnn_data(),
      cate = function() .make_data(),
      cate_frame = function() .make_df_data(),
      plqr = function() .make_plqr_data(),
      plqr_inference = function() .make_plqr_inference_data()
    )
    for (fixture in fixtures) {
      set.seed(271L)
      before <- .Random.seed
      first <- fixture()
      expect_identical(.Random.seed, before)
      set.seed(991L)
      expect_identical(fixture(), first)
      rm(".Random.seed", envir = .GlobalEnv)
      expect_identical(fixture(), first)
      expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    }
  })
})

test_that("test RNG scopes restore kinds and seed presence after errors", {
  .test_rng_scope({
    for (seeded in c(TRUE, FALSE)) {
      set.seed(271L)
      before <- .Random.seed
      kind <- RNGkind()
      if (!seeded) rm(".Random.seed", envir = .GlobalEnv)
      expect_error(.test_rng_scope({
        RNGkind("L'Ecuyer-CMRG", normal.kind = "Inversion")
        set.seed(99L)
        stop("fixture failure")
      }), "fixture failure")
      expect_identical(RNGkind(), kind)
      if (seeded) expect_identical(.Random.seed, before)
      else expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    }
  })
})
