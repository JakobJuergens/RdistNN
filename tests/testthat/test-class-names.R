test_that("public specifications use package-owned S3 identities", {
  metric <- new_distance("class_distance", eval = function(X, x0, ctx) abs(X[, 1L] - x0))
  rank <- new_rank_weight("class_rank", eval = function(ranks, ctx) rep(1, length(ranks)))
  score <- as_score(function(idxs, data) sum(data[idxs]), s = 1L)
  tuning <- new_tuning("class_tuning", eval = function(...) list(s = 2L))
  expect_s3_class(metric, "rdistnn_distance")
  expect_s3_class(rank, "rdistnn_rank_weight")
  expect_s3_class(score, "rdistnn_ustat_score")
  expect_s3_class(tuning, "rdistnn_tuning")
  expect_identical(rank$.__kind, "rank_weight")
  restored <- unserialize(serialize(score, NULL))
  expect_identical(score_eval(restored, 1L, c(2, 3)), 2)
})

test_that("current fitted and summary classes dispatch after serialization", {
  X <- matrix(seq_len(12L), ncol = 1L)
  np <- np_dnn(X = X, x0 = 0, Y = seq_len(12L), s = 3L, se = FALSE)
  cate <- cate_dnn(X = X, x0 = 0, pseudo_outcome = seq_len(12L), s = 3L, se = FALSE)
  for (kind in c("np", "cate")) {
    fit <- if (kind == "np") np else cate
    restored <- unserialize(serialize(fit, NULL))
    expect_s3_class(restored, paste0("rdistnn_", kind, "_result"))
    expect_identical(coef(restored), coef(fit))
    expect_s3_class(summary(restored), paste0("rdistnn_summary_", kind, "_result"))
    expect_output(print(summary(restored)), "Conditional|DNN|CATE|summary")
  }
})
