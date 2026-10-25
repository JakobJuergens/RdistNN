# Pure supplied-signal CATE scores; no estimator weights or nuisance fitting.
rdistnn_bench_cate <- function(n = 2000L, p = 5L, s = 32L, n_subsets = 512L,
                                reps = 5L, seed = 1L, include_cpp_direct = TRUE) {
  rdistnn_bench_require()
  set.seed(seed)
  demo <- cate_demo_data(n, p, seed = seed)
  score <- score_cate_rho(demo$x0, paste0("X", seq_len(p)), s = s,
                          rho_col = "aipw_oracle")
  idx <- rdistnn_bench_make_idx(n, n_subsets, s)
  rdistnn_bench_score(score, demo$data, idx, demo$data$aipw_oracle, reps, include_cpp_direct)
}
