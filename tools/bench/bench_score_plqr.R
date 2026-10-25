# Supplied-coefficient PLQR moment scores; no coefficient fitting or inference.
rdistnn_bench_plqr <- function(n = 2000L, p = 5L, p_d = 2L, s = 32L,
                                tau = 0.5, j = 1L, n_subsets = 512L,
                                reps = 5L, seed = 1L, include_cpp_direct = TRUE) {
  rdistnn_bench_require()
  set.seed(seed)
  demo <- plqr_demo_data(n, p, p_d, tau = tau, seed = seed)
  d <- demo$data
  D_cols <- paste0("D", seq_len(p_d))
  mD_cols <- paste0("mD", seq_len(p_d))
  theta <- demo$truth$theta_true
  score <- score_plqr_component(demo$x0, paste0("X", seq_len(p)), D_cols,
    g_col = "g_oracle", mD_cols = mD_cols, y_col = "Y", theta = theta,
    j = j, tau = tau, s = s)
  threshold <- as.vector(as.matrix(d[D_cols]) %*% theta) + d$g_oracle
  signal <- (tau - as.numeric(d$Y <= threshold)) * (d[[D_cols[j]]] - d[[mD_cols[j]]])
  idx <- rdistnn_bench_make_idx(n, n_subsets, s)
  rdistnn_bench_score(score, d, idx, signal, reps, include_cpp_direct)
}
