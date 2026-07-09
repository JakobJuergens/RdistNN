# Aliases for internal (unexported) low-level functions that the test suite
# exercises directly. These layers were removed from the public namespace in
# the pre-public API pruning (markdown_notes/rdistnn_public_release_cleanup_note.md,
# section 18); the tests keep covering them through these aliases.

dist_l1 <- RdistNN:::dist_l1
dist_l2 <- RdistNN:::dist_l2
dist_l1_scaled <- RdistNN:::dist_l1_scaled
dist_l2_scaled <- RdistNN:::dist_l2_scaled
order_by_l1 <- RdistNN:::order_by_l1
order_by_l2 <- RdistNN:::order_by_l2
order_by_l1_scaled <- RdistNN:::order_by_l1_scaled
order_by_l2_scaled <- RdistNN:::order_by_l2_scaled
ord_to_ranks <- RdistNN:::ord_to_ranks
tdnn_rank_table_exact <- RdistNN:::tdnn_rank_table_exact
tdnn_rank_table_asymp <- RdistNN:::tdnn_rank_table_asymp
dnn_weights_exact_full_from_ranks <- RdistNN:::dnn_weights_exact_full_from_ranks
dnn_weights_asymp_full_from_ranks <- RdistNN:::dnn_weights_asymp_full_from_ranks
dnn_weights_tdnn_exact_full_from_ranks <- RdistNN:::dnn_weights_tdnn_exact_full_from_ranks
dnn_weights_tdnn_asymp_full_from_ranks <- RdistNN:::dnn_weights_tdnn_asymp_full_from_ranks
