# Aliases for internal (unexported) low-level functions that the test suite
# exercises directly. These aliases provide internal test access without
# making the helpers part of the exported package API.

dist_l1 <- RdistNN:::dist_l1
dist_l2 <- RdistNN:::dist_l2
ord_to_ranks <- RdistNN:::ord_to_ranks
tdnn_rank_table_exact <- RdistNN:::tdnn_rank_table_exact
tdnn_rank_table_asymp <- RdistNN:::tdnn_rank_table_asymp
dnn_weights_exact_full_from_ranks <- RdistNN:::dnn_weights_exact_full_from_ranks
dnn_weights_asymp_full_from_ranks <- RdistNN:::dnn_weights_asymp_full_from_ranks
