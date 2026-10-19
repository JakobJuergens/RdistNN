# Adapter for the target-by-neighbor finite-support API. Use the same retained
# deletion primitive as arbitrary rank tables, then aggregate targets stably.
.jackknife_finite_accumulation_R <- function(ranks, a, values, n_ref, normalize) {
  targets <- nrow(ranks)
  s <- length(a)
  q <- dim(values)[3L]
  mass <- .compensated_sum_R(a)
  denominator <- if (normalize) mass else 1
  baseline <- matrix(0, targets, q)
  deletions <- array(0, c(targets, s + 1L, q))
  active <- .jackknife_active_ranks(a)
  for (t in seq_len(targets)) {
    for (component in seq_len(q)) {
      original <- values[t, , component]
      baseline[t, component] <- .compensated_sum_R(a * original[seq_len(s)]) / denominator
      level <- .numeric_signal_level(matrix(original[active], ncol = 1L))
      z <- numeric(s + 1L)
      z[active] <- original[active] - level
      deletions[t, , component] <- .compensated_deletions_R(z, a, level * mass) / denominator
    }
  }
  # Neighbor rows are in target order, matching the native adjacency order.
  ids <- as.vector(t(ranks[, seq_len(s), drop = FALSE]))
  which_target <- rep(seq_len(targets), each = s)
  which_rank <- rep(seq_len(s), targets)
  groups <- split(seq_along(ids), ids)
  affected <- vector("list", n_ref)
  affected[as.integer(names(groups))] <- unname(groups)
  theta0 <- numeric(q)
  loo <- matrix(0, n_ref, q)
  for (component in seq_len(q)) {
    state <- .compensated_state_R(baseline[, component])
    theta0[component] <- (state[1L] + state[2L]) / targets
    for (i in seq_len(n_ref)) {
      at <- affected[[i]]
      terms <- as.numeric(state)
      if (length(at)) {
        pairs <- rbind(-baseline[which_target[at], component],
          deletions[cbind(which_target[at], which_rank[at], component)])
        terms <- c(terms, as.vector(pairs))
      }
      loo[i, component] <- .compensated_sum_R(terms) / targets
    }
  }
  list(theta0 = theta0, theta_loo = loo,
    info = list(s = s, T = targets, q = q, normalize = normalize))
}
