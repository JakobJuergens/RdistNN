# Independent numeric sorting reference; never consult fitted signals or the
# production ordering helper when constructing a test expectation.
.tie_reference_order <- function(dx, X, row_ids = seq_len(nrow(X))) {
  do.call(order, c(list(dx), unname(as.list(as.data.frame(X))), list(row_ids)))
}

.tie_reference_weights <- function(X, X0, k, s) {
  sets <- utils::combn(nrow(X), s, simplify = FALSE)
  vapply(seq_len(nrow(X0)), function(j) {
    dx <- sqrt(rowSums(sweep(X, 2L, X0[j, ], "-")^2))
    w <- numeric(nrow(X))
    for (ids in sets) {
      ord <- .tie_reference_order(dx[ids], X[ids, , drop = FALSE], ids)
      take <- ids[ord][seq_len(k)]
      w[take] <- w[take] + 1 / (length(sets) * k)
    }
    w
  }, numeric(nrow(X)))
}

.tie_fixture <- function() {
  list(X = rbind(c(1, 0), c(-1, 0), c(0, 1), c(0, -1), c(2, 0), c(0, 2)),
    X0 = rbind(c(0, 0), c(.25, 0), c(0, 0)),
    y = c(7, -3, 9, 1, 4, -2), perm = c(6L, 1L, 4L, 3L, 2L, 5L))
}

.tie_fit <- function(kind, d, at = seq_len(nrow(d$X)), jackknife = "fast", ...) {
  common <- list(X = d$X[at, , drop = FALSE], X0 = d$X0,
    return = c("W", "ranks", "ord", "dx", "ess"), jackknife = jackknife)
  signal <- d$y[at]
  if (kind == "np") {
    do.call(np_dnn, c(common, list(Y = signal), list(...)))
  } else if (kind == "tdnn") {
    do.call(np_tdnn, c(common, list(Y = signal, s1 = 2L, s2 = 4L, dim = 2L), list(...)))
  } else if (kind == "cate") {
    do.call(cate_dnn, c(common, list(pseudo_outcome = signal), list(...)))
  } else if (kind == "clate") {
    A <- rep(0:1, length.out = nrow(d$X))[at]
    do.call(clate_dnn, c(common, list(Y = A * signal, D = A, A = A,
      e = rep(.5, length(at)), m_y0 = rep(0, length(at)), m_y1 = signal,
      m_d0 = rep(0, length(at)), m_d1 = rep(1, length(at))), list(...)))
  } else {
    stop("Unknown tie fixture estimator")
  }
}
