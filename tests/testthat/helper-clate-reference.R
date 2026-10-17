.clate_test_predictions <- function() {
  list(
    Y = c(-1, 2, 1, 0, 3, 2, 4, -2),
    D = c(0, 0, 1, 1, 0, 1, 1, 0),
    A = c(0, 1, 0, 1, 1, 0, 1, 0),
    e = c(.25, .5, .75, .4, .6, .3, .7, .5),
    m_y0 = seq(-.4, .3, length.out = 8L),
    m_y1 = seq(.6, 1.3, length.out = 8L),
    m_d0 = rep(.2, 8L), m_d1 = rep(.7, 8L)
  )
}

# Independent finite-subset oracle: no package weights or deletion helpers.
.clate_subset_average <- function(X, x0, signal, s) {
  subsets <- utils::combn(seq_len(nrow(X)), s)
  distance <- rowSums(sweep(X, 2L, x0, "-")^2)
  winners <- apply(subsets, 2L, function(idx) {
    idx[order(distance[idx], idx)[1L]]
  })
  colMeans(signal[winners, , drop = FALSE])
}

.clate_subset_reference <- function(X, x0, signal, s) {
  n <- nrow(X)
  full <- .clate_subset_average(X, x0, signal, s)
  loo <- t(vapply(seq_len(n), function(j) {
    .clate_subset_average(X[-j, , drop = FALSE], x0,
                          signal[-j, , drop = FALSE], s)
  }, numeric(2L)))
  pv <- sweep(-(n - 1) * loo, 2L, n * full, "+")
  centered <- sweep(pv, 2L, colMeans(pv), "-")
  list(components = full, loo = loo, pv = pv,
       vcov = crossprod(centered) / (n * (n - 1)))
}

.clate_fit_arguments <- function() {
  c(.clate_test_predictions(), list(
    X = cbind(x1 = seq_len(8L) / 8, x2 = sin(seq_len(8L))),
    X0 = rbind(c(.2, 0), c(.8, .2)), s = 3L, return = character()
  ))
}
