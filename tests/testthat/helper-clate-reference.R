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
.clate_subset_average <- function(X, x0, signal, s, k = 1L) {
  crossprod(.tie_reference_weights(X, matrix(x0, nrow = 1L), k, s), signal)[1L, ]
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

.ratio_fixture <- function(...) {
  args <- list(X = matrix(1:4, ncol = 1L), x0 = 0,
               numerator = c(2, 1, 4, 3), denominator = c(1, -1, 2, 1),
               s = 2L, return = character())
  extra <- list(...)
  args[names(extra)] <- extra
  do.call(affine_ratio_dnn, args)
}

.clate_test_frame <- function(a) {
  data.frame(x1 = a$X[, 1L], x2 = a$X[, 2L], a[names(.clate_test_predictions())])
}

.clate_frame_arguments <- function(a, data) {
  list(data = data, x_cols = c("x1", "x2"), X0 = a$X0, s = a$s,
       y_col = "Y", d_col = "D", a_col = "A", e_col = "e",
       m_y0_col = "m_y0", m_y1_col = "m_y1", m_d0_col = "m_d0", m_d1_col = "m_d1")
}
.clate_alignment_tuning <- function() {
  new_tuning(name = "clate_alignment", eval = function(...) list(s = 3L),
    eval_grid = function(...) data.frame(lambda = c(.5, 1, 1.1, 2, 4),
                                        s = c(2L, 3L, 3L, 6L, 8L)))
}
