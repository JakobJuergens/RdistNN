# Internal bounded-model reference learner, not a general automatic fit wrapper.
# The dictionary and ambient loss normalizer stay fixed across at most two
# deletions. Three small Gram solves per occupied cell supply five predictions.
.clate_ridge_predictions <- function(
  data, M, lambda, normalizer_n = nrow(data), delete = integer()
) {
  if (!is.data.frame(data)) .abort_bad_arg("data", "must be a data.frame")
  M <- .affine_integer(M, "M")
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda <= 0) {
    .abort_bad_arg("lambda", "must be a finite strictly positive numeric scalar")
  }
  normalizer_n <- .affine_integer(normalizer_n, "normalizer_n")
  J <- as.double(M)^4
  if (!is.finite(J) || 3 * J > .Machine$integer.max) {
    .abort_bad_arg("M", "produces a dictionary outside the R matrix index range")
  }
  X <- .dnn_col_mat(c("X1", "X2"), data, "X columns", out_arg = "X")
  Z <- .dnn_col_mat(paste0("Z", 1:4), data, "Z columns", out_arg = "Z")
  n <- nrow(X)
  if (any(X < 0 | X > 1) || any(Z < 0 | Z > 1)) {
    .abort_bad_arg("data", "must have reporting and adjustment covariates in [0, 1]")
  }
  .abort_if_na_or_nonfinite(X, "X", quoted = TRUE)
  .abort_if_na_or_nonfinite(Z, "Z", quoted = TRUE)
  Y <- .affine_column("Y", data, n, "Y")
  A <- .validate_binary01(.affine_column("A", data, n, "A"), "A")
  D <- .validate_binary01(.affine_column("D", data, n, "D"), "D")
  if (any(abs(Y) > 1)) .abort_bad_arg("Y", "must be bounded by one in the reference model")
  if (!is.numeric(delete) || !is.null(dim(delete)) || length(delete) > 2L ||
      anyNA(delete) || any(!is.finite(delete)) || any(delete != floor(delete)) ||
      any(delete < 1 | delete > n) || anyDuplicated(delete)) {
    .abort_bad_arg("delete", "must contain at most two distinct retained-sample indices in 1..n")
  }
  delete <- sort(as.integer(delete))
  retained <- !(seq_len(n) %in% delete)
  ids <- .affine_row_ids(data$row_id, n)
  cells <- as.integer(1 + pmin(floor(Z * M), M - 1) %*% (M^(0:3)))
  q <- cbind(1, sqrt(12) * (X[, 1L] - .5), sqrt(12) * (X[, 2L] - .5))
  labels <- c("e", "m_y0", "m_y1", "m_d0", "m_d1")
  beta <- matrix(0, 3 * J, 5L, dimnames = list(NULL, labels))
  raw <- matrix(0, n, 5L, dimnames = list(NULL, labels))
  for (cell in unique(cells)) {
    evaluate <- which(cells == cell)
    train <- which(cells == cell & retained)
    block <- matrix(0, 3L, 5L)
    for (arm in c(-1L, 0L, 1L)) {
      idx <- if (arm == -1L) train else train[A[train] == arm]
      if (!length(idx)) next
      P <- sqrt(J) * q[idx, , drop = FALSE]
      columns <- if (arm == -1L) 1L else if (arm == 0L) c(2L, 4L) else c(3L, 5L)
      response <- if (arm == -1L) matrix(A[idx], ncol = 1L) else cbind(Y[idx], D[idx])
      gram <- crossprod(P) / normalizer_n + lambda * diag(3L)
      rhs <- crossprod(P, response) / normalizer_n
      solution <- tryCatch(solve(gram, rhs), error = function(e) {
        .abort_bad_arg("lambda", "gives a numerically singular ridge system; increase regularization")
      })
      block[, columns] <- solution
    }
    beta[3L * (cell - 1L) + 1:3, ] <- block
    raw[evaluate, ] <- sqrt(J) * q[evaluate, , drop = FALSE] %*% block
  }
  used <- raw
  used[, 1L] <- pmax(.25, pmin(.75, raw[, 1L]))
  used[, 2:3] <- pmax(-1, pmin(1, raw[, 2:3]))
  used[, 4:5] <- pmax(0, pmin(1, raw[, 4:5]))
  as_vectors <- function(x) stats::setNames(lapply(seq_len(ncol(x)), function(j) x[, j]), labels)
  list(
    predictions = as_vectors(used), raw_predictions = as_vectors(raw),
    coefficients = beta, row_ids = ids,
    training = list(
      mechanism = "retained_ridge_series", M = M, dictionary_dim = as.integer(3 * J),
      lambda = as.numeric(lambda), normalizer_n = normalizer_n, delete = delete,
      retained_n = sum(retained), loss = "ambient_normalized_ridge"
    )
  )
}
