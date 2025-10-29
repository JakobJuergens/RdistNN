# R/wdist_gaussian_l2.R
# -----------------------------------------------------------------------------
# Preimplemented distance-based weight: Gaussian kernel over L2 distances
# -----------------------------------------------------------------------------

#' Gaussian kernel distance-weight for L2 distances
#'
#' Constructs a \code{dist_weight} that maps L2 distances to weights via
#' \deqn{w_i \propto \exp\{ - \|x_i - x_0\|^2 / (2 h^2) \}.}
#' By default, this expects **squared** L2 distances (pairs nicely with
#' \code{distance_l2(squared = TRUE)} to avoid an extra \code{sqrt}).
#'
#' If \code{normalize = TRUE} (default), columns are normalized to sum to 1.
#' In extreme underflow cases (all weights numerically zero), it falls back to
#' a 1-NN assignment (puts all mass on the nearest point per target column).
#'
#' @param h Positive bandwidth.
#' @param squared Logical; if \code{TRUE}, inputs \code{dx} are interpreted as
#'   squared L2 distances. If \code{FALSE}, inputs are plain L2 distances and
#'   the kernel uses \eqn{\exp\{-\|x_i-x_0\|^2/(2h^2)\}} internally.
#' @param normalize Logical; if \code{TRUE} (default) normalize weights to sum 1.
#' @return A \code{dist_weight} object usable with \code{dist_weight_eval()} or
#'   \code{dist_weight_eval_batch()}.
#' @examples
#' # Pair with distance_l2(squared=TRUE)
#' # X  <- matrix(rnorm(100), 50, 2)
#' # X0 <- matrix(rnorm(4),   2,  2)
#' # DX <- distance_eval_batch(distance_l2(squared=TRUE), X, X0)
#' # dw <- wdist_gaussian_l2(h = 0.5, squared = TRUE)
#' # W  <- dist_weight_eval_batch(dw, DX)  # 50 x 2, columns sum to 1
#' @export
wdist_gaussian_l2 <- function(h, squared = TRUE, normalize = TRUE) {
  if (!is.numeric(h) || length(h) != 1L || !is.finite(h) || h <= 0)
    stop("`h` must be a positive finite scalar.")

  name <- sprintf("gaussian_l2(h=%.6g,%s)", h, if (squared) "sq" else "lin")

  new_dist_weight(
    name = name,
    # single-target: dx is numeric(n)
    eval_dx = function(dx, ctx = NULL) {
      x <- as.numeric(dx)
      if (squared) {
        z <- - x / (2 * h * h)           # x = ||·||^2
      } else {
        z <- - (x * x) / (2 * h * h)     # x = ||·||
      }
      w <- exp(z)

      if (isTRUE(normalize)) {
        s <- sum(w)
        if (!is.finite(s) || s <= .Machine$double.xmin) {
          # underflow fallback: 1-NN assignment
          w[] <- 0
          w[ which.min(x) ] <- 1
        } else {
          w <- w / s
        }
      }
      w
    },
    # batch: DX is numeric(n x m)
    eval_dx_batch = function(DX, ctx = NULL) {
      M <- as.matrix(DX)
      if (squared) {
        Z <- - M / (2 * h * h)
      } else {
        Z <- - (M * M) / (2 * h * h)
      }
      W <- exp(Z)

      if (isTRUE(normalize)) {
        cs <- colSums(W)
        # handle extreme underflow columns by 1-NN assignment
        bad <- !is.finite(cs) | cs <= .Machine$double.xmin
        if (any(bad)) {
          for (j in which(bad)) {
            wj <- numeric(nrow(W))
            # use the raw distance matrix M for nearest neighbor
            wj[ which.min(M[, j]) ] <- 1
            W[, j] <- wj
          }
          cs <- colSums(W)
        }
        W <- sweep(W, 2L, cs, "/")
      }
      W
    }
  )
}
