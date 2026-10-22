#' Number of training observations
#'
#' @param object A fitted result from [np_dnn()], [np_tdnn()] or [cate_dnn()].
#' @param ... Not used.
#' @return The integer number of training observations. This is distinct from
#'   the number of reporting profiles and the local effective sample size.
#' @examples
#' fit <- np_dnn(X = matrix(1:12), Y = sin(1:12), x0 = 4, s = 3)
#' nobs(fit)
#' @importFrom stats nobs
#' @export
nobs.np_result <- function(object, ...) {
  .dnn_check_tie_model(object$.model)
  object$meta$n
}

#' @rdname nobs.np_result
#' @export
nobs.cate_result <- function(object, ...) {
  .dnn_check_tie_model(object$.model)
  object$meta$n
}
