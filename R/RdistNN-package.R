#' RdistNN: Distance-Based Building Blocks for Localized U-Statistics
#'
#' `RdistNN` provides a modular toolkit for distance-localized estimation based
#' on generalized U-statistics. It combines:
#' - distance specifications (`dnn_distance`) for ranking observations around a
#'   target;
#' - rank-weight specifications (`rank_weight`) for converting ranks to weights;
#' - score objects (`ustat_score`) for flexible kernel-based moment evaluation.
#'
#' @section Main User Workflows:
#' - Compute DNN weights with [weights_dnn()], [weights_dnn_exact()], and
#'   [weights_dnn_asymp()].
#' - Estimate local CATE targets with [cate_dnn()].
#' - Summarize/predict CATE fits with [summary.cate_result()],
#'   [confint.cate_result()], and [predict.cate_result()].
#' - Estimate local nonparametric means with [np_dnn()].
#' - Evaluate local PLQR orthogonal moments with [plqr_dnn()].
#' - Build custom generalized U-statistics with [new_score()],
#'   [score_eval_batch()], and [score_diff_path()].
#' - Estimate uncertainty with [jackknife_fast_rank()] and
#'   [jackknife_slow_rank()].
#'
#' @section Documentation Map:
#' - Start with `vignette("getting-started", package = "RdistNN")`.
#' - For step-by-step extension patterns for custom distances, rank weights,
#'   and scores, see `vignette("extending", package = "RdistNN")`.
#' - For distance/rank weighting details see
#'   `vignette("dnn-weights", package = "RdistNN")`.
#' - For estimator workflows see
#'   `vignette("conditional-z-estimation", package = "RdistNN")`.
#' - For jackknife details see
#'   `vignette("jackknife-gen-ustats", package = "RdistNN")`.
#'
#' @docType package
#' @name RdistNN-package
#' @aliases RdistNN
#' @seealso [weights_dnn()], [cate_dnn()], [np_dnn()], [plqr_dnn()],
#'   [summary.cate_result()], [predict.cate_result()],
#'   [new_score()], [distance_l2()]
#' @useDynLib RdistNN, .registration = TRUE
#' @importFrom stats pbeta plogis rbinom rnorm runif setNames
"_PACKAGE"
