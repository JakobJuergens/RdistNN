#' RdistNN: DNN Estimation for CATE and Nonparametric Regression
#'
#' The initial public release focuses on conditional average treatment effects
#' with [cate_dnn()] and nonparametric regression with [np_dnn()], using
#' distributional nearest-neighbor (DNN) weights and delete-1 jackknife standard
#' errors. Nuisance predictions for CATE are supplied externally.
#' The underlying toolkit combines:
#' - distance specifications (`dnn_distance`) for ranking observations around a
#'   target;
#' - rank-weight specifications (`rank_weight`) for converting ranks to weights;
#' - score objects (`ustat_score`) for flexible kernel-based moment evaluation.
#'
#' @section Main User Workflows:
#' - Compute DNN weights with [weights_dnn()], [weights_dnn_exact()], and
#'   [weights_dnn_asymp()].
#' - Resolve paper-default kernel orders with [tuning_dnn_paper()] and
#'   [tuning_eval()].
#' - Estimate local CATE targets with [cate_dnn()].
#' - Summarize/predict CATE fits with [summary.cate_result()],
#'   [confint.cate_result()], and [predict.cate_result()].
#' - Estimate local nonparametric means with [np_dnn()].
#' - Use two-scale signed-weight regression with [np_tdnn()]; see its help for
#'   order constraints and additional smoothness requirements.
#' - Inspect fitted local estimators with [diagnostics()].
#' - Plot CATE and regression estimates with [plot_estimates()] and inspect
#'   local weights with [plot_weights()] (requires ggplot2).
#'
#' @section Experimental and Advanced Workflows:
#' All PLQR functionality is **experimental**, including [plqr_dnn()], its
#' scores, demo helpers, methods, and plots. The current wrapper evaluates
#' moments at a supplied coefficient; it does not solve for local coefficients.
#' Local IV-GMM functionality is also experimental. These workflows are outside
#' the CATE and regression focus of the initial public release.
#' - Estimate conditional LATE from supplied nuisance predictions with
#'   [clate_dnn()], or localize two supplied signals with [affine_ratio_dnn()].
#'   [iv_augmented_signals()] constructs the augmented binary-IV signals.
#'   Reporting methods and [diagnostics()] retain components and statuses;
#'   [ratio_confidence_set()] provides separate set-valued moment inversion.
#'   [clate_demo_data()] supplies a bounded synthetic benchmark with oracle
#'   predictions. All conditional-LATE functionality is experimental.
#' - Build custom generalized U-statistics with [new_score()],
#'   [score_eval_batch()], and [score_diff_path()].
#' - Estimate uncertainty with [jackknife_fast_rank()] and
#'   [jackknife_slow_rank()].
#'
#' @section Documentation Map:
#' - Start with `vignette("getting-started", package = "RdistNN")`.
#' - For step-by-step extension patterns for custom distances, rank weights,
#'   and scores, see `vignette("extending", package = "RdistNN")`.
#' - For modular kernel-order tuning rules and registry usage, see
#'   `vignette("kernel-order-tuning", package = "RdistNN")`.
#' - For fitted-object diagnostics and reporting, see
#'   `vignette("diagnostics-reporting", package = "RdistNN")`.
#' - For distance/rank weighting details see
#'   `vignette("dnn-weights", package = "RdistNN")`.
#' - For estimator workflows see
#'   `vignette("conditional-z-estimation", package = "RdistNN")`, including
#'   experimental PLQR details.
#' - For supplied-prediction conditional LATE and paired inference, see
#'   `vignette("conditional-late", package = "RdistNN")` and
#'   `vignette("conditional-late-inference", package = "RdistNN")`.
#' - For jackknife details see
#'   `vignette("jackknife-gen-ustats", package = "RdistNN")`.
#'
#' @docType package
#' @name RdistNN-package
#' @aliases RdistNN
#' @seealso [weights_dnn()], [cate_dnn()], [np_dnn()], [plqr_dnn()],
#'   [diagnostics()], [summary.cate_result()], [predict.cate_result()],
#'   [new_score()], [distance_l2()]
#' @useDynLib RdistNN, .registration = TRUE
#' @importFrom stats pbeta plogis rbinom rnorm runif setNames
"_PACKAGE"
