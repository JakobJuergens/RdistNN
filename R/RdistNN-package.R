#' RdistNN: DNN Estimation for Conditional Effects and Regression
#'
#' Supported workflows include conditional average treatment effects with
#' [cate_dnn()], nonparametric regression with [np_dnn()], conditional LATE with
#' [clate_dnn()] and ordinary supplied-signal ratios with [affine_ratio_dnn()],
#' using distributional nearest-neighbor (DNN) weights and delete-1 jackknife
#' standard errors. Nuisance predictions for causal estimators are supplied
#' externally; the documented inference conditions remain estimator-specific.
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
#' - Report and reuse CATE fits with [summary.cate_result()],
#'   [confint.cate_result()], [vcov.cate_result()] and [predict.cate_result()].
#' - Estimate local nonparametric means with [np_dnn()].
#' - Use two-scale signed-weight regression with [np_tdnn()]; see its help for
#'   order constraints and additional smoothness requirements.
#'   Both regression classes share reporting, joint covariance and prediction
#'   methods documented in [np_result].
#' - Estimate conditional LATE from supplied nuisance predictions with
#'   [clate_dnn()], or localize two supplied signals with [affine_ratio_dnn()].
#'   [iv_augmented_signals()] constructs the augmented binary-IV signals.
#'   Reporting methods and [diagnostics()] retain components and statuses.
#'   [clate_demo_data()] supplies a bounded synthetic benchmark with oracle
#'   predictions.
#' - Inspect fitted local estimators with [diagnostics()].
#' - Plot CATE and regression estimates with [plot_estimates()] and inspect
#'   local weights with [plot_weights()] (requires ggplot2).
#'
#' @section Identical localization covariates:
#' Covariate-aware estimators average weights within exactly identical raw
#' covariate blocks, including every jackknife deletion's surviving members.
#' The package warns because this permutation-invariant modification is not
#' covered by the paper's theory. It can change observation weights, support
#' and effective sample size. See [weights_matrix()] for the rule; rank-only
#' utilities without covariates cannot infer duplicate groups. Refit tied saved
#' objects created before this policy before reporting or reusing them.
#'
#' @section Experimental and Advanced Workflows:
#' All PLQR functionality is **experimental**, including [plqr_dnn()], its
#' scores, demo helpers, methods, and plots. The current wrapper evaluates
#' moments at a supplied coefficient; it does not solve for local coefficients.
#' Local IV-GMM functionality is also experimental. These workflows are outside
#' the supported workflow focus. The optional [ratio_confidence_set()] API
#' remains experimental and provides separate set-valued moment inversion.
#' - Build custom generalized U-statistics with [new_score()],
#'   [score_eval_batch()], and [score_diff_path()].
#' - Estimate uncertainty with [jackknife_fast_rank()] and
#'   [jackknife_slow_rank()].
#'
#' @section Random numbers:
#' Explicit demo and simulation seeds restore `.Random.seed` and [RNGkind()];
#' `NULL` uses the current stream. R's saved state does not include the hidden
#' Box-Muller normal cache; use the default normal generator for saved-state
#' reproducibility.
#'
#' @section Native threads:
#' Compiled routines use at most two OpenMP threads by default. Set
#' `options(RdistNN.num_threads = 1L)` for one native thread, or supply a larger
#' positive integer to opt into more threads. Lower OpenMP runtime limits are
#' respected; builds without OpenMP run serially. This option does not change
#' other packages' thread settings.
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
#' - For experimental PLQR moments and coefficient uncertainty, see
#'   `vignette("conditional-z-estimation", package = "RdistNN")`.
#' - For supplied-prediction conditional LATE and paired inference, see
#'   `vignette("conditional-late", package = "RdistNN")` and
#'   `vignette("conditional-late-inference", package = "RdistNN")`.
#' - For jackknife details see
#'   `vignette("jackknife-gen-ustats", package = "RdistNN")`.
#' - For statistical background, assumptions and theoretical results, see
#'   the accompanying [paper](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z).
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
