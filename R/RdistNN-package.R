#' RdistNN: DNN Estimation for Conditional Effects and Regression
#'
#' Supported workflows include conditional average treatment effects with
#' [cate_dnn()], nonparametric regression with [np_dnn()], conditional LATE with
#' [clate_dnn()] and ordinary supplied-signal ratios with [affine_ratio_dnn()],
#' using distributional nearest-neighbor (DNN) weights and delete-1 jackknife
#' standard errors. Nuisance predictions for causal estimators are supplied
#' externally; the documented inference conditions remain estimator-specific.
#' The underlying toolkit combines:
#' - distance specifications (`rdistnn_distance`) for ranking observations around a
#'   target;
#' - rank-weight specifications (`rdistnn_rank_weight`) for converting ranks to weights;
#' - score objects (`rdistnn_ustat_score`) for flexible kernel-based moment evaluation.
#'
#' @section Main User Workflows:
#' - Compute DNN weights with [weights_dnn()], [weights_dnn_exact()], and
#'   [weights_dnn_asymp()].
#' - Resolve paper-default kernel orders with [tuning_dnn_paper()] and
#'   [tuning_eval()].
#' - Estimate local CATE targets with [cate_dnn()].
#' - Report and reuse CATE fits with [summary.rdistnn_cate_result()],
#'   [confint.rdistnn_cate_result()], [vcov.rdistnn_cate_result()] and [predict.rdistnn_cate_result()].
#' - Estimate local nonparametric means with [np_dnn()].
#' - Use two-scale signed-weight regression with [np_tdnn()]; see its help for
#'   order constraints and additional smoothness requirements.
#'   Both regression classes share reporting, joint covariance and prediction
#'   methods documented in [rdistnn_np_result].
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
#' @section Exact distance ties:
#' Distance-based estimators average rank weights within exactly equal computed
#' distances, separately for every target, including each deletion's surviving
#' members. Block mass is preserved for signed and unnormalized rules. Equality
#' uses represented values without a tolerance. Bare ranks/orders remain strict;
#' adding `dx` metadata enables averaging while preserving supplied permutations.
#' Calls with ties warn once: the accompanying theory assumes a continuous
#' distance distribution and does not establish inference under persistent ties.
#' During experimental development, refit saved models after package updates.
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
#' Built-in distance, rank-weight and accumulation routines preserve the caller's
#' RNG state in both R and compiled execution, including leaving `.Random.seed`
#' absent when it was absent. These deterministic calculations need no seed.
#' Explicit demo and simulation seeds restore `.Random.seed` and [RNGkind()];
#' `NULL` uses the current stream. R's saved state does not include the hidden
#' Box-Muller normal cache; use the default normal generator for saved-state
#' reproducibility.
#'
#' @section Native threads:
#' Eligible built-in calculations use C++ automatically when the corresponding
#' native routine is available. For engine comparisons or troubleshooting, set
#' `options(RdistNN.use_compiled = FALSE)` to select the package's R paths;
#' `TRUE` restores automatic native dispatch. Source installation still requires
#' a C++17 compiler. Native execution errors propagate rather than triggering an
#' automatic retry in R. During long native calculations, interruption can be
#' delayed until control returns to R.
#'
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
#'   [diagnostics()], [summary.rdistnn_cate_result()], [predict.rdistnn_cate_result()],
#'   [new_score()], [distance_l2()]
#' @useDynLib RdistNN, .registration = TRUE
#' @importFrom stats pbeta plogis rbinom rnorm runif setNames
"_PACKAGE"
