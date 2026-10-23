# RdistNN: DNN Estimation for Conditional Effects and Regression

Supported workflows include conditional average treatment effects with
[`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md),
nonparametric regression with
[`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md),
conditional LATE with
[`clate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/clate_dnn.md)
and ordinary supplied-signal ratios with
[`affine_ratio_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/affine_ratio_dnn.md),
using distributional nearest-neighbor (DNN) weights and delete-1
jackknife standard errors. Nuisance predictions for causal estimators
are supplied externally; the documented inference conditions remain
estimator-specific. The underlying toolkit combines:

- distance specifications (`dnn_distance`) for ranking observations
  around a target;

- rank-weight specifications (`rank_weight`) for converting ranks to
  weights;

- score objects (`ustat_score`) for flexible kernel-based moment
  evaluation.

## Main User Workflows

- Compute DNN weights with
  [`weights_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/weights_dnn.md),
  [`weights_dnn_exact()`](https://JakobJuergens.github.io/RdistNN/reference/weights_dnn_exact.md),
  and
  [`weights_dnn_asymp()`](https://JakobJuergens.github.io/RdistNN/reference/weights_dnn_asymp.md).

- Resolve paper-default kernel orders with
  [`tuning_dnn_paper()`](https://JakobJuergens.github.io/RdistNN/reference/tuning_dnn_paper.md)
  and
  [`tuning_eval()`](https://JakobJuergens.github.io/RdistNN/reference/tuning_eval.md).

- Estimate local CATE targets with
  [`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md).

- Report and reuse CATE fits with
  [`summary.cate_result()`](https://JakobJuergens.github.io/RdistNN/reference/summary.cate_result.md),
  [`confint.cate_result()`](https://JakobJuergens.github.io/RdistNN/reference/confint.cate_result.md),
  [`vcov.cate_result()`](https://JakobJuergens.github.io/RdistNN/reference/vcov.cate_result.md)
  and
  [`predict.cate_result()`](https://JakobJuergens.github.io/RdistNN/reference/predict.cate_result.md).

- Estimate local nonparametric means with
  [`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md).

- Use two-scale signed-weight regression with
  [`np_tdnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_tdnn.md);
  see its help for order constraints and additional smoothness
  requirements. Both regression classes share reporting, joint
  covariance and prediction methods documented in
  [np_result](https://JakobJuergens.github.io/RdistNN/reference/np_result.md).

- Estimate conditional LATE from supplied nuisance predictions with
  [`clate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/clate_dnn.md),
  or localize two supplied signals with
  [`affine_ratio_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/affine_ratio_dnn.md).
  [`iv_augmented_signals()`](https://JakobJuergens.github.io/RdistNN/reference/iv_augmented_signals.md)
  constructs the augmented binary-IV signals. Reporting methods and
  [`diagnostics()`](https://JakobJuergens.github.io/RdistNN/reference/diagnostics.md)
  retain components and statuses.
  [`clate_demo_data()`](https://JakobJuergens.github.io/RdistNN/reference/clate_demo_data.md)
  supplies a bounded synthetic benchmark with oracle predictions.

- Inspect fitted local estimators with
  [`diagnostics()`](https://JakobJuergens.github.io/RdistNN/reference/diagnostics.md).

- Plot CATE and regression estimates with
  [`plot_estimates()`](https://JakobJuergens.github.io/RdistNN/reference/plot_estimates.md)
  and inspect local weights with
  [`plot_weights()`](https://JakobJuergens.github.io/RdistNN/reference/plot_weights.md)
  (requires ggplot2).

## Exact distance ties

Distance-based estimators average rank weights within exactly equal
computed distances, separately for every target, including each
deletion's surviving members. Block mass is preserved for signed and
unnormalized rules. Equality uses represented values without a
tolerance. Bare ranks/orders remain strict; adding `dx` metadata enables
averaging while preserving supplied permutations. Calls with ties warn
once: the accompanying theory assumes a continuous distance distribution
and does not establish inference under persistent ties. Earlier retained
fit schemas require refitting before reporting or reuse.

## Experimental and Advanced Workflows

All PLQR functionality is **experimental**, including
[`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md),
its scores, demo helpers, methods, and plots. The current wrapper
evaluates moments at a supplied coefficient; it does not solve for local
coefficients. Local IV-GMM functionality is also experimental. These
workflows are outside the supported workflow focus. The optional
[`ratio_confidence_set()`](https://JakobJuergens.github.io/RdistNN/reference/ratio_confidence_set.md)
API remains experimental and provides separate set-valued moment
inversion.

- Build custom generalized U-statistics with
  [`new_score()`](https://JakobJuergens.github.io/RdistNN/reference/new_score.md),
  [`score_eval_batch()`](https://JakobJuergens.github.io/RdistNN/reference/score_eval_batch.md),
  and
  [`score_diff_path()`](https://JakobJuergens.github.io/RdistNN/reference/score_diff_path.md).

- Estimate uncertainty with
  [`jackknife_fast_rank()`](https://JakobJuergens.github.io/RdistNN/reference/jackknife_fast_rank.md)
  and
  [`jackknife_slow_rank()`](https://JakobJuergens.github.io/RdistNN/reference/jackknife_slow_rank.md).

## Random numbers

Explicit demo and simulation seeds restore `.Random.seed` and
[`RNGkind()`](https://rdrr.io/r/base/Random.html); `NULL` uses the
current stream. R's saved state does not include the hidden Box-Muller
normal cache; use the default normal generator for saved-state
reproducibility.

## Native threads

Compiled routines use at most two OpenMP threads by default. Set
`options(RdistNN.num_threads = 1L)` for one native thread, or supply a
larger positive integer to opt into more threads. Lower OpenMP runtime
limits are respected; builds without OpenMP run serially. This option
does not change other packages' thread settings.

## Documentation Map

- Start with
  [`vignette("getting-started", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/getting-started.md).

- For step-by-step extension patterns for custom distances, rank
  weights, and scores, see
  [`vignette("extending", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/extending.md).

- For modular kernel-order tuning rules and registry usage, see
  [`vignette("kernel-order-tuning", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/kernel-order-tuning.md).

- For fitted-object diagnostics and reporting, see
  [`vignette("diagnostics-reporting", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/diagnostics-reporting.md).

- For distance/rank weighting details see
  [`vignette("dnn-weights", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/dnn-weights.md).

- For experimental PLQR moments and coefficient uncertainty, see
  [`vignette("conditional-z-estimation", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/conditional-z-estimation.md).

- For supplied-prediction conditional LATE and paired inference, see
  [`vignette("conditional-late", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/conditional-late.md)
  and
  [`vignette("conditional-late-inference", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/conditional-late-inference.md).

- For jackknife details see
  [`vignette("jackknife-gen-ustats", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/jackknife-gen-ustats.md).

- For statistical background, assumptions and theoretical results, see
  the accompanying
  [paper](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z).

## See also

[`weights_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/weights_dnn.md),
[`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md),
[`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md),
[`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md),
[`diagnostics()`](https://JakobJuergens.github.io/RdistNN/reference/diagnostics.md),
[`summary.cate_result()`](https://JakobJuergens.github.io/RdistNN/reference/summary.cate_result.md),
[`predict.cate_result()`](https://JakobJuergens.github.io/RdistNN/reference/predict.cate_result.md),
[`new_score()`](https://JakobJuergens.github.io/RdistNN/reference/new_score.md),
[`distance_l2()`](https://JakobJuergens.github.io/RdistNN/reference/distance_l2.md)

## Author

**Maintainer**: Jakob R. Juergens <info@jakobjuergens.com>

Authors:

- Jakob R. Juergens <info@jakobjuergens.com>

Other contributors:

- Zakareya Chehime (Quality assurance, testing and validation)
  \[contributor\]
