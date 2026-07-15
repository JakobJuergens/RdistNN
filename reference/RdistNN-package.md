# RdistNN: Distance-Based Building Blocks for Localized U-Statistics

`RdistNN` provides a modular toolkit for distance-localized estimation
based on generalized U-statistics. It combines:

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

- Summarize/predict CATE fits with
  [`summary.cate_result()`](https://JakobJuergens.github.io/RdistNN/reference/summary.cate_result.md),
  [`confint.cate_result()`](https://JakobJuergens.github.io/RdistNN/reference/confint.cate_result.md),
  and
  [`predict.cate_result()`](https://JakobJuergens.github.io/RdistNN/reference/predict.cate_result.md).

- Estimate local nonparametric means with
  [`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md).

- Estimate local PLQR coefficient surfaces and orthogonal moments with
  [`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md).

- Summarize/predict PLQR fits with
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html), and
  [`predict()`](https://rdrr.io/r/stats/predict.html) methods.

- Inspect fitted local estimators with
  [`diagnostics()`](https://JakobJuergens.github.io/RdistNN/reference/diagnostics.md).

- Build custom generalized U-statistics with
  [`new_score()`](https://JakobJuergens.github.io/RdistNN/reference/new_score.md),
  [`score_eval_batch()`](https://JakobJuergens.github.io/RdistNN/reference/score_eval_batch.md),
  and
  [`score_diff_path()`](https://JakobJuergens.github.io/RdistNN/reference/score_diff_path.md).

- Estimate uncertainty with
  [`jackknife_fast_rank()`](https://JakobJuergens.github.io/RdistNN/reference/jackknife_fast_rank.md)
  and
  [`jackknife_slow_rank()`](https://JakobJuergens.github.io/RdistNN/reference/jackknife_slow_rank.md).

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

- For estimator workflows see
  [`vignette("conditional-z-estimation", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/conditional-z-estimation.md).

- For jackknife details see
  [`vignette("jackknife-gen-ustats", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/jackknife-gen-ustats.md).

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

**Maintainer**: Jakob Juergens <info@jakobjuergens.com>

Authors:

- Jakob Juergens <info@jakobjuergens.com>
