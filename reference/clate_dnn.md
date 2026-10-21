# Conditional LATE from externally supplied nuisance predictions

An augmented binary-IV estimator with rich adjustment and
reduced-dimensional reporting. Nuisance functions are supplied as
aligned predictions; the estimator fits no learners.

## Usage

``` r
clate_dnn(
  data = NULL,
  x_cols = NULL,
  y_col = NULL,
  d_col = NULL,
  a_col = NULL,
  e_col = NULL,
  m_y0_col = NULL,
  m_y1_col = NULL,
  m_d0_col = NULL,
  m_d1_col = NULL,
  X = NULL,
  X0 = NULL,
  x0 = NULL,
  Y = NULL,
  D = NULL,
  A = NULL,
  e = NULL,
  m_y0 = NULL,
  m_y1 = NULL,
  m_d0 = NULL,
  m_d1 = NULL,
  clip = 1e-08,
  one_sided = FALSE,
  distance = "l2",
  method = c("exact", "asymp"),
  k = 1L,
  s = NULL,
  tuning = NULL,
  tuning_mode = c("default", "grid"),
  normalize = TRUE,
  se = TRUE,
  jackknife = c("fast", "loo", "none"),
  denominator_tol = sqrt(.Machine$double.eps),
  row_ids = NULL,
  training = "supplied",
  ctx = NULL,
  return = c("W")
)
```

## Arguments

- data:

  Optional data frame or tibble. Do not combine with X or direct sample
  vectors.

- x_cols:

  Reporting-covariate columns in data, in target-coordinate order.

- y_col, d_col, a_col, e_col:

  Data-frame columns for outcome, treatment, instrument and instrument
  propensity.

- m_y0_col, m_y1_col:

  Data-frame outcome predictions by instrument arm.

- m_d0_col, m_d1_col:

  Data-frame treatment predictions by instrument arm. The zero-arm
  prediction may be omitted with one_sided = TRUE.

- X:

  Numeric reporting matrix, observations in rows.

- X0:

  Numeric target matrix, targets in rows.

- x0:

  Numeric vector for a single target. Supply exactly one of X0,x0.

- Y:

  Finite numeric outcome vector.

- D:

  Numeric binary treatment vector (0/1).

- A:

  Numeric binary instrument vector (0/1).

- e:

  Instrument propensity predictions in \\\[0,1\]\\, before clipping.

- m_y0, m_y1:

  Outcome predictions conditional on instrument arm 0 and 1.

- m_d0, m_d1:

  Treatment predictions conditional on instrument arm 0 and 1, in
  \\\[0,1\]\\. m_d0 may be omitted only with one_sided = TRUE.

- clip:

  Propensity clipping threshold in (0,0.5), with `1 - clip < 1` in
  double precision.

- one_sided:

  Explicit structural restriction that treatment is impossible under
  instrument 0. Requires D \<= A and zero m_d0 if supplied.

- distance:

  Existing distance specification, registered name or function. Geometry
  is frozen across deletions; input position resolves distance ties.

- method:

  One of "exact" or "asymp", as in
  [`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md)
  and
  [`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md).

- k:

  Number of neighbors averaged within each size-s subsample, 1 \<= k \<=
  s.

- s:

  Deterministic integer order. Covariance requires s \< n.

- tuning:

  Optional existing tuning specification, mutually exclusive with s. If
  neither is provided, the package paper rule resolves the order.

- tuning_mode:

  "default" retains a single fit; "grid" also evaluates feasible tuning
  candidates using the same frozen signals and ranks.

- normalize:

  Logical normalization control passed to the asymptotic weight rule, as
  in
  [`weights_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/weights_dnn.md).
  Exact weights are normalized by construction.

- se:

  Whether to compute covariance.

- jackknife:

  "fast" uses the existing prefix/suffix engine; "loo" uses explicit
  frozen-rank deletions; "none" disables covariance.

- denominator_tol:

  Nonnegative numerical cancellation tolerance. Changes a diagnostic
  flag only, never the ratio or covariance.

- row_ids:

  Optional aligned unique character/numeric identities.

- training:

  Descriptive declaration: supplied, oracle, full_sample, oob or
  cross_fit. It is not verified from predictions.

- ctx:

  Optional context passed to the existing distance machinery.

- return:

  Optional extras: W, ranks, ord, dx, ess, signals, deletions,
  pseudo_values, neighbors or tuning. Unknown values warn and are
  ignored.

## Value

A clate_result inheriting affine_ratio_result. Component roles are
reduced_form and first_stage. Optional nuisance output is the complete
iv_signals construction record, including raw/used propensity.

## Details

The documented sampling framework uses independent observations; cluster
and survey-design inference are not supplied. Input row order breaks
distance ties, so permuting tied observations can change estimates.

Nuisance predictions use all pre-instrument adjustment variables, while
weights use the reporting coordinates X. The estimate is the ratio of
the adjusted local reduced form and first stage, using common weights.

The accompanying
[paper](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z)
gives the IV identifying assumptions and conditions for pointwise
inference. A training label or a computational choice of k or weight
method does not establish those conditions. Inspect the local first
stage alongside the ratio.

Covariance keeps paired variation of both frozen signals. The
residualized IV-GMM workflow uses a different moment and can have a
different target. Reporting uses the common generics; overlap panels
describe the instrument and influence panels retain both augmented
signals.

## Examples

``` r
d <- data.frame(
  x = 1:6, Y = c(1, 2, 0, 3, 2, 1),
  D = c(0, 1, 0, 1, 0, 1), A = c(0, 1, 1, 1, 0, 1),
  e = rep(.5, 6), y0 = rep(1, 6), y1 = rep(2, 6), d1 = rep(.6, 6)
)
fit <- clate_dnn(
  data = d, x_cols = "x", y_col = "Y", d_col = "D", a_col = "A",
  e_col = "e", m_y0_col = "y0", m_y1_col = "y1", m_d1_col = "d1",
  one_sided = TRUE, x0 = 3, s = 2
)
fit$components
#>        numerator denominator
#> [1,] -0.06666667   0.5733333
fit$theta
#> [1] -0.1162791
clate_dnn(
  X = matrix(d$x, ncol = 1), Y = d$Y, D = d$D, A = d$A,
  e = d$e, m_y0 = d$y0, m_y1 = d$y1, m_d1 = d$d1,
  one_sided = TRUE, x0 = 3, s = 2
)$theta
#> [1] -0.1162791
```
