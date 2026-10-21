# Paired frozen-score inference for conditional ratios

This guide explains the covariance and status fields returned by
[`affine_ratio_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/affine_ratio_dnn.md)
and
[`clate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/clate_dnn.md).
Start with the [CLATE
workflow](https://JakobJuergens.github.io/RdistNN/articles/conditional-late.md)
for inputs and a complete example. For identification, nuisance-learning
conditions and asymptotic inference, see the accompanying
[paper](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z).

## A small paired-signal example

``` r

fit <- affine_ratio_dnn(
  X = matrix(1:4, ncol = 1), x0 = 0,
  numerator = c(2, 1, 4, 3), denominator = c(1, -1, 2, 1), s = 2
)
fit$components
#>      numerator denominator
#> [1,]         2         0.5
summary(fit)$coefficients
#>          estimate std.error statistic   p.value  conf.low conf.high
#> target_1        4  5.830952 0.6859943 0.4927167 -7.428456  15.42846
```

The numerator and denominator use the same weights, and their paired
covariance enters ratio uncertainty. Dividing the numerator SE by the
first stage or discarding the component cross-covariance gives a
different calculation. Deletions hold supplied signals and geometry
fixed; nuisance models are not refitted.

## Keep the variance scales distinct

| Returned field | Meaning |
|----|----|
| `components[, "numerator"]` | Reduced form for a CLATE fit |
| `components[, "denominator"]` | First stage for a CLATE fit |
| `component_vcov[[j]]` | Paired component covariance on the estimator scale |
| `theta_vcov[[j]]` | Ratio variance on the estimator scale |
| `theta_omega[[j]]` | `neff[j] * theta_vcov[[j]]`, local-limit covariance estimate |
| `se[j]` | Square root of the ratio variance |
| `neff[j]` | Inverse sum of squared reporting weights |

`theta_vcov` already includes the local sample-size scaling. Do not
divide it by `neff` again. Covariance lists are marginal across
reporting targets; they do not supply covariance between different
targets or simultaneous confidence bands.

``` r

summary(fit)
#>   localization: s = 2, s/n = 0.5
#> Components: numerator = numerator; denominator = denominator
#> Ratio statuses: ok=1
#>          estimate std.error statistic p.value conf.low conf.high
#> target_1        4     5.831     0.686  0.4927   -7.428     15.43
confint(fit)
#>           conf.low conf.high
#> target_1 -7.428456  15.42846
diagnostics(fit)
#> <dnn_diagnostics>
#>   targets: 1
#>   n_eff  : min=2.5714, median=2.5714, max=2.5714
#>   weights: max(top5 share)=1, max(weight)=0.5
#>   dist   : available for 1/1 target(s)
#>   instrument 1 mass: unavailable
#>   instrument 0 mass: unavailable
#>   clipped propensity mass: unavailable
#>   first-stage statistic: 0.577 to 0.577
#>   top-five squared jackknife-contribution share: 1 to 1
```

Effective sample size describes rank-weight concentration. It is not a
count of compliers. Physical radius, instrument-arm support, score
dispersion and first-stage strength can vary even when `neff` is the
same. The largest centered pseudo-value-square shares attribute fitted
variance to deletion indices, including rank changes; they differ from
raw weighted score contributions.

## Signed signals and unavailable intervals

Augmented treatment signals can be negative even with one-sided
compliance; they are not individual compliance indicators. A finite
negative aggregate denominator gives a signed ratio, without replacement
or trimming.

Inspect these fields before reporting a CLATE result:

| Field | What to check |
|----|----|
| `ratio_status` | `ok`, zero denominator, component overflow or ratio overflow |
| `negative_denominator` | A finite negative first stage |
| `near_zero_denominator` | Numerical cancellation; not a calibrated weak-IV test |
| `covariance_status` | Available, missing, zero-score or numerically failed covariance |
| `wald_available` | Whether the reported SE supports the computational Wald output |

An exactly zero first stage gives an undefined ratio while retaining
component estimates and available covariance. Exact zero score variance
gives SE zero and unavailable Wald bounds. A finite ratio may also have
unavailable uncertainty. These cases stay visible in summaries and
plots; inspect the first stage and status rather than interpreting
missing bounds as precision.

CLATE generally differs from an average of stratum-specific LATEs or an
ordinary residualized-IV coefficient. Use the paper’s estimand
definition when choosing the appropriate workflow.

## Ratio confidence sets

The experimental
[`ratio_confidence_set()`](https://JakobJuergens.github.io/RdistNN/reference/ratio_confidence_set.md)
helper uses the paired estimator-scale covariance. It can return bounded
intervals, singletons, half-lines, two rays, all real numbers or empty
sets. It remains separate from
[`confint()`](https://rdrr.io/r/stats/confint.html), which reports
pointwise Wald intervals.

``` r

ratio_confidence_set(fit$components, fit$component_vcov)
#> Paired ratio confidence sets
#>   target 1: two_rays (-Inf, -0.60543] U [1.3053, Inf) (ok)
sets <- ratio_confidence_set(
  rbind(c(2, 1), c(3, 0), c(0, 0)),
  list(diag(c(.5, 0)), diag(2), matrix(0, 2, 2))
)
sets
#> Paired ratio confidence sets
#>   target 1: bounded [0.6141, 3.3859] (ok)
#>   target 2: two_rays (-Inf, -1.1588] U [1.1588, Inf) (ok)
#>   target 3: all_real (-Inf, Inf) (ok)
vapply(sets, function(x) x$type, character(1))
#> [1] "bounded"  "two_rays" "all_real"
```

Finite endpoints are included. An `unresolved` result reports a
numerical status without asserting a set geometry; it differs from a
mathematically empty set. Materially indefinite or unavailable
covariance is rejected. See
[`?ratio_confidence_set`](https://JakobJuergens.github.io/RdistNN/reference/ratio_confidence_set.md)
for input validation and returned fields.

Neither a finite Wald interval nor this experimental inversion
establishes uniform coverage under weak first stages. The paper
describes the scope of its inference results. Numerical safeguards and
successful array tests do not establish the statistical assumptions for
a particular application.
