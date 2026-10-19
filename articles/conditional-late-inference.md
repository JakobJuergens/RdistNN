# Paired frozen-score inference for conditional ratios

The supported
[`affine_ratio_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/affine_ratio_dnn.md)
and
[`clate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/clate_dnn.md)
workflows average two signals using common DNN weights, then take their
ratio. Exact weights and one neighbor per subsample are the defaults;
asymptotic DNN and general `k` use the same paired calculation. The
worked example below is exact. This article explains the paired
covariance calculation and its numerical diagnostics. Start with
[Conditional LATE from supplied nuisance
predictions](https://JakobJuergens.github.io/RdistNN/articles/conditional-late.md)
for the economic target and input roles.

## A four-observation calculation

The following fixed arrays isolate the computation from nuisance
learning. The signed denominator signal is intentional.

``` r

fit <- affine_ratio_dnn(
  X = matrix(1:4, ncol = 1), x0 = 0,
  numerator = c(2, 1, 4, 3), denominator = c(1, -1, 2, 1), s = 2,
  return = c("W", "signals", "deletions", "pseudo_values")
)
fit$W
#>           [,1]
#> [1,] 0.5000000
#> [2,] 0.3333333
#> [3,] 0.1666667
#> [4,] 0.0000000
fit$components
#>      numerator denominator
#> [1,]         2         0.5
coef(fit)
#> target_1 
#>        4
```

The weights are $`(1/2,1/3,1/6,0)`$, the components are $`(2,1/2)`$, and
the fitted ratio is four. Here `s` counts observations competing to be
the nearest in each subsample. It does not count the nonzero weights.

For each deletion, keep the signals, target point, metric and order `s`
fixed, and recompute exact weights on the retained observations. The
fitted ratio remains fixed when forming the score. This calculation
requires `s < n` and involves no nuisance refitting.

``` r

fit$deletions[, , 1]
#>          [,1]      [,2]
#> [1,] 2.000000 0.0000000
#> [2,] 2.666667 1.3333333
#> [3,] 1.666667 0.3333333
#> [4,] 1.666667 0.3333333
P <- fit$pseudo_values[, , 1]
P
#>              [,1] [,2]
#> [1,] 2.000000e+00    2
#> [2,] 4.440892e-16   -2
#> [3,] 3.000000e+00    1
#> [4,] 3.000000e+00    1
centered <- sweep(P, 2, colMeans(P), "-")
V <- crossprod(centered) / (4 * 3)
V
#>      [,1] [,2]
#> [1,]  0.5 0.50
#> [2,]  0.5 0.75
stopifnot(isTRUE(all.equal(unname(V), matrix(c(.5, .5, .5, .75), 2))))
```

The pseudo-values are $`P_Y=(2,0,3,3)`$ and $`P_D=(2,-2,1,1)`$. The
first deleted denominator is zero, but all component deletion averages
and pseudo-values are finite. A jackknife of deleted ratios would
perform an undefined division at this deletion; it is a different
statistic.

With $`\hat\theta=4`$, score pseudo-values are
$`Q_j=\hat\theta P_{D,j}-P_{Y,j}`$. Their centered variation gives the
same variance as the paired component covariance:

``` math
\widehat V_\theta=
\frac{V_{YY}+\hat\theta^2V_{DD}-2\hat\theta V_{YD}}
     {\hat\Delta_D^2}.
```

``` r

Q <- fit$theta * P[, 2] - P[, 1]
Q
#> [1]  6 -8  1  1
variance <- sum((Q - mean(Q))^2) / (4 * 3 * fit$components[1, 2]^2)
c(paired_variance = variance,
  omitting_cross_covariance = (V[1, 1] + fit$theta^2 * V[2, 2]) /
    fit$components[1, 2]^2)
#>           paired_variance.denominator omitting_cross_covariance.denominator 
#>                                    34                                    50
fit$theta_vcov
#> [[1]]
#>       theta
#> theta    34
fit$neff
#> [1] 2.571429
fit$theta_omega
#> [[1]]
#>          theta
#> theta 87.42857
stopifnot(abs(variance - 34) < 1e-10,
          abs(fit$neff - 18 / 7) < 1e-10,
          abs(fit$theta_omega[[1]][1, 1] - 612 / 7) < 1e-10)
```

Discarding the cross-covariance gives variance 50 instead of 34 in this
example. Dividing the numerator standard error by the first stage also
omits denominator uncertainty.

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
```

Effective sample size describes rank-weight concentration. It is not a
count of compliers. Physical radius, instrument-arm support, score
dispersion and first-stage strength can vary even when `neff` is the
same. The largest centered pseudo-value-square shares attribute fitted
variance to deletion indices, including rank changes; they differ from
raw weighted score contributions.

## Singular covariance and signed signals

With perfect compliance, treatment equals instrument assignment and the
augmented first-stage signal equals one. Its estimated variance is zero,
while the outcome signal and ratio can still have positive variance.
Positive definiteness of the full two-component covariance is
unnecessary for this scalar calculation.

``` r

A <- c(0, 1, 0, 1)
sig <- iv_augmented_signals(
  Y = c(0, 1, .2, 1.2), D = A, A = A, e = rep(.5, 4),
  m_y0 = rep(0, 4), m_y1 = rep(.5, 4),
  m_d0 = rep(0, 4), m_d1 = rep(1, 4)
)
constant_first_stage <- affine_ratio_dnn(
  X = matrix(1:4), x0 = 0, numerator = sig$phi_y,
  denominator = sig$phi_d, s = 2
)
sig$phi_d
#> [1] 1 1 1 1
constant_first_stage$component_vcov
#> [[1]]
#>             numerator denominator
#> numerator        0.18           0
#> denominator      0.00           0
```

Even in a one-sided compliance model, eligible nonparticipants can have
negative augmented first-stage signals:

``` r

iv_augmented_signals(
  Y = c(0, 1), D = c(0, 1), A = c(1, 1), e = c(.5, .5),
  m_y0 = c(0, 0), m_y1 = c(1, 1), m_d1 = c(.6, .6), one_sided = TRUE
)$phi_d
#> [1] -0.6  1.4
```

The output is $`(-0.6,1.4)`$. Treatment predictions may lie in
$`[0,1]`$, but augmented treatment signals need not. The estimator keeps
signed signals and finite negative denominators. At an exactly zero
aggregate denominator it retains both components and their covariance,
while reporting an undefined ratio and unavailable Wald interval.

The `near_zero_denominator` flag is a numerical cancellation diagnostic;
it does not change the estimator. The first-stage standardized magnitude
is a different, descriptive diagnostic. No universal weak-instrument
cutoff is attached to it. An exactly zero score norm gives a
computational SE of zero and `wald_available = FALSE`, with `NA`
interval bounds.

## Which IV effect is being reported?

Within a reporting profile, consider two equally common adjustment
strata:

``` r

strata <- data.frame(compliance = c(1 / 4, 3 / 4),
                     effect = c(0, 2), propensity = c(1 / 2, 1 / 4))
with(strata, c(
  CLATE = sum(compliance * effect) / sum(compliance),
  average_stratum_LATE = mean(effect),
  residualized_IV = sum(propensity * (1 - propensity) * compliance * effect) /
    sum(propensity * (1 - propensity) * compliance)
))
#>                CLATE average_stratum_LATE      residualized_IV 
#>             1.500000             1.000000             1.384615
```

The three answers are $`3/2`$, $`1`$, and $`18/13`$. Residualized IV
introduces instrument-variance weights. Its target can coincide with
CLATE under additional restrictions, but equality should be checked. The
[local IV-GMM
article](https://JakobJuergens.github.io/RdistNN/articles/local-iv-estimation.md)
documents the package’s separate residualized moment and covariance
workflow.

## What justifies full-sample inference?

For response $`R=Y,D`$, holding candidate nuisance functions fixed gives
the exact conditional mean error

``` math
(e-e_0)\left\{\frac{m_{R,1}-m_{R,1,0}}e+
\frac{m_{R,0}-m_{R,0,0}}{1-e}\right\}.
```

This product removes first-order nuisance mean bias. Individual scores
still depend on observations used for training. The affine theorem also
requires expected prediction convergence over the specified retained
fits, quantitative deletion stability, bounded envelopes and
exceptional-event control, together with oracle moments, local design,
positive score variance, fixed relevance and point-target
undersmoothing. Proof-only deleted nuisance fits are distinct from the
frozen-array deletions used for standard errors.

The bounded
[`clate_demo_data()`](https://JakobJuergens.github.io/RdistNN/reference/clate_demo_data.md)
construction and its designated growing ridge-series training rule
provide a concrete sufficient route developed in the research notes. The
internal ridge adapter is a reference study tool;
[`clate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/clate_dnn.md)
still accepts predictions rather than fitting learners. Array parity
tests verify the algorithm. They do not prove asymptotic learner rates
or causal assumptions. Practical forests, unbounded asset outcomes and
data-selected geometry need their own justification.

## Ratio confidence sets

The experimental
[`ratio_confidence_set()`](https://JakobJuergens.github.io/RdistNN/reference/ratio_confidence_set.md)
helper inverts the candidate moment inequality

``` math
(\widehat\Delta_Y-\theta\widehat\Delta_D)^2
\le z^2(V_{YY}-2\theta V_{YD}+\theta^2V_{DD}).
```

This helper uses the same paired estimator-scale covariance. It can
return bounded intervals, singletons, half-lines, two rays, all real
numbers or empty sets. It remains separate from
[`confint()`](https://rdrr.io/r/stats/confint.html), which reports
pointwise Wald intervals for fitted ratio objects.

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

Finite endpoints are included; infinite endpoints use open notation. An
`unresolved` result has a numerical status and no asserted set geometry,
and differs from a mathematically `empty` result. Materially indefinite
covariance is rejected. Power-of-two scaling and stable roots limit
overflow and cancellation; near-boundary or unrepresentable calculations
remain explicit in the returned record.

With a deterministic nonzero denominator, the helper computes the
interval width directly from the numerator SE. This preserves small
positive uncertainty that can disappear when quadratic coefficients are
expanded. An apparent repeated root of rounded coefficients alone cannot
certify a singleton; an uncertifiable result is `unresolved`. If the
estimator’s component covariance underflows despite nonzero component
variation, it returns unavailable covariance with `numerical_failure`
status. The helper rejects that unavailable input rather than
interpreting it as zero variance.

At a fixed true candidate, inclusion probability is the acceptance
probability of that one candidate’s score test. Simultaneous validity of
all candidate tests is unnecessary for this pointwise statement. The
test still needs its own local-score CLT, consistent variance, nuisance
and bias conditions, including positive variance. Algebraic
zero-variance branches are not automatically covered by that theorem.
Uniform coverage over changing models with increasingly weak first
stages requires an additional sequence-specific argument; the quadratic
calculation does not establish it.
