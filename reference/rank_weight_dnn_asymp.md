# Rank-weight: DNN asymptotic (full U, WOR)

Rank-weight: DNN asymptotic (full U, WOR)

## Usage

``` r
rank_weight_dnn_asymp(k, s, normalize = TRUE)
```

## Arguments

- k:

  Integer; 1 \<= k \<= s

- s:

  Integer; U-statistic order

- normalize:

  logical; renormalize to sum 1 (default TRUE)

## Value

rank_weight spec (list)

## Numerical limits

Asymptotic intervals use host Beta probabilities. At subnormal scales,
endpoint rounding can erase an interval or change its zero/nonzero
support; R and the Stata companion can disagree. Signed TDNN
combinations can amplify component rounding. These extreme-tail
limitations remain unresolved. Normalization or switching the R/compiled
deletion engine does not recover probability mass already lost during
table construction. Exact and asymptotic weights are different
constructions; changing the method changes the localization rule.

## Examples

``` r
rw <- rank_weight_dnn_asymp(k = 1, s = 3, normalize = TRUE)
rank_weight_eval(rw, c(2L, 1L, 3L, 4L))
#> [1] 0.296875 0.578125 0.109375 0.015625
```
