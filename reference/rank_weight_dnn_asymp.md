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

## Examples

``` r
rw <- rank_weight_dnn_asymp(k = 1, s = 3, normalize = TRUE)
rank_weight_eval(rw, c(2L, 1L, 3L, 4L))
#> [1] 0.296875 0.578125 0.109375 0.015625
```
