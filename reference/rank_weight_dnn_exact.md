# Rank-weight: DNN exact (full U, WOR)

Rank-weight: DNN exact (full U, WOR)

## Usage

``` r
rank_weight_dnn_exact(k, s)
```

## Arguments

- k:

  Integer; 1 \<= k \<= s

- s:

  Integer; U-statistic order

## Value

rank_weight spec (list)

## Examples

``` r
rw <- rank_weight_dnn_exact(k = 1, s = 3)
rank_weight_eval(rw, c(2L, 1L, 3L, 4L))
#> [1] 0.25 0.75 0.00 0.00
```
