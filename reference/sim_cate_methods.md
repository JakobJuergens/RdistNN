# List available CATE simulation methods

Returns the names of all estimator methods registered for use with
[`sim_cate_one_rep()`](https://JakobJuergens.github.io/RdistNN/reference/sim_cate_one_rep.md)
and
[`sim_cate_run()`](https://JakobJuergens.github.io/RdistNN/reference/sim_cate_run.md).

## Usage

``` r
sim_cate_methods()
```

## Value

Character vector of registered method names.

## Details

|  |  |
|----|----|
| Name | Description |
| `dnn_oracle` | Proposed DNN Z-estimator, L2 distance, oracle nuisances |
| `dnn_oracle_mahal` | Proposed DNN Z-estimator, Mahalanobis distance, oracle nuisances (`dgp5_correlated`) |
| `dnn_fullsample` | Proposed DNN Z-estimator, L2 distance, full-sample estimated AIPW nuisances |
| `dnn_crossfit` | Proposed DNN Z-estimator, L2 distance, K-fold cross-fitted AIPW nuisances |
| `local_constant` | Nadaraya-Watson smoother on the oracle AIPW score |
| `ols_interact` | OLS with D\*X interactions, delta-method SE at x0 |
| `causal_forest` | GRF causal forest, oracle nuisances (requires **grf**) |
| `llf` | GRF causal forest with local-linear corrections (requires **grf**) |

## See also

[`sim_cate_run()`](https://JakobJuergens.github.io/RdistNN/reference/sim_cate_run.md),
[`sim_cate_one_rep()`](https://JakobJuergens.github.io/RdistNN/reference/sim_cate_one_rep.md)

## Examples

``` r
sim_cate_methods()
#> [1] "dnn_oracle"       "dnn_oracle_mahal" "dnn_fullsample"   "dnn_crossfit"    
#> [5] "local_constant"   "ols_interact"     "causal_forest"    "llf"             
```
