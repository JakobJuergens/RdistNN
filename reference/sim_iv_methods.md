# List available IV simulation methods

Returns the names of all estimator methods registered for use with
[`sim_iv_one_rep()`](https://JakobJuergens.github.io/RdistNN/reference/sim_iv_one_rep.md)
and
[`sim_iv_run()`](https://JakobJuergens.github.io/RdistNN/reference/sim_iv_run.md).

## Usage

``` r
sim_iv_methods()
```

## Value

Character vector of registered IV simulation method names.

## Examples

``` r
sim_iv_methods()
#> [1] "dnn_iv_oracle"       "dnn_iv_oracle_mahal" "dnn_iv_fullsample"  
#> [4] "dnn_iv_crossfit"     "iv_oracle"           "2sls_global"        
#> [7] "2sls_interact"       "local_2sls"         
```
