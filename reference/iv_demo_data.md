# Minimal IV-GMM demo data

Generates a small synthetic local IV-GMM dataset for experimenting with
oracle doubly-robust IV moments. The output is a thin, test-friendly
wrapper around
[`sim_iv_data()`](https://JakobJuergens.github.io/RdistNN/reference/sim_iv_data.md)
with shorter design names.

## Usage

``` r
iv_demo_data(
  n = 200,
  p = 2,
  q = 2,
  design = c("homog", "heterog", "weak"),
  sigma = 0.5,
  clip = 0.001,
  seed = NULL,
  x0 = NULL,
  inject = NULL,
  inject_col = "Y",
  idx_inject = 3L
)
```

## Arguments

- n:

  Integer sample size.

- p:

  Integer feature dimension, at least 2.

- q:

  Integer number of instruments, between 1 and 3.

- design:

  Character; one of `"homog"`, `"heterog"`, or `"weak"`.

- sigma:

  Noise standard deviation in `Y`.

- clip:

  Clipping bound in `(0, 0.5)` for binary instrument probabilities.

- seed:

  Optional integer seed, used locally without changing the caller's RNG
  state. `NULL` uses the current stream.

- x0:

  Optional numeric length-`p` target point.

- inject:

  Optional, one of `NULL`, `"na"`, or `"inf"`; if provided, injects a
  single non-finite value into `inject_col` at `idx_inject`.

- inject_col:

  Character injection target.

- idx_inject:

  Integer row to corrupt.

## Value

A list with the same structure as
[`sim_iv_data()`](https://JakobJuergens.github.io/RdistNN/reference/sim_iv_data.md).

## Examples

``` r
demo <- iv_demo_data(n = 120, p = 2, q = 2, design = "homog", seed = 1)
str(demo$data)
#> 'data.frame':    120 obs. of  11 variables:
#>  $ X1          : num  -0.626 0.184 -0.836 1.595 0.33 ...
#>  $ X2          : num  -0.506 1.343 -0.215 -0.18 -0.1 ...
#>  $ Y           : num  1.516 -0.111 0.494 1.608 1.316 ...
#>  $ D           : int  1 0 1 0 0 1 1 1 0 0 ...
#>  $ Z1          : num  1 0 1 0 1 0 1 1 0 0 ...
#>  $ Z2          : num  1 1 0 0 0 1 1 1 0 1 ...
#>  $ piZ1_oracle : num  0.44 0.396 0.404 0.611 0.498 ...
#>  $ piZ2_oracle : num  0.489 0.445 0.453 0.658 0.547 ...
#>  $ muY_oracle  : num  0.614 0.762 0.216 3.43 1.771 ...
#>  $ muD_oracle  : num  0.417 0.373 0.375 0.629 0.489 ...
#>  $ theta_oracle: num  -0.0291 0.6707 -0.1485 1.1838 0.8003 ...
demo$theta_x0
#> [1] 1.211738
```
