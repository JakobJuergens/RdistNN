# Paired component and first-stage diagnostics

Paired component and first-stage diagnostics

## Usage

``` r
# S3 method for class 'rdistnn_affine_ratio_result'
diagnostics(object, candidates = "default", ...)

# S3 method for class 'rdistnn_clate_result'
diagnostics(object, candidates = "default", ...)
```

## Arguments

- object:

  An rdistnn_affine_ratio_result or rdistnn_clate_result.

- candidates:

  "default", "all", or candidate ids from a grid-mode fit. Candidate
  diagnostics do not provide candidate-selection inference.

- ...:

  Not used.

## Value

A rdistnn_diagnostics data frame. Component covariance is
estimator-scale. component_covariance reports the paired
numerator-denominator covariance; geometry_rule and tie_rule identify
the frozen geometry and tie convention. Zero denominator variance gives
NA standardized magnitude/correlation. Disabled covariance gives
unavailable uncertainty/concentration fields. physical_radius_90 uses
original-coordinate Euclidean distance; radius_90 uses the fitted
localization metric. Both enclose 90 percent of weight.

## Details

Reports component uncertainty, denominator diagnostics, geometry and
centered score pseudo-value concentration. Instrument and propensity
fields are available only from CLATE observations/predictions. Training
declarations are descriptive. Effective information is not a complier
count.

## Examples

``` r
f <- affine_ratio_dnn(
  X = matrix(1:4, ncol = 1), x0 = 0,
  numerator = c(2, 1, 4, 3), denominator = c(1, -1, 2, 1), s = 2
)
diagnostics(f)
#> <rdistnn_diagnostics>
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
