# Reporting and prediction for paired affine-ratio results

These methods apply to affine_ratio_result and its clate_result
subclass. Intervals are pointwise normal-reference intervals with
frozen-score standard errors. A zero or numerically unavailable ratio
variance gives NA interval bounds, not a zero-width confidence interval.

## Usage

``` r
# S3 method for class 'affine_ratio_result'
coef(object, ...)

# S3 method for class 'affine_ratio_result'
vcov(object, target = NULL, ...)

# S3 method for class 'affine_ratio_result'
confint(object, parm = NULL, level = 0.95, ...)

# S3 method for class 'affine_ratio_result'
summary(object, level = 0.95, ...)

# S3 method for class 'affine_ratio_result'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'summary.affine_ratio_result'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'affine_ratio_result'
predict(
  object,
  newdata = NULL,
  X0 = NULL,
  x0 = NULL,
  se.fit = FALSE,
  interval = c("none", "confidence"),
  level = 0.95,
  ...
)
```

## Arguments

- object, x:

  A fitted affine_ratio_result or clate_result.

- ...:

  Not used.

- target:

  Optional single target index for vcov.

- parm:

  Optional numeric target indices for confint.

- level:

  Confidence level in (0,1).

- digits:

  Printing precision.

- newdata:

  Optional reporting-target matrix/data frame. For a data-frame fit,
  only the retained reporting columns are extracted.

- X0, x0:

  Target matrix or single target vector, alternatives to newdata.

- se.fit:

  Whether prediction returns standard errors.

- interval:

  Either "none" or "confidence".

## Value

coef returns a named vector; vcov returns one marginal matrix or a list
of marginal matrices; confint returns a two-column matrix. summary
returns a classed list with coefficients, level, has_se and meta, plus
component and per-target status tables. The coefficient table follows
the CATE estimate/std.error/statistic/p.value/conf.low/conf.high schema.
predict returns a vector, a fit/se.fit list or a fit/lwr/upr matrix.
Printing returns its input invisibly.

## Details

Covariance matrices are on the estimator scale. Prediction changes
reporting targets only, holding the observations, signals and order
fixed. It never trains nuisance models. A supplied-training label does
not verify the statistical or causal assumptions underlying inference.

## Examples

``` r
f <- affine_ratio_dnn(
  X = matrix(1:4, ncol = 1), x0 = 0,
  numerator = c(2, 1, 4, 3), denominator = c(1, -1, 2, 1), s = 2
)
coef(f)
#> target_1 
#>        4 
vcov(f)
#>       theta
#> theta    34
confint(f)
#>           conf.low conf.high
#> target_1 -7.428456  15.42846
summary(f)
#> Components: numerator = numerator; denominator = denominator
#> Ratio statuses: ok=1
#>          estimate std.error statistic p.value conf.low conf.high
#> target_1        4     5.831     0.686  0.4927   -7.428     15.43
predict(f, x0 = 5, se.fit = TRUE)
#> $fit
#> [1] 3
#> 
#> $se.fit
#> [1] 1.414214
#> 
```
