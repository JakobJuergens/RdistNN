# Doubly robust IV GMM model

Creates the scalar local IV model used by the GMM IV simulation. The
observation-level moment is \$\$ (Z_j - \pi_j(X))\\Y - \mu_Y(X) - (D -
\mu_D(X))\theta\\, \$\$ evaluated for each instrument `j`.

## Usage

``` r
gmm_model_dr_iv(y_col = "Y", d_col = "D", z_cols, pi_cols, muY_col, muD_col)
```

## Arguments

- y_col, d_col:

  Column names or indices for outcome and endogenous treatment.

- z_cols:

  Column names or indices for instruments.

- pi_cols:

  Column names or indices for \\E\[Z_j \mid X\]\\.

- muY_col:

  Column name or index for \\E\[Y \mid X\]\\.

- muD_col:

  Column name or index for \\E\[D \mid X\]\\.

## Value

A `gmm_model` object for use with
[`gmm_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/gmm_dnn.md).

## Details

No nuisance models are fit internally. All nuisance values must be
supplied as columns in `data` or `model_data`.

## Examples

``` r
dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 80, kx = 2, q = 2, seed = 1)
model <- gmm_model_dr_iv(
  z_cols = c("Z1", "Z2"),
  pi_cols = c("piZ1_oracle", "piZ2_oracle"),
  muY_col = "muY_oracle",
  muD_col = "muD_oracle"
)
is_gmm_model(model)
#> [1] TRUE
```
