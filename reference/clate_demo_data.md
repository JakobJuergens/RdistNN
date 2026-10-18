# Bounded endogenous-treatment data for conditional LATE examples

Generates a binary instrument/treatment model with known nuisance
functions. The latent compliance rank also affects untreated outcomes,
so treatment is endogenous conditional on the observed adjustment
variables.

## Usage

``` r
clate_demo_data(
  n = 200L,
  seed = NULL,
  heterogeneous = FALSE,
  keep_latent = FALSE
)
```

## Arguments

- n:

  Positive integer number of observations.

- seed:

  Optional nonnegative integer seed. Like the other package demos,
  supplying a seed resets the R random-number stream.

- heterogeneous:

  Whether effects also vary with an adjustment coordinate. This makes
  CLATE differ from CATE by 3/1024.

- keep_latent:

  Whether to return the simulated latent rank, potential treatments and
  potential outcomes for structural checks.

## Value

A list with data (observations and explicitly named oracle prediction
columns), x0, truth (theta/cate/first_stage functions and core
covariance), configuration and an optional latent audit frame. Truth
functions accept a length-two vector or a matrix with two reporting
columns in the unit square. Latent variables are excluded from the
default adjustment columns.

## Details

X has two uniform reporting coordinates and Z has four independent
uniform adjustment coordinates. With \\S_j=\sin(2\pi Z_j)\\, instrument
propensity is \\1/2+S_4/8\\, always-taker probability is
\\a=1/4+S_1/32\\ and compliance probability is \\c=1/3+S_2/32\\. The
independent latent rank H gives D(0)=1(H\<=a), D(1)=1(H\<=a+c). Outcomes
are \$\$Y(d)=g(X,Z)+\tau(X,Z)d+(H-1/2)/8+\varepsilon,\$\$ with
equiprobable errors +/-1/8, \\g=(S_1+S_4+X_1S_2+X_2S_3)/8\\ and
\\\tau=(X_1+X_2)/8\\. The heterogeneous variant adds \\S_2/16\\ to the
treatment effect. The raw outcome is bounded by one in both models.

Core-model score variance is 13/(48 sqrt(15)); ratio local-limit
variance is 39/(16 sqrt(15)). Those constants do not carry over to the
heterogeneous variant: its sigma/omega metadata are explicitly NA.

## Examples

``` r
z <- clate_demo_data(n = 50, seed = 1)
z$truth$theta(z$x0)
#> [1] 0.125
fit <- clate_dnn(
  data = z$data, x_cols = c("X1", "X2"), y_col = "Y", d_col = "D", a_col = "A",
  e_col = "e_oracle", m_y0_col = "m_y0_oracle", m_y1_col = "m_y1_oracle",
  m_d0_col = "m_d0_oracle", m_d1_col = "m_d1_oracle",
  x0 = z$x0, s = 8, training = "oracle"
)
summary(fit)
#> Components: numerator = reduced_form; denominator = first_stage
#> Ratio statuses: ok=1
#>          estimate std.error statistic p.value conf.low conf.high
#> target_1   0.2433     0.347    0.7011  0.4832  -0.4368    0.9233
```
