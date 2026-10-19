# Construct augmented reduced-form and first-stage signals

Constructs two aligned signals from externally supplied nuisance
predictions. This interface fits no nuisance models.

## Usage

``` r
iv_augmented_signals(
  Y,
  D,
  A,
  e,
  m_y0,
  m_y1,
  m_d0 = NULL,
  m_d1,
  clip = 1e-08,
  one_sided = FALSE,
  row_ids = NULL,
  training = "supplied"
)
```

## Arguments

- Y:

  Finite numeric outcome vector.

- D:

  Numeric binary treatment vector (0/1).

- A:

  Numeric binary instrument vector (0/1).

- e:

  Instrument propensity predictions in \\\[0,1\]\\, before clipping.

- m_y0, m_y1:

  Outcome predictions conditional on instrument arm 0 and 1.

- m_d0, m_d1:

  Treatment predictions conditional on instrument arm 0 and 1, in
  \\\[0,1\]\\. m_d0 may be omitted only with one_sided = TRUE.

- clip:

  Propensity clipping threshold in (0,0.5), with `1 - clip < 1` in
  double precision.

- one_sided:

  Explicit structural restriction that treatment is impossible under
  instrument 0. Requires D \<= A and zero m_d0 if supplied.

- row_ids:

  Optional unique, nonmissing character or numeric row identities. These
  record alignment; a single ID vector cannot detect independently
  misordered prediction arrays.

- training:

  Descriptive declaration: supplied, oracle, full_sample, oob or
  cross_fit. It is not verified from predictions.

## Value

An iv_signals list with phi_y, phi_d, row_ids, propensity (raw, used,
clipped), nuisance (m_y0, m_y1, m_d0, m_d1) and meta.

## Details

For either response \\R\\, the signal is
\$\$\phi_R=m\_{R1}-m\_{R0}+A(R-m\_{R1})/e -(1-A)(R-m\_{R0})/(1-e).\$\$
Both responses use the same clipped propensity. Arm subscripts refer to
instrument assignment, not treatment received. Treatment predictions are
validated, not silently clipped. Outcomes and outcome predictions are
not capped. Finite unbounded-outcome computations do not certify the
bounded learner conditions of the affine inference theorem.

Individual first-stage signals can be negative even under one-sided
compliance. They remain unchanged. The CLATE equation averages both
signals before dividing; do not average observation-level signal ratios.

## Examples

``` r
z <- iv_augmented_signals(
  Y = c(1, 2, 0, 3), D = c(0, 1, 0, 1), A = c(0, 1, 1, 1),
  e = rep(.5, 4), m_y0 = rep(1, 4), m_y1 = rep(2, 4),
  m_d1 = rep(.6, 4), one_sided = TRUE
)
z$phi_y
#> [1]  1  1 -3  3
z$phi_d
#> [1]  0.6  1.4 -0.6  1.4
```
