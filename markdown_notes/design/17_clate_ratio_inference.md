# Conditional LATE and paired-ratio inference

`clate_dnn()` constructs binary-IV augmented signals and delegates localization
and uncertainty to `affine_ratio_dnn()`. Core supplied-prediction CLATE and
ordinary supplied-signal ratios are supported workflows under the contracts
below. The optional `ratio_confidence_set()` inversion API remains experimental.
Executable user examples are in [Conditional LATE](../../vignettes/conditional-late.Rmd)
and [Paired frozen-score inference](../../vignettes/conditional-late-inference.Rmd).
This document records the numerical and scientific contracts for maintenance.

## Target and signal construction

Let `A` be a binary instrument, `D` binary treatment and `V = (X,Z)` the full
pre-instrument adjustment state. Reporting weights use `X`. For `R = Y,D`,
write `m_Ra(V) = E[R | A=a,V]` and `e(V) = P(A=1 | V)`. The augmented signal is

$$
\phi_R=m_{R1}-m_{R0}+\frac A e(R-m_{R1})
                         -\frac{1-A}{1-e}(R-m_{R0}).
$$

The five supplied nuisance predictions are indexed by instrument arm and use
the same row order as the observations. Both signals use the same instrument
propensity. `iv_augmented_signals()` retains their construction and provenance.
Unique row IDs document alignment; one shared ID vector cannot detect an
independently permuted prediction array. Align external predictions explicitly.

Propensities are clipped according to the documented `clip` control, retaining
raw/used values and clipping flags. Treatment-mean predictions must be in
`[0,1]`; outcome predictions must be finite. `one_sided = TRUE` selects
`D(0)=0`, permits an omitted zero-arm treatment prediction and checks `D <= A`.
The observed check does not establish the counterfactual restriction.

The standardized contrasts integrate over the same law of `Z | X=x`:

$$
\Delta_R(x)=E[m_{R1}(V)-m_{R0}(V)\mid X=x],\qquad
\theta(x)=\Delta_Y(x)/\Delta_D(x).
$$

Consistency, exclusion, conditional instrument independence given `V`,
monotonicity, instrument overlap and positive local compliance give this ratio
its conditional complier-effect interpretation. Without those assumptions it
is an observed-data Wald functional wherever its denominator is nonzero.
It generally differs from population CATE and from residualized IV-GMM.

## Localization and frozen state

The solution of the affine score `theta * phi_d - phi_y` is a ratio of the two
weighted signal averages. At each target, observations with exactly equal
represented selected-metric distance share the mean of their assigned rank
weights, using the same groups for both signals. Each deletion removes the
actual observation before averaging its remaining group's `n-1` weights.
This preserves paired covariance and aligned permutation invariance; groups
use localization distances, not instrument, treatment or nuisance values.
The paper assumes a continuous distance distribution, under which the earlier
and current rules agree almost surely. The warning records that inference for
persistent distance ties is not established.
Never average observation-level signal ratios,
floor the aggregate denominator or clip negative augmented first stages.

Exact/asymptotic single-scale DNN, general `1 <= k <= s`, normalization,
distance, tuning and optional extras follow the common estimator grammar.
Exact DNN and `k = 1L` remain defaults. Covariance requires `s < n`; point-only
fits can use `s = n`. Two-scale TDNN is not a CATE/CLATE option.

Full and deleted calculations use the same order and weight rule. Targets,
signals, distance state and full-sample tie-broken geometry remain frozen;
deletion removes the row and updates ranks among retained observations.
`predict()` changes reporting targets and preserves the selected weight
controls. It does not train nuisances or predict outcomes for new households.

Grid sensitivity reuses the [common grid engine](14_integrated_tuning_grid_reuse.md).
Candidate views rebuild paired components, covariance and requested extras
from retained signals/ranks while preserving IV state. They are scalar-style
views of a specified candidate, not automatic order selection.

## Paired covariance and scales

Let `L_i = (L_Yi,L_Di)` be the pair of deleted component averages, with the
same deleted row in both signals, and let `L_bar` be their average. The
estimator-scale component covariance is

$$
V=\frac{n-1}{n}\sum_i(L_i-\bar L)(L_i-\bar L)' .
$$

At the full fitted ratio the scalar variance is

$$
V_\theta=\{V_{YY}+\widehat\theta^2V_{DD}
                  -2\widehat\theta V_{YD}\}/\widehat\Delta_D^2.
$$

Keep the cross-covariance. Compute this primarily as the squared norm of the
centered deleted-score combination divided by the full denominator, retaining
the quadratic formula as an independent reference. This avoids cancellation
of large terms. Deleted denominators may be zero without making this paired
calculation undefined. Deleted-ratio covariance, retrained-pipeline jackknife
and a weighted score-cross-product meat are different calculations.

| Field | Meaning |
| --- | --- |
| `components` | Target-by-two matrix, numerator then denominator |
| `component_vcov[[j]]` | Paired 2-by-2 estimator covariance |
| `theta_vcov[[j]]` | 1-by-1 ratio variance, already on the estimator scale |
| `theta_omega[[j]]` | `neff[j] * theta_vcov[[j]]` |
| `neff[j]` | Inverse squared-weight sum, not a complier count |
| `se[j]` | Square root of ratio variance; no further division by `neff` |

Covariances are marginal across targets. Optional deletions and pseudo-values
have dimensions `n x 2 x m`. With covariance disabled the covariance fields
are `NULL`. Summary coefficient columns follow CATE; component and status
tables remain separate. The [result architecture](13_estimator_and_result_architecture.md)
specifies generic dispatch and common plotting semantics.

## Numerical statuses and confidence sets

Retain every target, including a finite negative denominator or an undefined
ratio. `ratio_status` distinguishes `ok`, `zero_denominator`,
`component_overflow` and `ratio_overflow`. Near-zero and negative-denominator
flags do not replace or filter estimates. `covariance_status` separately
identifies unavailable, failed or zero-score covariance calculations.
Nonfinite computed component averages are returned as `NA` with
`component_overflow`; finite components and the original retained signals stay
available. Status and covariance are computed before this reporting conversion.

Nonzero component variation that underflows its variance is a numerical
failure, not a deterministic covariance. Zero score variance has computational
SE zero but unavailable Wald statistics and bounds. `wald_available` controls
rowwise inference; an invalid target must not suppress another target's
available inference or manufacture a zero-width interval.

`ratio_confidence_set()` separately inverts the fixed-candidate paired moment
inequality. It returns interval unions, including rays, singleton, empty and
all-real sets, and a distinct unresolved numerical status. Preserve covariance
uncertainty through coefficient formation and verify roots against the input
components/covariance; rounded polynomial coefficients cannot certify a false
singleton. This helper does not replace `confint()` or establish uniform
weak-instrument coverage.

Covariance validation uses the correlation calculation away from its PSD
boundary and the exact sign of the represented two-by-two determinant near
absolute correlation one. Rounded square roots can otherwise hide a negative
determinant or reject an exactly singular PSD matrix. At `tol = 0` a negative
represented determinant is rejected; a small negative determinant within the
relative tolerance gives `near_psd_boundary` and unresolved intervals. The
near-boundary sign check reuses the shared exact sum-of-products arithmetic.
`R/utils-exact-products.R` owns that arithmetic for both this check and the
PLQR comparison adapter; ratio inference must not copy the expansion logic.

Diagnostics expose both components and their covariance, first-stage variation,
instrument support, clipping, localization and centered score concentration.
The geometry and tie rules remain explicit. No automatic weak-instrument
cutoff or target filtering is implied.

## Statistical scope

See the [paper](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z)
for pointwise inference assumptions: local/oracle regularity, nonzero population
first stage, bias control and learner accuracy/deletion stability. Supplied
training labels and numerical tests do not certify them. The covariance freezes
predictions; theoretical one-/two-deleted nuisance fits are proof devices.

At a valid fixed target the limit variance is
`Var(theta * phi_d - phi_y | X=x) / Delta_D(x)^2`. General-k/asymptotic
computational support does not extend an exact-k=1 theorem. Automatic nuisance
fitting, adaptive selection, simultaneous bands, multiple instruments and
weak-IV-uniform inference require separate justification.

## Independent numerical references

For `X = matrix(1:4)`, target zero, `s = 2`, numerator `(2,1,4,3)` and
denominator `(1,-1,2,1)`, exact weights are `(1/2,1/3,1/6,0)`. Components are
`(2,1/2)`, the ratio is `4`, paired covariance is
`rbind(c(1/2,1/2), c(1/2,3/4))`, ratio variance is `34`, `neff = 18/7`, and
`theta_omega = 612/7`. The first deleted denominator is zero.

Asymptotic weights for the same arrays are `(7,5,3,1)/16`; components are
`(17/8,9/16)` and the ratio is `34/9`. General-k exact weights are independently
checked by complete-subset enumeration. Under perfect compliance with
`phi_d = 1`, estimates and uncertainty reduce to CATE.

`clate_demo_data()` supplies a bounded latent-compliance benchmark with known
nuisances. Its core first stage is `1/3`; CLATE and CATE are `(X1+X2)/8`.
The heterogeneous option adds `sin(2*pi*Z2)/16` to the treatment effect and
gives CLATE minus CATE equal to `3/1024`. Latent variables are optional audit
output, not default learner covariates. Core oracle variance constants must
not be reused for the heterogeneous model.

The internal ridge-series reference adapter fixes its dictionary, ambient
normalizer, penalty and clipping across deletions. It is tested against dense
normal equations and explicit one-/two-deleted fits, including empty cells.
These checks validate computation, not arbitrary learner rates. Small generated
examples are not repeated-sampling coverage studies or manuscript evidence.

See [the test-suite map](../../tests/README.md) for independent fixtures,
engine parity, numerical geometry, interfaces, serialization and rendering.
