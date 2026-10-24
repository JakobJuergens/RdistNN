# Package scope and extension boundaries

RdistNN provides modular distance, rank-weight, score and jackknife primitives
with high-level estimation workflows. Supported workflows include CATE,
nonparametric regression, supplied-prediction conditional LATE and ordinary
paired affine ratios. Their documented inference requirements remain part of
their use. Optional ratio confidence-set inversion, PLQR and local IV-GMM
remain experimental; the package-wide lifecycle remains experimental.

The [design reference](design/README.md) specifies maintained contracts.
[CONTRIBUTING.md](../CONTRIBUTING.md) owns development commands, and
[RELEASING.md](../RELEASING.md) owns release checks. The public notation hub is
`vignettes/notation-and-conventions.Rmd`.

## Current capabilities

- High-level localization averages assigned rank weights within exact represented
  selected-metric distance groups, separately for every target and deletion's
  surviving observations. Cached groups and the selected policy are validated
  across fitting and reuse; legacy fits require refitting. The paper assumes a
  continuous distance distribution, under which old and new policies agree almost
  surely. Persistent distance-tie inference and a uniform MSE improvement are
  not established. Version 0.1.3 uses this exact-distance policy; released
  0.1.2 and scientific dependency pins remain immutable.

- Matrix/vector and data-frame interfaces resolve to common numeric inputs.
- Exact/asymptotic DNN supports shared tuning sensitivity, diagnostics and
  presentation. CATE/CLATE use single-scale DNN; TDNN remains NP-regression
  functionality.
- NP DNN and TDNN expose conditional-mean summaries, coefficients, pointwise
  intervals, full joint target covariance and prediction from retained training
  state, including serialized fits without optional weight/outcome extras.
- CATE accepts AIPW nuisance inputs or supplied signals. With fitted uncertainty
  and an active jackknife, `vcov()` supplies joint target covariance from
  retained model state, including serialized fits without optional extras.
  Common deletions follow the fitted jackknife: `fast` freezes geometry and
  `loo` refits deletion geometry.
  Covariance supplies the arithmetic for prespecified contrasts; valid contrast
  inference requires the corresponding joint assumptions. See the
  [covariance reference](https://JakobJuergens.github.io/RdistNN/reference/vcov.rdistnn_cate_result.html).
- CLATE constructs paired augmented binary-IV signals from supplied nuisance
  predictions, retaining components, estimator-scale covariance and numerical
  statuses. Its summaries, plots and candidate views use the shared interfaces.
- PLQR evaluates moments at supplied coefficients. Coefficient-level inference
  needs the documented density-backed sandwich; componentwise moment SEs do
  not supply coefficient inference.
- Local IV-GMM solves a distinct moment system with its own covariance and
  convergence diagnostics. It is not an alias for the augmented CLATE ratio.
  Targets use independent optimizer starts and retain their fitted geometry
  for later weight reconstruction.
- The installed affine replication bundle pins scientific source, production
  configurations and accepted evidence independently of the general package
  estimators. Its manifest identifies the exact scientific snapshot.

## Extension boundaries

A unified score-localize-solve-report conditional Z-estimation front end is not
implemented. Current wrappers take externally fitted nuisance predictions;
keeping nuisance training outside estimation is a deliberate design boundary.
The bounded CLATE generator and independent supplied-signal references provide
numerical verification. Automatic nuisance learning is outside the estimator
interface.

Adaptive candidate selection, post-selection inference, simultaneous bands,
weak-instrument-uniform guarantees and new empirical/production campaigns
require separate scientific and interface justification. Successful tests,
training labels and computational weighting choices do not establish them.
Maintain a shared usage and presentation grammar where quantities agree in
meaning; preserve method-specific scores, covariance and scientific limits.

## Compensated accumulation

Asymptotic Beta-tail construction has a known numerical limitation at
subnormal scales, including R/Mata support differences and lost masses after
averaging or signed TDNN blending. Normalization and native deletion engines
do not resolve it. See the [weight guide](../vignettes/dnn-weights.Rmd#numerical-limits)
and [numerical reference](design/19_weighted_accumulation.md); this remains an
open numerical boundary rather than a relaxed parity requirement.

Rank-weight normalization, weighted signal means and frozen-rank deletions
use shared compensated arithmetic across native and R engines and the Stata
companion. Prefix/suffix calculations retain correction state until the final
merge. The finite-support jackknife adapter reuses the same primitive;
explicit deletion references remain independently testable.

Support-aware centering, frozen geometry, signed TDNN weights and paired
covariance retain their existing definitions. R/Mata reverse sweeps use bounded
merge blocks rather than full-length vector temporaries for every addition.
These finite double-precision safeguards do not promise arbitrary-precision
results or recover information already lost in multiplication. See
[weighted accumulation](design/19_weighted_accumulation.md) for ownership,
workspace, engine controls and numerical verification.
