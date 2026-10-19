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

- Matrix/vector and data-frame interfaces resolve to common numeric inputs.
- Exact/asymptotic DNN supports shared tuning sensitivity, diagnostics and
  presentation. CATE/CLATE use single-scale DNN; TDNN remains NP-regression
  functionality.
- NP DNN and TDNN expose conditional-mean summaries, coefficients, pointwise
  intervals, full joint target covariance and prediction from retained training
  state, including serialized fits without optional weight/outcome extras.
- CLATE constructs paired augmented binary-IV signals from supplied nuisance
  predictions, retaining components, estimator-scale covariance and numerical
  statuses. Its summaries, plots and candidate views use the shared interfaces.
- PLQR evaluates moments at supplied coefficients. Coefficient-level inference
  needs the documented density-backed sandwich; componentwise moment SEs do
  not supply coefficient inference.
- Local IV-GMM solves a distinct moment system with its own covariance and
  convergence diagnostics. It is not an alias for the augmented CLATE ratio.
- The installed affine replication bundle pins scientific source, production
  configurations and accepted evidence independently of the general package
  estimators. Its manifest identifies the exact scientific snapshot.

## Extension boundaries

A unified score-localize-solve-report conditional Z-estimation front end is not
implemented. Current wrappers take externally fitted nuisance predictions;
keeping nuisance training outside estimation is a deliberate design boundary.
The bounded CLATE generator and internal ridge reference adapter provide
numerical verification, not a general automatic learner interface.

Adaptive candidate selection, post-selection inference, simultaneous bands,
weak-instrument-uniform guarantees and new empirical/production campaigns
require separate scientific and interface justification. Successful tests,
training labels and computational weighting choices do not establish them.
Maintain a shared usage and presentation grammar where quantities agree in
meaning; preserve method-specific scores, covariance and scientific limits.

## Planned compensated accumulation

Compensated sensitive weighted and prefix/suffix sums are planned, not
implemented. Evaluate a shared compensated accumulator for rank-weight mass,
signed TDNN weighted estimates and the forward/backward deletion sums in
`R/jackknife-rank-full.R` and `src/jackknife_rank_full.cpp`, alongside the
corresponding Stata weight and rank-deletion helpers. Preserve support-aware
centering, frozen geometry, normalization choices and paired covariance.

Use constant correction state per active sum, reusing existing prefix/suffix
buffers rather than adding another deletion array. Compare candidate schemes
against independent high-accuracy and explicit-deletion references on long
sums, mixed magnitudes/signs, near cancellation, constant signals and unused
extremes. Run the same fixtures through C++, R and Mata, including returned
deletions, covariance, candidates and saved prediction. Check exact discrete
decisions separately from justified numeric tolerances.

Adoption requires measured runtime and peak-memory results for representative
and adversarial workloads. In particular, replacing vectorized R/Mata
`cumsum`/`runningsum` or dot products with interpreted scalar loops may outweigh
the numerical gain. Keep this work separate from the stable TDNN coefficient,
Beta probability-increment and exceptional L2 norm implementations; those
changes do not establish compensated-sum coverage. The
[compiled-code policy](design/10_rcpp_and_compiled_code_policy.md) defines the
cross-engine numerical contract.
