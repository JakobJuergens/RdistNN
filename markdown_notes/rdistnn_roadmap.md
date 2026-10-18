# Package scope and extension boundaries

RdistNN provides modular distance, rank-weight, score and jackknife primitives
with high-level estimation workflows. The initial public release centers on
CATE and nonparametric regression. The development version additionally
supports conditional LATE and paired affine ratios. CLATE, PLQR and local
IV-GMM are experimental; their documented inference requirements remain part
of their use.

The [design reference](design/README.md) specifies maintained contracts.
[CONTRIBUTING.md](../CONTRIBUTING.md) owns development commands, and
[RELEASING.md](../RELEASING.md) owns release checks. The public notation hub is
`vignettes/notation-and-conventions.Rmd`.

## Current capabilities

- Matrix/vector and data-frame interfaces resolve to common numeric inputs.
- Exact/asymptotic DNN supports shared tuning sensitivity, diagnostics and
  presentation. CATE/CLATE use single-scale DNN; TDNN remains NP-regression
  functionality.
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
