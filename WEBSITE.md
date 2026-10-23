# RdistNN website copy

Prepared for a software entry or research-page addition on
[jakobjuergens.com](https://jakobjuergens.com/). The text below describes the
production source version `0.1.3`, published on 2026-10-22. The package
documentation is live at its canonical pkgdown URL. The current version-pinned
release is [`v0.1.3`](https://github.com/JakobJuergens/RdistNN/releases/tag/v0.1.3),
with a verified installable source archive and checksum file.

## Software entry

**RdistNN — local treatment effects and nonparametric regression in R**

RdistNN implements distributional nearest-neighbor estimation for conditional
average treatment effects and nonparametric regression. It combines local
estimates with delete-1 jackknife standard errors, tuning comparisons, and
diagnostics for weights and treatment overlap. Supported CLATE and paired-signal
ratios use paired frozen-score covariance and the same reporting/plotting
grammar. Causal nuisance predictions are supplied by the user; valid pointwise
inference requires the documented first-stage, nuisance, localization, stability
and bias conditions. Ratio confidence-set inversion, PLQR and local IV-GMM
remain experimental.

[Source and installation](https://github.com/JakobJuergens/RdistNN) ·
[Documentation](https://JakobJuergens.github.io/RdistNN/) ·
[Getting started](https://JakobJuergens.github.io/RdistNN/articles/getting-started.html) ·
[Accompanying paper](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z)

## 0.1.3 production-update announcement draft

RdistNN 0.1.3 averages localization weights over observations at exactly equal
represented distance under the selected metric. The rule covers distinct
covariate profiles at the same distance, preserves each group's assigned mass,
and recomputes survivor weights after every jackknife deletion. Shared R/C++
helpers apply it to fitting, prediction, candidates, diagnostics and weight
plots. Refit saved objects from earlier tie policies.

The paper assumes a continuous distance distribution, under which the old and
new rules agree almost surely. Inference with persistent distance ties remains
unestablished. Numerical symmetry does not establish the underlying inference
conditions. Shared reporting now gives bounded, actionable errors and warnings
while preserving unexpected failures and the caller's warning policy.

The [weight guide](https://JakobJuergens.github.io/RdistNN/articles/dnn-weights.html#numerical-limits)
describes remaining numerical limits, including asymptotic-tail rounding and
independently fitted Mahalanobis geometry.

Task-oriented guides, runnable vignettes and focused function help support
conditional-effect and regression workflows. Seeded demos and simulations
preserve caller RNG state; native calculations use two threads by default with
a configurable opt-in. Version 0.1.3 is available from production
`master` and the version-pinned `v0.1.3` release.

[Source updates](https://github.com/JakobJuergens/RdistNN/blob/master/NEWS.md) ·
[Tagged releases](https://github.com/JakobJuergens/RdistNN/releases) ·
[Getting started](https://JakobJuergens.github.io/RdistNN/articles/getting-started.html)

Use the [verified release URL](https://github.com/JakobJuergens/RdistNN/releases/tag/v0.1.3)
for version-pinned downloads.
This repository copy is prepared text for the personal website; package-site
deployment and personal-site publication are separate steps.
