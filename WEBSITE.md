# RdistNN website copy

Prepared for a software entry or research-page addition on
[jakobjuergens.com](https://jakobjuergens.com/). The text below
describes the production source version `0.1.4`, published on
2026-10-23. The package documentation is live at its canonical pkgdown
URL. The current version-pinned release is
[`v0.1.4`](https://github.com/JakobJuergens/RdistNN/releases/tag/v0.1.4),
with a verified installable source archive and checksum file.

## Software entry

**RdistNN — local treatment effects and nonparametric regression in R**

RdistNN implements distributional nearest-neighbor estimation for
conditional average treatment effects and nonparametric regression. It
combines local estimates with delete-1 jackknife standard errors, tuning
comparisons, and diagnostics for weights and treatment overlap.
Supported CLATE and paired-signal ratios use paired frozen-score
covariance and the same reporting/plotting grammar. Causal nuisance
predictions are supplied by the user; valid pointwise inference requires
the documented first-stage, nuisance, localization, stability and bias
conditions. Ratio confidence-set inversion, PLQR and local IV-GMM remain
experimental.

[Source and installation](https://github.com/JakobJuergens/RdistNN) ·
[Documentation](https://JakobJuergens.github.io/RdistNN/) · [Getting
started](https://JakobJuergens.github.io/RdistNN/articles/getting-started.html)
· [Accompanying
paper](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z)

## 0.1.4 production-update announcement draft

RdistNN 0.1.4 strengthens fitted-object validation, captures extension
configuration at construction and reuses retained localization groups
across supported workflows. Current reporting gives bounded, actionable
conditions while preserving unexpected failures and the caller’s warning
policy. The release also improves numerical behavior for very small
treatment-arm mass and paired cross-covariances, and preserves RNG state
in deterministic native calculations. Earlier saved objects must be
refitted.

Exact represented distance ties continue to receive averaged
localization weights, including after jackknife deletion. The paper
assumes a continuous distance distribution; inference with persistent
distance ties remains unestablished.

The [weight
guide](https://JakobJuergens.github.io/RdistNN/articles/dnn-weights.html#numerical-limits)
describes remaining numerical limits, including asymptotic-tail rounding
and independently fitted Mahalanobis geometry.

Task-oriented guides, runnable vignettes and focused function help
support conditional-effect and regression workflows. Seeded demos and
simulations preserve caller RNG state; native calculations use two
threads by default with a configurable opt-in. Version 0.1.4 is
available from production `master` and the version-pinned `v0.1.4`
release.

[Source
updates](https://github.com/JakobJuergens/RdistNN/blob/master/NEWS.md) ·
[Tagged releases](https://github.com/JakobJuergens/RdistNN/releases) ·
[Getting
started](https://JakobJuergens.github.io/RdistNN/articles/getting-started.html)

Use the [verified release
URL](https://github.com/JakobJuergens/RdistNN/releases/tag/v0.1.4) for
version-pinned downloads. This repository copy is prepared text for the
personal website; package-site deployment and personal-site publication
are separate steps.
