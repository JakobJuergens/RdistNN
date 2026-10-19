# RdistNN website copy

Prepared for a software entry or research-page addition on
[jakobjuergens.com](https://jakobjuergens.com/). The text below
introduces the local `0.1.0.9003` development snapshot. Publish an
announcement only after reviewed source, checked release assets and
documentation links are live.

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

## Announcement draft

RdistNN 0.1.0.9003 combines conditional average treatment effects,
nonparametric regression and supplied-prediction conditional LATE with
distributional nearest-neighbor weights. Shared result methods provide
pointwise uncertainty, tuning comparisons, diagnostics and saved-fit
reuse. The latest maintenance update consolidates shared components and
strengthens tests for candidate uncertainty, numerical parity and
package installation. Estimator interfaces and inference limits are
unchanged. The guides demonstrate data-frame and matrix interfaces and
distinguish supported workflows from experimental confidence sets, PLQR
and local IV-GMM.

[Source and release
history](https://github.com/JakobJuergens/RdistNN/releases) · [Getting
started](https://JakobJuergens.github.io/RdistNN/articles/getting-started.html)

Replace the release-history link with the checked release URL when
publishing. Local documentation builds do not publish a tag, release or
website update.
