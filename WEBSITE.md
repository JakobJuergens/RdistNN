# RdistNN website copy

Prepared for a software entry or research-page addition on
[jakobjuergens.com](https://jakobjuergens.com/). The text below describes the
`0.1.0.9005` documentation and maintenance update for production `master`. The package
documentation is live at its canonical pkgdown URL. The latest tagged release
remains `v0.1.0`; a tagged-release announcement should link to its verified
source archive and release assets.

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

## Source-update announcement draft

RdistNN 0.1.0.9005 makes conditional-effect and regression workflows easier to
learn, with shorter guides, runnable vignettes and focused function help.
It adds training-sample counts through `nobs()`, preserves the caller's saved
RNG state for seeded demos and simulations, and caps native calculations at
two threads by default with a configurable opt-in. Existing pointwise
uncertainty, diagnostics, tuning comparisons and saved-fit reuse retain their
documented inference scope. The accompanying paper supplies the theory.

[Source updates](https://github.com/JakobJuergens/RdistNN/blob/master/NEWS.md) ·
[Tagged releases](https://github.com/JakobJuergens/RdistNN/releases) ·
[Getting started](https://JakobJuergens.github.io/RdistNN/articles/getting-started.html)

For a future tagged-release announcement, add the verified release URL.
This repository copy is prepared text for the personal website; package-site
deployment and personal-site publication are separate steps.
