# RdistNN website copy

Prepared for a software entry or research-page addition on
[jakobjuergens.com](https://jakobjuergens.com/). The text below describes the
planned release version `0.1.2`. This revision is unpublished. The package
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

## 0.1.2 release announcement draft

RdistNN 0.1.2 averages weights across exactly identical
reporting-covariate rows and redistributes their block mass after every jackknife
deletion. Estimates and covariance are invariant to aligned row permutations up
to numerical precision. Shared R/C++ helpers apply the rule to fitting,
prediction, candidates, diagnostics and weight plots. Individual weights,
support and effective sample size can change; old tied saved fits require
refitting. The package warns that this practical modification is not covered
by the accompanying paper's theory. The paper supplies the underlying inference
conditions; numerical invariance does not establish those conditions.

The guides and function help describe the updated callback shape, saved-state
and independent GMM-start contracts.
The [weight guide](https://JakobJuergens.github.io/RdistNN/articles/dnn-weights.html#numerical-limits)
also describes the remaining asymptotic-tail rounding limitation.

Task-oriented guides, runnable vignettes and focused function help support
conditional-effect and regression workflows. Seeded demos and simulations
preserve caller RNG state; native calculations use two threads by default with
a configurable opt-in. This release candidate has not been published.

[Source updates](https://github.com/JakobJuergens/RdistNN/blob/master/NEWS.md) ·
[Tagged releases](https://github.com/JakobJuergens/RdistNN/releases) ·
[Getting started](https://JakobJuergens.github.io/RdistNN/articles/getting-started.html)

For a future tagged-release announcement, add the verified release URL.
This repository copy is prepared text for the personal website; package-site
deployment and personal-site publication are separate steps.
