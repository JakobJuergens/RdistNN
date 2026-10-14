# RdistNN website copy

Prepared for a software entry or research-page addition on
[jakobjuergens.com](https://jakobjuergens.com/). The text below
introduces the 0.1.0 release; publish it after the release and
documentation links are live.

## Software entry

**RdistNN — local treatment effects and nonparametric regression in R**

RdistNN implements distributional nearest-neighbor estimation for
conditional average treatment effects and nonparametric regression. It
combines local estimates with delete-1 jackknife standard errors, tuning
comparisons, and diagnostics for weights and treatment overlap. CATE
nuisance predictions are supplied by the user; valid inference requires
the conditions developed in the accompanying research. The initial
release focuses on CATE and regression. PLQR and local IV-GMM
functionality are experimental.

[Source and installation](https://github.com/JakobJuergens/RdistNN) ·
[Documentation](https://JakobJuergens.github.io/RdistNN/) · [Getting
started](https://JakobJuergens.github.io/RdistNN/articles/getting-started.html)
· [Accompanying
paper](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z)

## Release announcement

RdistNN 0.1.0 is available: an R package for conditional average
treatment effects and nonparametric regression with distributional
nearest-neighbor weights. The package includes jackknife standard
errors, localization tuning, and weight and overlap diagnostics. The
getting-started guide demonstrates both data-frame and matrix interfaces
and explains the inference conditions. PLQR and local IV-GMM remain
experimental.

[Release and source
package](https://github.com/JakobJuergens/RdistNN/releases/tag/v0.1.0) ·
[Getting
started](https://JakobJuergens.github.io/RdistNN/articles/getting-started.html)
