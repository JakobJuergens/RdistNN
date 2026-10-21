# Reproduce the affine paper results

The installed `paper-replication/affine` bundle provides three
workflows: recreate manuscript displays from saved evidence, rerun the
simulation campaign, and fit the SIPP application. It is a
study-specific companion to this
[paper](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z).
For general estimation, start with [Getting
started](https://JakobJuergens.github.io/RdistNN/articles/getting-started.md).

## Locate the installed bundle

Install current production source using the
[README](https://github.com/JakobJuergens/RdistNN#installation). The
historical `v0.1.0` release predates this bundle.

``` r

bundle <- system.file("paper-replication/affine", package = "RdistNN")
list.files(bundle)
#>  [1] "configs"           "evidence.zip"      "helpers.R"        
#>  [4] "manifest.json"     "README.md"         "render_evidence.R"
#>  [7] "run_simulations.R" "run_sipp.R"        "source"           
#> [10] "verify.R"
```

Read its `README.md` for the full protocol and dependencies.
`manifest.json` pins the paper/source identity, accepted configurations
and evidence hashes. No paper checkout or original run cache is needed
to render the saved results.

## Recreate tables and figures without refitting

Install optional rendering dependencies once:

``` r

install.packages(c("digest", "jsonlite", "ggplot2", "patchwork"))
```

Run the portable launcher with a new output directory:

``` r

system2(file.path(R.home("bin"), "Rscript"), c(
  shQuote(file.path(bundle, "render_evidence.R")),
  "--output-dir=paper-displays"
))
```

It verifies the saved evidence, recreates **21 table fragments and six
PDF figures**, and compares tables with the accepted references. It
performs no simulation or estimation. Output manifests record provenance
and comparison results. Historical transformed SIPP estimates remain
saved evidence because the supplement reports that distinct target; the
current fitting workflow uses raw dollars.

Rendering uses Cairo PDF. CRAN R on macOS may require
[XQuartz](https://www.xquartz.org/), as explained in the [R for macOS
notes](https://mac.r-project.org/bin/macosx/).

## Validate the simulation workflow

A smoke run uses two small preflight draw families. It checks execution
and outputs; it does not reproduce the production study or provide
manuscript evidence.

``` r

install.packages(c("digest", "jsonlite", "future", "future.apply"))
```

``` r

system2(file.path(R.home("bin"), "Rscript"), c(
  shQuote(file.path(bundle, "run_simulations.R")),
  "--mode=smoke", "--output-dir=paper-smoke"
))
```

The full campaign is substantial and requires Linux, Python 3 and the
resource budget in the installed README. Use `--mode=production` with a
new directory; `--profiles=...` selects profiles and `--resume` requires
an identical runtime freeze. The launcher records the executed package,
seed/configuration identities and process receipts. Incomplete campaigns
are not production evidence; fresh results require review before
manuscript promotion.

## Validate or rerun the SIPP application

The application studies raw-dollar net financial assets in 9,915
complete households, with full-fit, own-arm OOB and held-out five-fold
nuisance predictions. The installed README defines the accepted
specification, profiles, orders, seed and export layout.

``` r

install.packages(c("hdm", "grf", "digest", "jsonlite"))
```

``` r

system2(file.path(R.home("bin"), "Rscript"), c(
  shQuote(file.path(bundle, "run_sipp.R")),
  "--mode=smoke", "--output-dir=sipp-smoke"
))
```

Smoke uses a small sample and validation-labelled outputs. Use
`--mode=production` with a new directory for the full application.
Output directories are not silently overwritten.

These practical forests, unbounded outcomes and data-selected
profiles/scales are outside the bounded-series primitive theorem. Frozen
deletions do not refit nuisance learners. Order sensitivity does not
estimate localization bias or select an optimal order; full-fit/OOB
agreement does not verify a stability rate.

For package-wide argument and covariance conventions, see [Notation and
conventions](https://JakobJuergens.github.io/RdistNN/articles/notation-and-conventions.md).
