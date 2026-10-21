# Conditional Z-estimation paper replication

This bundle accompanies [*Orthogonal Inference for Conditional Z-Estimation
without Cross-Fitting: A Distributional Nearest-Neighbor Approach*](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z).
Its `manifest.json` identifies the accepted paper snapshot, configurations,
source hashes and saved evidence, refreshed October 14, 2026.

Choose the task before launching a run:

| Entry point | Purpose |
| --- | --- |
| `render_evidence.R` | Recreate accepted tables and figures from saved evidence, without fitting models. |
| `run_simulations.R` | Recompute the study with the actual installed RdistNN package. |
| `run_sipp.R` | Fit and export the raw-dollar SIPP application with the installed package. |

Locate the installed bundle in R:

```r
bundle <- system.file("paper-replication/affine", package = "RdistNN")
list.files(bundle)
```

In the shell commands below, replace `/path/to/bundle` with that path. Quote
paths containing spaces. Every launcher accepts `--help`, works from any
working directory and requires a **new external output directory**; it never
writes into the installed package.

## Recreate the accepted paper displays

Install `digest`, `jsonlite`, `ggplot2` and `patchwork`, then run:

```sh
Rscript /path/to/bundle/verify.R
Rscript /path/to/bundle/render_evidence.R --output-dir=paper-displays
```

Rendering requires R's Cairo PDF device. With CRAN R on macOS, install
[XQuartz](https://www.xquartz.org/) for the required X11 libraries; see the
[R for macOS notes](https://mac.r-project.org/bin/macosx/).

The renderer verifies the evidence archive, rebuilds eight simulation and
thirteen SIPP table fragments, compares them with their accepted references,
and produces six PDF figures. Output includes hashes, display manifests and
`replication.json`. Historical transformed SIPP estimates remain saved evidence
for the supplement's distinct target; their deprecated fitting workflow is
not included. No learner fitting or native run caches are needed for rendering.

## Run the paper simulations

Install `RdistNN`, `digest`, `jsonlite`, `future` and `future.apply`:

```sh
Rscript /path/to/bundle/run_simulations.R --mode=smoke --output-dir=paper-smoke
Rscript /path/to/bundle/run_simulations.R --mode=production --output-dir=paper-production
```

Smoke checks both affine applications, nuisance-training routes and DNN/kernel
localization using two draw families at n=128. It validates execution rather
than recreating production evidence or its large-dictionary endpoint.

Production runs four accepted profiles (`dnn`, `kernel_comparison`, `M3`,
`settings_comparison`) with 10,000 draw families and 28 cells. The M3 endpoint
has n=531,441. The full launch requires Linux and Python 3 for process
monitoring; it uses four workers, one inner thread, eight-family chunks,
per-profile time limits and a sampled 40-GiB RSS cap. This is a substantial
campaign. Smoke and evidence rendering also run on other R platforms.

Select a subset with `--profiles=dnn,kernel_comparison`, for example. Use
`--resume` only with the same output directory and matching source, dependency,
seed and configuration identities. Chunks, failures, resource measurements,
installed-package hashes and completion receipts are retained. A partial run
is not accepted production. The study's order `floor(n^(3/4))` (multiplier 0.4
in C2 sensitivity cells) overrides the general package default.

## Run the current SIPP application

Install `hdm`, `grf`, `RdistNN`, `digest` and `jsonlite`. The public
`hdm::pension` data require no private input:

```sh
Rscript /path/to/bundle/run_sipp.R --mode=smoke --output-dir=sipp-smoke
Rscript /path/to/bundle/run_sipp.R --mode=production --output-dir=sipp-production
```

Production uses 9,915 complete households and raw net financial assets in
1991 dollars. It compares full-fit, own-arm OOB and held-out predictions at
seven age-40 income profiles, using common final DNN weights and an eight-order
sensitivity grid. Smoke uses 300 households, 50 trees and two folds and labels
its exports as validation. Fits, predictions, diagnostics, paired covariance
and exports are saved under `runs/paper-current/`, with a log for each stage.

These forests, outcomes and selected profiles lie outside the paper's
bounded-series primitive theorem. Frozen deletions do not refit learners;
order sensitivity does not estimate spatial bias or select an optimal order.
Full-fit/OOB agreement does not verify a stability rate. Saved fits and nuisance
arrays contain observation-level data, inherit their sharing restrictions and
are excluded from the public evidence bundle.

## Reproducibility and maintenance

Re-rendering saved evidence and recomputing an estimator are different tasks.
The accepted simulations and SIPP run have distinct historical executable
identities. A current installation records its own identity; matching version
labels or a clean source checkout do not establish that it generated the
accepted evidence. A historical-source rebuild is a rebuild unless its
executable identity also matches. Review new results before manuscript promotion.

The [replication contract](https://github.com/JakobJuergens/RdistNN/blob/master/markdown_notes/design/16_affine_replication_contract.md)
records the executable pins and bundle-refresh procedure. The paper repository
owns scientific changes; package maintainers synchronize accepted snapshots
without editing frozen scientific files.
