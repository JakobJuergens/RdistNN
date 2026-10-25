# Conditional Z-estimation paper replication

This bundle accompanies [*Orthogonal Inference for Conditional Z-Estimation
without Cross-Fitting: A Distributional Nearest-Neighbor Approach*](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z).
Its `manifest.json` identifies the accepted paper snapshot, configurations,
source hashes and saved evidence. The October 25, 2026 synchronization vendors
paper snapshot `122ad9c` and the accepted SIPP packet
`rdistnn-0.1.4-20261024` (fingerprint
`5f7203f70e52646d134089a09564a91a2aea45391df8bf1ef2dc8b93bfc836a2`).
The manifest separately records the manuscript-content snapshot.

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
its exports as validation. Exact represented distance groups average all
assigned rank weights, including zeros; every deletion recomputes that average
over its surviving groups. Nearby unequal distances stay separate.
Fits, predictions, diagnostics, paired covariance and exports are saved under `runs/paper-current/`, with a log for each stage.

These forests, outcomes and selected profiles lie outside the paper's
bounded-series primitive theorem. Frozen deletions do not refit learners;
order sensitivity does not estimate spatial bias or select an optimal order.
Full-fit/OOB agreement does not verify a stability rate. Inference under
persistent distance ties is not established. Saved fits and nuisance
arrays contain observation-level data, inherit their sharing restrictions and
are excluded from the public evidence bundle.

## Reproducibility and maintenance

Saved-evidence PDF rendering requires Cairo support in R (`capabilities("cairo")`).
The three SIPP manuscript figures are redrawn from portable summaries with embedded fonts; accepted historical evidence remains unchanged.

Re-rendering saved evidence and recomputing an estimator are different tasks.
The accepted simulations retain their original executable and configuration
identities. The October 24 raw-dollar SIPP evidence used released RdistNN
0.1.4 (`4acade538f7ac8a88fe7f837b1673ebc17daa17a`); older 0.1.0 executable
records belong to historical SIPP evidence. The accepted packet's preparatory
`input row order` label is stale: its acceptance and derivation receipts verify
actual `exact_distance_block_averaging_v1`. The separate dated clarification
preserves the accepted objects, fingerprint and receipt chain.
Newly vendored source requires complete retained geometry for current fits. A current installation records its own identity; matching version
labels or a clean source checkout do not establish that it generated the
accepted evidence. A historical-source rebuild is a rebuild unless its
executable identity also matches. Review new results before manuscript promotion.

The [replication contract](https://github.com/JakobJuergens/RdistNN/blob/master/markdown_notes/design/16_affine_replication_contract.md)
records the executable pins and bundle-refresh procedure. The paper repository
owns scientific changes; package maintainers synchronize accepted snapshots
without editing frozen scientific files.
