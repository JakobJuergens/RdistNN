# Current affine paper replication

This bundle represents the affine variant of *Orthogonal Inference for
Conditional Z-Estimation without Cross-Fitting: A Distributional
Nearest-Neighbor Approach*, refreshed October 14, 2026.
`manifest.json` pins the exact paper commit, manuscript and supplement source
snapshot, accepted production run, SIPP run,
source hashes and saved evidence. The independent coupled-score manuscript
is outside this bundle. Superseded study launchers are not included.

Locate these files after installation:

```r
bundle <- system.file("paper-replication/affine", package = "RdistNN")
list.files(bundle)
```

Run the entry points with `Rscript`, passing an explicit **new output
directory**. They work from any working directory and never write into the
installed package. In the commands below, replace `/path/to/affine` with
`bundle`. Paths containing spaces should be quoted.

## Recreate the accepted paper displays

Install the optional dependencies `digest`, `jsonlite`, `ggplot2` and
`patchwork`. Then run:

```bash
Rscript /path/to/affine/verify.R
Rscript /path/to/affine/render_evidence.R --output-dir=affine-paper-displays
```

The renderer verifies and decompresses the accepted evidence, generates eight
simulation table fragments and thirteen SIPP table fragments, and produces six
PDF figures (three per study). It runs both upstream verification routines and
compares every generated table with its accepted reference. It does not draw
simulation data, read native run caches, or fit learners. The output includes
the input hashes, display manifests and `replication.json`.

Historical transformed SIPP estimates are saved evidence because the current
supplement reports them as a distinct target. Their deprecated fitting
workflow is not included. The newly overhauled application uses raw dollars.

## Run the current affine simulations

Required packages: `RdistNN`, `digest`, `jsonlite`, `future`, `future.apply`.

```bash
Rscript /path/to/affine/run_simulations.R --mode=smoke --output-dir=affine-smoke
Rscript /path/to/affine/run_simulations.R --mode=production --output-dir=affine-production
```

Smoke uses two preflight draw families at n=128. It checks CATE, residualized
slope, full fitting, pooled five-fold fitting, oracle nuisances, DNN and
population-ESS-matched Gaussian localization. It is **validation only** and
does not recreate the large-dictionary endpoint or production evidence.

Production uses the four unchanged accepted configurations:

| Profile | Cells | Main replications per cell | Main + control rows |
|---|---:|---:|---:|
| DNN C0/C1/C2 | 14 | 10,000 | 896,000 |
| DNN/kernel comparison | 9 | 10,000 | 1,080,000 |
| M=3 endpoint, n=531,441 | 1 | 10,000 | 60,400 |
| C0H/C2S settings | 4 | 10,000 | 480,000 |
| Total | 28 | | 2,516,400 |

There are 10,000 common independent draw families; 480,000 DNN rows are copies
across profiles. Both applications have point target 1/8 at X=(1/2,1/2).
Adjustment uses four additional Z coordinates and clipped growing ridge
series. The default order is floor(n^(3/4)), with multiplier 0.4 in the C2
sensitivity cells. This overrides the general package tuning default for this
study. The variance routes are frozen-score jackknife, direct moment and
direct squared-weight estimates, including joint slope coefficient covariance.

The full launch uses the paper's Linux/Python3 process monitor, four workers,
one inner thread, eight-family chunks, per-profile CPU/wall limits and a sampled
40-GiB RSS cap. Budget roughly the paper's measured campaign resources; the M3
profile is substantial. Smoke and evidence rendering also work on other R
platforms. Production writes immutable chunks, actual installed package/DLL
hashes, seed maps, failures, process measurements and summaries. After all four
profiles finish, it runs the upstream production audit and saved-data analysis.

Select profiles with `--profiles=dnn,kernel_comparison,M3,settings_comparison`.
Resume the same directory explicitly with `--resume`; the source, dependency,
seed and configuration freeze must match. Audit and analysis stages are
published only after their completion receipts and hashes pass; interrupted
directories are preserved and rerun. R-version differences in configuration
serialization are reconciled only after exact equality of the settings.
A stopped or partial campaign is not
accepted production. The proposal's `not_started` field is a historical
prelaunch snapshot; completed run manifests determine execution status.

The historical dependency pin is RdistNN commit
`40060f155cc4d84d59ee5af119f6a2534c569bab`. Fresh runs use and record the actual
installed package. Matching version labels alone do not establish executable
identity. The original run's R/package metadata are retained in the evidence
packet. Fresh outputs require review before manuscript promotion.

## Run the current SIPP application

Required packages: `hdm`, `grf`, `RdistNN`, `digest`, `jsonlite`.
The data are the public `hdm::pension` sample; no private data are required.

```bash
Rscript /path/to/affine/run_sipp.R --mode=smoke --output-dir=sipp-smoke
Rscript /path/to/affine/run_sipp.R --mode=production --output-dir=sipp-production
```

Production uses all 9,915 complete households, raw accumulated net financial
assets in 1991 dollars, eligibility as treatment, seven age-40 income profiles,
baseline adjustment, sample-scaled age/asinh-income reporting coordinates,
2,000-tree forests, seed 20260709 and five folds. It compares literal all-tree
full-fit predictions, own-arm OOB predictions from the same full forests, and
held-out predictions. Final DNN weights are common across mechanisms. Default
s=397; the fitted order grid is 298,397,496,745. Exploratory localization uses
99,149,198,248,298,397,496,745 and retains all 168 profile/mechanism/order rows.

The launcher fits, saves models and prediction arrays, computes frozen-score
diagnostics and paired component covariance, exports the affine application,
and exports localization sensitivity. Each stage has a log. Results live in
`runs/affine-current/`; the run cannot overwrite an earlier output. Smoke uses
300 households, 50 trees and two folds, with explicit validation export labels.

The empirical forests, unbounded dollars, selected scales/profiles and observed
geometry are outside the bounded-series primitive theorem. Frozen deletions
do not refit learners; order sensitivity does not estimate spatial bias or
select an optimal order. Full-fit/OOB agreement is not a stability-rate check.

## Maintain the bundle

The paper repository remains the scientific source of truth. From the package
checkout, run `tools/replication/sync_affine.py --paper PAPER_CHECKOUT --campaign
ACCEPTED_CAMPAIGN_ROOT` after an accepted paper update. It copies unchanged
source and configurations, preserves evidence bytes in a deterministic ZIP,
and updates the manifest and reference table hashes. All launchers are covered
by the installed-package replication tests.
The shared paper bootstrap is replaced by a small explicit-path adapter so
the unused older CATE orchestration is not installed. All nine frozen affine
scientific files retain their accepted hashes.
The rendering adapter supplies the pinned source identity and uses forward
slashes and literal path prefixes for portable output directories.
