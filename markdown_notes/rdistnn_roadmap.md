# RdistNN Live Roadmap

**Updated:** 2026-10-17
**Purpose:** the single live list of open package work. Completed plans are
archived under `markdown_notes/archive/`; the design policy corpus is
`markdown_notes/design/`. Verify any item against current code before
treating it as still open.

## Current State

The initial public release focuses on CATE, nonparametric regression, and
supporting DNN weights, tuning, diagnostics, and plotting. All PLQR functions
and methods are labeled experimental, as is local IV-GMM. The October release
polish is recorded in `release/2026-10-14-public-release-polish.md`;
`RELEASING.md` describes packaging and publication, and `WEBSITE.md` contains
copy for the author's website. These scope changes do not alter the API or
computational behavior.

The pre-public overhaul is implemented (archived roadmap in
`archive/completed_plans_2026-07-11/`): paper-aligned PLQR names
(`g`, `theta_eval`, `theta_se`, `moment_se`, `theta_vcov`, `theta_omega`),
tuning constant `lambda`/`lambda_grid`, jackknife-meat coefficient sandwich
(`theta_meat = c("analytic", "jackknife")`), supplied-theta S3 semantics,
moment-first PLQR plotting, hardened direct sandwich, density/Jacobian
diagnostics, oracle demo density, API tier table, repositioned experimental
IV-GMM, and untracked `GH_logs/`/`archives/`. `devtools::check()` was clean
and the full CI matrix was dispatched on the release candidate.

On 2026-07-10/11 the CATE simulation layer additionally gained theory-facing
diagnostics for the paper's evidence contract (realized localized targets,
physical-radius diagnostics, localization-drift/centered-remainder
decomposition support consumed by
`Cond_Z_Estimation/markdown_notes/simulation/main_power_sipp_theory_alignment_2026-07-11.md`).

The public notation and package-conventions vignette now consolidates the
package-to-paper symbol map, CATE and PLQR nuisance names, effective-size
conventions, fitted-object fields, and the distinction between moment
diagnostics and coefficient inference. No pre-public documentation item
remains open in this roadmap.

## Open: Long-Term Framework (not release blockers)

1. **General conditional Z-estimation front end.** One polished entry point
   that takes a score/moment description, localizes, solves, and reports
   inference with diagnostics. The score API provides the modular pieces;
   the unified front end does not exist yet.
2. **Nuisance-fitting convenience.** Wrappers still expect externally fitted
   nuisances (`g`, `mD`, propensities). Keeping this external is a deliberate
   current design choice; revisit only after the front-end question.

The current affine simulation and SIPP replication workflows are now bundled
under `inst/paper-replication/affine/`, with unchanged paper source, production
configurations, and compressed accepted evidence. The paper repository remains
the scientific source of truth. See design note 16 for the refresh contract.

## Conditional LATE implementation

The experimental supplied-prediction CLATE feature now has a signal constructor,
paired exact-DNN ratio engine, dual interfaces, result methods and diagnostics.
The bounded generator, internal ridge-series reference adapter and separate
ratio confidence-set helper are implemented. The [CLATE collection](clate_implementation/README.md)
owns the design contracts and requirements-to-evidence record. Package and
hosted CI acceptance are recorded there as they complete.

Automatic practical learners, adaptive/grid selection and plots remain future
contracts. CLATE production simulations and SIPP participation analysis are
separate study decisions. Existing local IV-GMM keeps its distinct score and
covariance route; the frozen paper replication bundle remains independent.

## Rule

Keep this note short. When an item lands, move its detail into a dated
archive packet (or the design corpus if it becomes policy) and delete it
here.
