# RdistNN Live Roadmap

**Updated:** 2026-07-22
**Purpose:** the single live list of open package work. Completed plans are
archived under `markdown_notes/archive/`; the design policy corpus is
`markdown_notes/design/`. Verify any item against current code before
treating it as still open.

## Current State

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
3. **Empirical-application packaging.** The SIPP 401(k) workflow lives in the
   paper repository as manuscript replication code. Package-level packaging
   is needed only if the application becomes part of the package story.

## Rule

Keep this note short. When an item lands, move its detail into a dated
archive packet (or the design corpus if it becomes policy) and delete it
here.
