# Agent Instructions for RdistNN

Read this file before modifying code. It is the canonical instruction
file for LLM agents and automated coding assistants working in this
repository.

## Mandatory read order

Before making changes, read these files in order:

1.  `markdown_notes/design/01_interface_standardization_decisions.md`
2.  `markdown_notes/design/02_code_design_principles.md`
3.  `markdown_notes/design/03_input_validation_contract.md`
4.  `markdown_notes/design/04_api_shapes_and_input_formats.md`
5.  `markdown_notes/design/05_error_messages_and_warnings.md`
6.  `markdown_notes/design/13_existing_codebase_standardization_plan.md`

For package-engineering, lifecycle, compiled-code, documentation, or CI
work, also read as relevant:

- `markdown_notes/design/06_namespace_and_dependency_policy.md`
- `markdown_notes/design/07_lifecycle_and_deprecation_policy.md`
- `markdown_notes/design/08_comments_and_function_docs.md`
- `markdown_notes/design/09_testing_and_regression_policy.md`
- `markdown_notes/design/10_rcpp_and_compiled_code_policy.md`
- `markdown_notes/design/11_examples_vignettes_and_pkgdown_policy.md`
- `markdown_notes/design/12_release_and_ci_policy.md`

For tuning-grid work read
`markdown_notes/design/14_integrated_tuning_grid_reuse.md`; for PLQR
inference work read
`markdown_notes/design/15_plqr_variance_estimation_rationale.md`.

## Markdown note layout

- `markdown_notes/design/` is the canonical, ordered design-policy
  corpus.
- `markdown_notes/rdistnn_roadmap.md` is the compact live roadmap; read
  it before starting new package work. The pre-public overhaul is
  implemented; the public notation page is implemented in
  `vignettes/notation-and-conventions.Rmd`. Consult the roadmap for
  remaining long-term framework and packaging work.
- Completed or superseded implementation plans are archived with status
  headers under `markdown_notes/archive/` (packets dated 2026-07-09 and
  2026-07-11, including the implemented pre-public roadmap and the
  historical package audit note); treat them as provenance, not as live
  instructions — they largely predate the 2026-07-08 API renames.

## Repository-wide non-negotiable rules

- User-facing high-level estimators may accept either:
  - explicit matrices and vectors; or
  - `data.frame` inputs, including tibbles.
- After validation and interface resolution, internal computation must
  proceed on canonical numeric matrices and vectors.
- Do not add [`stopifnot()`](https://rdrr.io/r/base/stopifnot.html) to
  exported functions or to helpers on user-facing paths.
- Validation must happen early and produce actionable, argument-specific
  messages.
- Diagnostics should use clean user-facing
  [`stop()`](https://rdrr.io/r/base/stop.html) /
  [`warning()`](https://rdrr.io/r/base/warning.html) messages, generally
  with `call. = FALSE`.
- Unknown `return` values must be handled consistently across wrappers.
- Do not hand-edit generated files such as `NAMESPACE` or
  `RcppExports.*`.
- New compiled paths must include fallback-equivalence tests.
- New user-visible validation or behavior changes must include tests.
- Do not introduce new interface conventions without updating the design
  notes.
- Mathematica / Wolfram Language is available through `wolframscript`.
  Use it when helpful for symbolic derivations, finite combinatorial
  checks, or estimator-formula sanity checks, and record any substantive
  conclusions back into code, tests, or design notes.

## Local CLI Tooling

Resolve tools in the environment actually running the command. In WSL,
use `command -v TOOL`; in Windows PowerShell, use `Get-Command TOOL`. A
tool installed on one machine is not evidence that it is on another
machine’s PATH.

- Use `rg` for source search; `jq` or Python for structured metadata.
- Source installation needs a C++17 toolchain and Rcpp; use the platform
  requirements in `README.md`. Regenerate `README.md` from `README.Rmd`
  with `devtools::build_readme()` when its source changes.
- Use `Rscript` for R work, `actionlint` for workflows, and
  `wolframscript` when a mathematical check calls for it. Follow the
  task-specific checks below.
- Poppler (`pdfinfo`, `pdftotext`, `pdftoppm`) and `qpdf` support
  document inspection. Windows Ghostscript is commonly `gswin64c`; use
  `magick` for ImageMagick. Check availability before promising a build
  or check.
- For Literature documents, use the shared catalog/extraction workflow
  below; do not replace its preserved evidence with ad hoc or in-place
  OCR.

## Working rules for agents

When editing existing code:

1.  Prefer standardizing toward the design notes over local convenience.
2.  Preserve the public API unless there is a documented reason to
    change it.
3.  Keep the matrix/data.frame dual-interface policy at the outer
    boundary only.
4.  Reuse or extend shared validators instead of duplicating validation
    logic.
5.  Keep error messages specific about what was expected, what was
    received, and what the user should do instead when helpful.
6.  If a code path changes user-visible behavior, update tests and docs
    in the same pass.

## In case of conflict

If code, comments, and design notes conflict, treat the numbered design
notes as current source of truth unless there is a stronger
repository-level instruction added after this file.

## Literature and source lookup

For local literature work, use `literature-vault`; use
`obsidian-paper-notes` for Obsidian paper/book notes. If a skill is not
exposed, read its maintained `SKILL.md` under the vault’s `_ai/skills/`.
Read the vault’s `AGENTS.md` and `_ai/AGENTS.md` before source or
citation work.

The WSL vault is `/home/jakob/research/Literature`; Windows uses
`D:\Jakob_Clouds\Dropbox\Literature` on the desktop and
`C:\Users\jakob\Dropbox\Literature` on the laptop. Preserve both Windows
layouts and resolve the current machine’s paths. Project-local
`01_Literature` collections and extraction caches are separate from this
full vault.

Use the vault’s `_bibliography/scripts/registry.py` by absolute path:
`search` for passages, `search --find` for discovery, and `lookup` for
citation identity. Preserve existing project keys and unresolved
identities. Check `status`’s `catalog_schema_version` and
`variants DOCUMENT_SHA`: baseline-only catalogs use original sidecars
plus `supplemental` MinerU reads; registered v5 catalogs support
preferred `read` and `--variant all --compare`. Follow the vault’s
`_extractions/README.md` and `_bibliography/QUALITY.md` for exact
commands and review scope. Preferred text is not verification; inspect
relied-on mathematics, tables, figures, and uncertain exact wording
unless a current scoped visual review covers it. Keep extraction
identities and physical/printed pages distinct. The personal vault is a
research resource, not a repository build dependency.
