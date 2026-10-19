# Agent Instructions for RdistNN

Read this file before modifying code. It is the canonical instruction file for
LLM agents and automated coding assistants working in this repository.

## Development and production branches

Develop on the local `development` branch. Keep it unpublished until the
maintainer authorizes publication; branches in a public repository are public
when pushed. Production is `master`. Submit production changes through a pull
request and use a squash merge after review and required checks. Do not push
development commits directly to `master` unless explicitly instructed to do so.
See `CONTRIBUTING.md` for the workflow.

## Mandatory read order

Before making changes, read these files in order:

1. `markdown_notes/design/01_interface_standardization_decisions.md`
2. `markdown_notes/design/02_code_design_principles.md`
3. `markdown_notes/design/03_input_validation_contract.md`
4. `markdown_notes/design/04_api_shapes_and_input_formats.md`
5. `markdown_notes/design/05_error_messages_and_warnings.md`
6. `markdown_notes/design/13_estimator_and_result_architecture.md`

For package-engineering, lifecycle, compiled-code, documentation, or CI work,
also read as relevant:

- `markdown_notes/design/06_namespace_and_dependency_policy.md`
- `markdown_notes/design/07_lifecycle_and_deprecation_policy.md`
- `markdown_notes/design/08_comments_and_function_docs.md`
- `markdown_notes/design/09_testing_and_regression_policy.md`
- `markdown_notes/design/10_rcpp_and_compiled_code_policy.md`
- `markdown_notes/design/11_examples_vignettes_and_pkgdown_policy.md`
- `markdown_notes/design/12_release_and_ci_policy.md`

For tuning-grid work read `markdown_notes/design/14_integrated_tuning_grid_reuse.md`;
for PLQR inference work read `markdown_notes/design/15_plqr_variance_estimation_rationale.md`.
For paper-replication work read `markdown_notes/design/16_affine_replication_contract.md`.
For CLATE/ratio work read `markdown_notes/design/17_clate_ratio_inference.md`.

## Maintainer documentation

- `markdown_notes/design/README.md` indexes the maintained interface, numerical
  and engineering contracts. Keep public documentation focused on current
  decisions, scientific limits and extension requirements.
- `markdown_notes/rdistnn_roadmap.md` records current capabilities and extension
  boundaries. Check it before starting new package work.
- Keep session instructions, implementation task lists, private source-library
  procedures and execution receipts outside tracked public documentation.
  Local `archives/` is ignored by Git and excluded from source packages; verify
  that any private working record is untracked before committing.

## Repository-wide non-negotiable rules

- User-facing high-level estimators may accept either:
  - explicit matrices and vectors; or
  - `data.frame` inputs, including tibbles.
- After validation and interface resolution, internal computation must proceed on
  canonical numeric matrices and vectors.
- Do not add `stopifnot()` to exported functions or to helpers on user-facing
  paths.
- Validation must happen early and produce actionable, argument-specific
  messages.
- Diagnostics should use clean user-facing `stop()` / `warning()` messages,
  generally with `call. = FALSE`.
- Unknown `return` values must be handled consistently across wrappers.
- Do not hand-edit generated files such as `NAMESPACE` or `RcppExports.*`.
- Preserve numerical parity between compiled and R engines and, for shared
  methods, the Stata companion. Follow design notes 09 and 10 for exact
  discrete decisions, justified floating-point tolerances and tests that
  establish which engine actually ran. New compiled paths need such tests.
- New user-visible validation or behavior changes must include tests.
- Do not introduce new interface conventions without updating the design notes.
- Mathematica / Wolfram Language is available through `wolframscript`. Use it
  when helpful for symbolic derivations, finite combinatorial checks, or
  estimator-formula sanity checks, and record any substantive conclusions back
  into code, tests, or design notes.

## Local CLI Tooling

Resolve tools in the environment actually running the command. In WSL, use
`command -v TOOL`; in Windows PowerShell, use `Get-Command TOOL`. A tool
installed on one machine is not evidence that it is on another machine's PATH.

- Use `rg` for source search; `jq` or Python for structured metadata.
- Source installation needs a C++17 toolchain and Rcpp; use the platform
  requirements in `README.md`. Regenerate `README.md` from `README.Rmd` with
  `devtools::build_readme()` when its source changes.
- Use `Rscript` for R work, `actionlint` for workflows, and `wolframscript`
  when a mathematical check calls for it. Follow the task-specific checks below.
- Poppler (`pdfinfo`, `pdftotext`, `pdftoppm`) and `qpdf` support document
  inspection. Windows Ghostscript is commonly `gswin64c`; use `magick` for
  ImageMagick. Check availability before promising a build or check.
- For external sources, use an available local source collection or public
  references. A personal library is not a package dependency.

## Working rules for agents

When editing existing code:

1. Prefer standardizing toward the design notes over local convenience.
2. Preserve the public API unless there is a documented reason to change it.
3. Keep the matrix/data.frame dual-interface policy at the outer boundary only.
4. Reuse or extend shared validators instead of duplicating validation logic.
5. Keep error messages specific about what was expected, what was received, and
   what the user should do instead when helpful.
6. If a code path changes user-visible behavior, update tests and docs in the
   same pass.

## In case of conflict

If code, comments, and design notes conflict, treat the numbered design notes as
current source of truth unless there is a stronger repository-level instruction
added after this file.

## External sources

Use the package bibliography, public reference documentation and source papers
for scientific claims. Check the relevant source before changing a mathematical
contract, and preserve citation keys and the distinction between proved,
conditional and illustrative claims. Personal library paths and machine-specific
access instructions belong in private configuration. Package use, tests and
maintainer documentation must remain understandable without that configuration.
