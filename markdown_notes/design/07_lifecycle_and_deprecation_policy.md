# Lifecycle and deprecation policy

## Current scope

RdistNN retains its package-wide experimental lifecycle. Supported workflows
include CATE, DNN/TDNN nonparametric regression, supplied-prediction binary-IV
CLATE (including `iv_augmented_signals()`), supplied-signal affine ratios and
shared weights, tuning, reporting, diagnostics, plots and frozen-model reuse.
Computational support does not extend the corresponding pointwise theorem.

`ratio_confidence_set()` is experimental and does not establish weak-IV-uniform
coverage. All PLQR constructors, helpers, methods and shared-generic dispatch,
and local IV-GMM, remain experimental. Mark lifecycle status in help, articles
and pkgdown without changing computation or adding runtime warnings solely for
that label. [Package scope](../rdistnn_roadmap.md) records extension boundaries.

## Public changes

Version 0.1.1 records the shared interface baseline. Version 0.1.2 entered
production `master` on 2026-10-22. Version 0.1.3 replaces its identical-covariate
averaging with exact-distance averaging. The README installation section and
release tooling own the current version-pinned release; earlier tagged releases
remain immutable.
Refit saved objects from earlier tie policies; no optional legacy policy is
provided. See [RELEASING.md](../../RELEASING.md) for tag and asset publication
gates. The first development commit after production `0.1.x` establishes
`0.1.x.9000`; subsequent development commits increment the fourth component,
including small changes. Uncommitted drafts do not advance this sequence. For
example, development after production `0.1.2` starts at `0.1.2.9000`, then
`0.1.2.9001`. These are development revisions,
not public releases. Each production merge advances `0.1.x` to `0.1.(x+1)` unless
the maintainer chooses otherwise; the next development cycle appends `.9000`
to the selected production version. Read the current version from `DESCRIPTION`
and the current `NEWS.md` heading; keep them synchronized. Current fits must still
survive serialization and support the documented prediction, candidate and reporting paths.

Names, arguments and their meanings/defaults, result fields/classes, methods
and consequential warning/error behavior are public contracts. Evolve them
deliberately with documentation and regression tests; experimental status does
not excuse undocumented breakage.

## Compatibility during early development

The current policy is deliberate: support only the current documented
interfaces and complete current fitted-state schemas. Remove compatibility
layers for older package versions, including retired aliases and wrappers,
saved-state upgraders and migration shims, unless the maintainer explicitly
requests an exception. Missing or inconsistent retained state must produce a
clear refit error rather than reconstructing or silently repairing an older
fit. Keep one current implementation and preserve current-fit serialization
and documented reuse.

An earlier production release or version tag does not create a compatibility
commitment at this stage. Preserve released changelogs, immutable tags/assets
and frozen scientific evidence. Revisit a stable-release compatibility and
deprecation policy when the packages stabilize.

Keep core result fields/classes and their methods aligned. Document intentional
changes and update roxygen, affected examples/articles and regression tests in
the same change as public behavior. Any explicitly requested compatibility
exception must record its scope, purpose and intended removal.
