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

Version 0.1.1 records the shared interface baseline. Version 0.1.2 is the next
planned shared, version-pinned release; it remains unpublished while prepared
on `development`. See [RELEASING.md](../../RELEASING.md) for the publication
gates. A fourth component starting at 9000 marks development revisions rather
than public releases. Earlier development versions
have no external users, so their retired interfaces and fit schemas are removed
without deprecation shims. Current fits must still survive serialization and
support the documented prediction, candidate and reporting paths.

Names, arguments and their meanings/defaults, result fields/classes, methods
and consequential warning/error behavior are public contracts. Evolve them
deliberately with documentation and regression tests; experimental status does
not excuse undocumented breakage.

After the first shared release, prefer adding the improved interface, preserving compatibility where feasible,
warning on deprecated use, updating examples, then removing the old interface
once the transition is visible and tested. Warnings name the deprecated use,
replacement and removal timing/condition. Use a hard error when the old behavior
is harmful/misleading, compatibility creates serious inconsistency and the
migration is clear.

Keep core result fields/classes and their methods aligned. Do not silently
rename fields. A short-lived shim must record why it exists, its replacement
and intended removal. Update roxygen, affected examples/articles and tests in
the same change as public behavior.
