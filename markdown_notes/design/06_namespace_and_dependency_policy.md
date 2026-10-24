# Namespace and dependency policy

Export primary user entry points, documented extension APIs, S3 methods that
need visibility and lower-level primitives with a clear external use. A public
function needs a user story, roxygen documentation and public-interface tests,
and must remain meaningful outside its current module. Keep validation,
wrapper-specific plumbing, transitional helpers and likely-to-change internals
private. Avoid unnecessary export renaming/removal.

Constructors and registries intentionally expose extension interfaces; this
does not make their implementation helpers public. Index every public topic in
pkgdown and give primary workflows a clear article/README route.

## Class and registry identities

Owned S3 classes use an `rdistnn_` prefix; summary classes use
`rdistnn_summary_` followed by the estimator class's meaningful stem. Keep
subclass relationships explicit and preserve distinct estimands. The exported
constructor/generic names and semantic argument or extension-kind tokens are
separate contracts, not automatically renamed with a class. See
[Advanced R, S3 classes](https://adv-r.hadley.nz/s3.html).

Use existing specification constructors and validation owners. Revalidating or
registering a current specification preserves its class identity without adding
duplicate class entries. Estimator
adapters each own final result assembly; a new single-use wrapper that merely
sets a class does not consolidate invariants. Shared retained-state validation
checks the completed scientific payload before reporting or reuse. Do not
introduce a common superclass just to reuse display code.

Every accepted registry key must be discoverable through its public name
lister, including dot-prefixed and nonsyntactic keys. Registration, lookup and
discovery share the same name contract; do not infer a hidden-name convention
from `ls()` defaults. Keep listing and available-name diagnostics at the common
registry owner.

## Generated files

Edit roxygen in source and regenerate `NAMESPACE` and Rd topics through the
package tooling. Do not hand-edit generated namespace content. Treat drift as
an ownership/toolchain problem; use the declared generator version in
[CONTRIBUTING.md](../../CONTRIBUTING.md).

## Dependencies

Prefer base R and a small, stable dependency footprint. Add Imports only when
normal package functionality requires them and the benefit justifies install,
CI, portability and maintenance costs. Examples, tests, vignettes, benchmarks
and optional comparisons belong in Suggests.

Check Suggested packages at their use sites, give actionable unavailable-feature
messages and guard documentation examples. An absent optional package must not
disable unrelated functionality. Audit both the dependency's use and the code
that requests it before removing it.
