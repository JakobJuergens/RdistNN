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
