# Examples, vignettes and pkgdown

Documentation is part of the public interface and should teach the preferred
current API. Give each document a task:

| Surface | Role |
| --- | --- |
| README | Installation, a first interpretable result and links onward |
| Vignette | One coherent estimation, diagnostics, extension or replication workflow |
| Function help | Complete syntax, shapes, returned quantities and relevant constraints |
| Paper | Theory, derivations and statistical assumptions |
| Maintainer reference | Implementation and verification contracts |

Keep practical interpretation limits beside outputs. Separate paper-specific
replication from general onboarding. Where both input routes are supported,
explain matrix/vector and data-frame usage without repeating equivalent fits
throughout the collection.

## Executable examples

Use small, fast, seeded examples with minimal dependencies. Do not rely on
ambient RNG state. Guard Suggested packages at the point of use, explain
optional workflows and keep unrelated portions usable when a dependency is
absent. Richer workflows belong in focused vignettes rather than long roxygen
examples. Execute ordinary lightweight chunks during rendering; explicitly
identify costly installation/production steps that readers run separately.

## Site and maintenance

Pkgdown articles and reference groups should reflect current core primitives,
user estimators, extension APIs and replication tools. Keep every public topic
discoverable. Update affected examples, article text/code and site grouping
alongside user-visible behavior. Verify deterministic execution, optional
package guards, links and generated output using the commands in
[CONTRIBUTING.md](../../CONTRIBUTING.md).
