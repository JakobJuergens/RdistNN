# Comments and function documentation

Comments explain purpose, invariants, non-obvious validation and design choices.
Use short conceptual blocks and descriptive names rather than narrating trivial
assignments. Internal helpers need a clear purpose; nontrivial behavior warrants
a brief contract comment even without full roxygen.

Map mathematical code to the relevant score, weighting rule or deletion path
with a short explanation and source reference. Link to the paper for derivations;
keep formulas needed to understand the implementation beside that code. Explain
non-obvious guards such as the delete-one requirement `s < n`.

Exported roxygen should cover purpose, supported input formats, constraints,
method choices and returned quantities. Prioritize user tasks and precise units
rather than exhaustive implementation detail. Examples are small, seeded,
representative and current; follow
[the examples policy](11_examples_vignettes_and_pkgdown_policy.md).

Update comments/help with behavior. Move unresolved design questions to an
issue or design note rather than leaving speculative "maybe optimize later"
comments in production code. Use direct, concrete prose that makes the next
read easier.
