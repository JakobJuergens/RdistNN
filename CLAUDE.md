# Claude Instructions

Read `AGENTS.md` first.

`AGENTS.md` is the canonical repository instruction file for coding agents in
this repo. The package design standards, interface rules, validation contract,
and cleanup plan are defined there and in the linked files under
`markdown_notes/design/`.

Before changing code, especially in exported functions or shared helpers, follow
`AGENTS.md` and its required read order.

For mathematical derivations, symbolic checks, or finite combinatorial
verification, Mathematica / Wolfram Language is available through
`wolframscript` and should be used when helpful; record nontrivial conclusions
in code, tests, or design notes rather than relying on an ephemeral session.
