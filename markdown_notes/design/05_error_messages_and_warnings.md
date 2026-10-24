# Error messages and warnings

This note records the reporting contract for implementation and maintenance.
Shared argument validators, localization notices, retained-state boundaries,
grid recovery and supported ratio summaries use the condition owners described
below. Owned validation uses a semantic condition family, including one-off
checks; an existing sentence does not justify bypassing that identity. Host
errors, opaque callback exceptions and native failures retain their identities.
The [testing policy](09_testing_and_regression_policy.md)
defines the required behavioral evidence.

## Shared ownership and language conventions

Extend the existing shared validators and fitted-state/extra accessors rather
than adding estimator-specific message systems. Prefer base R conditions and
a small dependency footprint; richer console formatting alone does not justify
a new runtime dependency. Public wrappers own validation and user context;
numerical kernels own computation and return numerical facts/statuses.

R and Stata should agree on a problem's meaning, consequences and permitted
recovery. They retain their own argument syntax, display conventions and
handling mechanisms. Stata uses its standard return codes and distinguishes
error output from suppressible text notices; R uses its condition protocols.

## Severity and timing

| Situation | Reporting rule |
| --- | --- |
| Invalid input or infeasible requested main order | Error before expensive computation, naming the input and constraint |
| Incomplete/inconsistent current scientific state | Error identifying the state problem and requesting refitting |
| Valid fit lacks a required optional extra | Distinct missing-extra error with the exact retention/refit action |
| Requested optional feature lacks a dependency | Feature-specific error naming the dependency and installation action |
| Automatic acceleration is unavailable | Use the documented equivalent portable path without warning by default |
| Explicitly required engine is unavailable | Error explaining the selector and available alternative |
| Exact represented distance ties | One qualified warning per public call |
| Weak localization or ignored unknown extras | One warning per issue family; explain the consequence and useful action |
| Infeasible optional candidates | Preserve candidate rows/reasons and give one bounded summary |
| Expected partial estimates or inference | Preserve estimator-specific statuses and rows; summarize requested unavailable output |
| Negative first stage with a defined ratio | Descriptive status/note; do not invent an identification failure |
| Uncertainty deliberately disabled | Description, not a warning about missing uncertainty |
| Unexpected invariant, resource failure or interrupt | Propagate/stop with appropriate context; do not silently discard or retry |

Warn when computation remains defined but the user needs to understand a
consequence or qualification. Issue ignored-input warnings before expensive
work. Shared reporting must preserve different estimator contracts: ordinary
means reject nonfinite estimates, whereas paired ratios can legitimately return
status-bearing unavailable ratios, including at every target. A common reporter
must not turn all such cases into the same generic failure.

Preserve existing mathematical criteria when consolidating messages. Signed
TDNN concentration differs from probability-weight ESS; give method-specific
advice. Changes to inference assumptions or a score kernel's supported semantics
need their own justification, not a cosmetic severity change.

## Informative, bounded copy

Write for a user who knows the call but not private implementation names.
State the problem, the useful expected/received fact, and a remedy when known.
Name arguments verbatim, usually in backticks. For example:

- "`k` must satisfy k <= s; got k = 5 and s = 3."
- "Delete-one inference requires s < n; got s = 100 and n = 100. Reduce `s`,
  or use `se = FALSE` if uncertainty is not required."
- "Column 'm1_oracle', selected by `m1_col`, was not found in `data`. Check
  the selector against the available column names."

Distinguish type, shape, missingness, nonfinite values, absent fields and
semantic conflicts. Give actual lengths/dimensions when useful. Reserve
internal-inconsistency wording for genuine package invariants, with a useful
report-a-reproducible-example action; ordinary input errors are not package bugs.
A fresh computed range failure needs rescaling guidance, not an older-fit error.

Use complete sentences and consistent plain punctuation. Show at most five
example locations/names and the omitted count. Abbreviate individual names or
values longer than 80 characters explicitly. Keep the remedy and the route to
full status details visible; aim for primary messages under 800 characters
rather than relying on the host's truncation. Arbitrary user-set display limits
remain outside the package's control.

Do not dump full arrays, outcomes, nuisance predictions, fitted objects or call
environments. Store bounded type/shape/scalar summaries in conditions. Escape
control characters in displayed names. In Stata, render dynamic text literally
rather than interpreting it as SMCL; keep owned help markup separate. Native
format strings must treat dynamic strings as data.

## Condition identities and selective recovery

The shared R constructors use small base condition objects with `message`
and `call`, package classes `rdistnn_error` or `rdistnn_warning`, and a semantic
`code`. Useful bounded fields include `argument`, `expected`, `actual`, `hint`
and a bounded parent message/class summary when wrapping a known failure. Full
parent conditions may contain large or private inputs and are not retained.
Recurring
recoverable families need a specific identity; one-off sentences do not each
need a new subclass. In particular, distinguish `rdistnn_extra_unavailable`
from invalid retained state and use `rdistnn_candidate_failure` only for
explicitly declared optional-candidate recovery. Those classes also inherit
the ordinary package/error condition classes.

The shared vocabulary separates `argument_invalid`, `state_invalid`,
`extra_unavailable`, `dependency_missing`, `engine_unavailable`,
`numerical_range`, `candidate_failure`, `callback_invalid` and
`internal_invariant`. Classify by cause: omitted optional uncertainty is an
extra, corrupted retained state is a state error, invalid callback output is a
callback error, and freshly computed range failure is numerical. Do not use a
generic bad-state helper for these different situations. Warning
families include `distance_ties`, `weak_localization`, `unknown_return`,
`infeasible_candidates`, `grid_single_value`, `partial_result`, `method_failed`
and the developer score's `score_symmetry`.
Keep the small mapping with the shared constructors; do not read a catalogue
from disk on each signal or introduce a parallel public Stata numeric code space.

Use `call. = FALSE` for ordinary string-based user errors/warnings unless the
call adds useful context. For an owned condition object, set `call = NULL` in
the object and pass it alone to `stop()` or `warning()`. A bounded `operation`
field and message prefix identify the outer public package function, including
registered S3 methods, without retaining its argument expressions or frame.
Resolve that name only on the signaling path. Shared choice validation preserves
base R's default and partial-match semantics while naming the public argument.
Preserve base handling,
including `muffleWarning`, user handlers and `options(warn = -1/0/1/2)`; do not
change global options to control package output.

Recover only from an explicitly understood condition. A dashboard may omit an
optional panel whose retained extra is absent; malformed state, invalid explicit
plot options and unexpected rendering errors propagate. An optional candidate
may retain a declared failure row; the main/default candidate still fails.
Warnings promoted to errors, interrupts, allocation failures and unexpected
callback errors must not silently become dropped panels or failed-candidate
rows. A numerical-range class alone does not authorize recovery everywhere.

This rule also applies to descriptive plot captions and tuning-grid rebuilding.
Captions consume validated retained localization facts without invoking full
diagnostics just to format concentration text. This works without optional
weights and across fitted classes with different diagnostic interfaces. Other
optional panels can omit a recognized missing extra; corrupt state or unexpected
errors propagate. Plotting preserves unrelated warnings, calling handlers and
promotion rather than blanket-suppressing the calculation. Validation aggregators catch only owned input-validation conditions;
they must not relabel resource errors as bad arguments. A solver's singularity
precheck is a numerical fact, not permission to catch every solver exception.

## Call scope and grouping reuse

Exactly equal represented localization distances at any target trigger the
shared warning that inference with persistent distance ties is not established.
The paper assumes a continuous distance distribution, under which the earlier
and current rules agree almost surely. Preserve this qualification without
suggesting jitter or a larger `s` as a remedy for exact ties.

Extend the existing distance/reuse scope at the public orchestration boundary.
Emit each recurring warning family at most once per public call, including
nested reconstruction, candidates and multi-panel plots. Independent calls
warn independently. Do not persist a session-global or fitted-object
"already warned" flag. Scope holds a small set of emitted-family flags;
complete affected counts come from the owning batch/status table, not by
adding overlapping per-panel counts.

Reuse already-computed groups and their derived tie information. Reporting
must not recompute distances, re-sort, reconstruct groups or rescan the same
group columns for every panel or candidate. Validate current retained state
before deriving call-local facts; never trust an unchecked serialized tie flag.
Deletion-refitted geometry has its own actual groups. Frozen geometry still
uses actual surviving block sizes for arithmetic. Keep reporting context
separate from scientific distance `ctx` and out of observation/deletion loops.

Fresh distance preparation derives `groups` and `tied_targets` in the same
pass. Current retained-state validation checks equality against the cached
distance ordering and derives those facts without reconstructing distance
groups. A bounded call-local cache holds the current groups/facts and validated
model/report references. The largest sample preparation and one scratch entry
keep deletion refits from evicting full-sample facts. Candidate kernels receive
facts directly, without repeated scope lookup. Unchanged nested consumers reuse them; changed state
is revalidated. This relies on ordinary R copy-on-modify semantics, not on a
persisted trust flag. Native extensions must not mutate shared retained arrays
in place. Scope references disappear when the public operation exits.

Fits and standalone weight results expose `meta$target_notices`, with original
`target_id`, `distance_ties`, `weight_concentration` and `weak_localization`
columns, and `meta$n_tied_targets` / `meta$n_weak_targets` totals. Concentration
is the existing inverse sum of squared final weights; for signed TDNN it is
not probability-weight ESS. These facts remain available without optional
weight extras or visible warnings. Candidates update them for the chosen order.
Retained-state validation checks target alignment, flags against validated
groups and stored concentration, and aggregate counts. Missing current fields
require refitting. These are diagnostic facts, never warning-history flags.

Simulation runners recover only explicitly declared method failures or
adapter-produced unavailable-result statuses. Opaque external learner, data
generation and tuning errors propagate, since an unclassed error may represent
a resource failure. Warnings promoted by the caller remain errors, and warning
handlers retain their normal chance to muffle the original warning. Expected
method failures retain method, reason and seed details, with one batch summary.
This does not certify a numerical failure as an identification diagnosis.

Construct text only when it will be signaled. Reuse existing diagnostic/status
passes, with bounded examples and output proportional to issue families rather
than observations/targets. Do not remove safety checks to make reporting cheap.
Do not flush queued warnings during fatal-error cleanup: `warn = 2` could then
replace the original error. Native host signaling follows
[the compiled-code policy](10_rcpp_and_compiled_code_policy.md).

## Sources and validation

The design uses [base R conditions](https://stat.ethz.ch/R-manual/R-patched/library/base/html/conditions.html),
[R warning semantics](https://stat.ethz.ch/R-manual/R-patched/library/base/html/warning.html),
[message-writing guidance](https://style.tidyverse.org/errors.html),
[Advanced R's custom conditions](https://adv-r.hadley.nz/conditions.html#custom-conditions),
[R Packages' condition tests](https://r-pkgs.org/testing-basics.html#testing-errors),
[public error context](https://rlang.r-lib.org/reference/topic-error-call.html)
and [Stata display conventions](https://www.stata.com/manuals/pdisplay.pdf).
Condition identity, copy, recovery and cost are separate test contracts;
see [the testing policy](09_testing_and_regression_policy.md). Validation order
remains owned by [the input contract](03_input_validation_contract.md).
