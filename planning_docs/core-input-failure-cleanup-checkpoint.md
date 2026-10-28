# Input, failure and operational checkpoint — 2026-10-25

Core graph acceptance is complete. This checkpoint continues categories 2, 3 and 7
without reopening assembly cuts, purging optimization or biological rearrangement
decisions. The source/tool rationale is now maintained in
[Pangenome tool changes](../docs/pangenome-tool-changes.md).

## Completed cleanup

- Remove the five completed pangenome experiment families: assignment replay,
  deletion-policy replay, grouping comparison, gate trace and chunk comparison.
  Remove their 20 modules/adapters/fixtures, five CLI switches, selectors, workflow
  channels and test-harness branches. Historical evidence remains in assessment
  documents, run packets and Git history.
- Keep production guarded filtering, corrected-tool build, native fragment probe,
  interval accounting, graph/node audits and their focused regression tests.
  Ownership/round-trip tests now import the permanent production helper directly.
- Do not change the accepted graph process command, its staged production helpers,
  vendor source directory, build environment or effective resource requests.
  Normal resume should reuse expensive graph work; do not delete caches or old data.

## Bounded verification matrix

The actual input parser, reference process, production graph workflow and auditors
run in local Nextflow integration tests. External assembly/graph executables are
doubles; those tests do not establish biological quality of a short-read pilot.

| Contract | Evidence/check |
|---|---|
| HiFi-only, shotgun-only and mixed cohort | Actual parser selects hifiasm/SPAdes, expected haplotype counts and capabilities |
| HiFi plus shotgun reads in one sample | HiFi assembly route remains selected |
| Shotgun plus Hi-C | Accepted shotgun route; no assertion that unimplemented short-read Hi-C scaffolding ran |
| One/two Hi-C libraries | Stable sample/library/readset identities and complete pairs |
| Hi-C-only rows | Explicit skipped assembly; not an invented successful assembly |
| Missing files/mates | Invalid row recorded in assembly-only mode; graph mode rejects an invalid cohort input |
| Malformed table/unsafe identifiers | Duplicate headers, field-count mismatch and unsafe sample names fail |
| Duplicate biological/read identities | Duplicate sample, duplicated file and repeated Hi-C readset fail |
| Conflicting Hi-C sources | Inline plus library-table reads fail rather than double include |
| Species/reference separation | No cross-species pooling; singleton reference independent of graph admission; overrides validated |
| Disabled graph/no admitted species | Actual top-level graph workflow closes with no graph reports/tasks |
| Graph/tool failures | Nonzero tool exit, fragment-header incompatibility and missing required graph roles fail explicitly |
| Organelle disabled/no complete mitogenome | Nuclear sample continues with an explicit empty bait set |
| Failed organelle executable | Failure is not recast as successful no-recovery |
| Cache | Unchanged production tool, graph and source audits resume normally |

Run `python3 tests/run_core_contract_checks.py` in the configured local test
environment to execute the 20 integration scenarios. The production tool build
requires its C++ compiler; normal Crest runs use the existing pinned conda environment.

Verification completed: all 20 scenarios passed (18 matrix cases plus the two
top-level disabled/empty checks), followed by successful production cache reuse.
The Windows focused suite passed 95 Python tests: 72 pangenome, 19 eligibility and
four shared-reference tests; one additional GNU-sed-specific test was skipped on
Windows. The unchanged upstream sed expression also executed in the Linux native
header-probe integration. Logs are in `comparisons/core-cleanup-validation/`.

## Historical remaining-work list and current acceptance

**2026-10-28 acceptance:** categories 2 and 3 are now accepted from the completed
local matrix and native HiFi/shotgun/identity-reference checks. Items 1 and 2 below
are complete within their supported execution/contract scope. Items 3 and 4 remain
category 7, and item 5 remains category 6. The final combined core integration
check is still required. See
[native identity/reference acceptance](core-identity-reference-01-assessment.md).

The subsequent [organelle/taxonomy pass](organelle-failure-contracts-2026-10-25.md)
now covers the enabled-route and missing-context contracts below with actual
workflow modules, deterministic external CLI doubles and real GenBank conversion.
Native biological pilots and name-to-taxid lookup acceptance remain distinct work.

**2026-10-27 update:** the native short-read execution pilot and final completion
resume are now accepted from run 1553560: three new/15 cached tasks, exit 0,
eligibility/reference/summary present. This closes the native-execution requirement
in item 1 below. The reduced output does not establish biological assembly quality.
See [completion assessment](core-shortread-native-01-completion-assessment.md).
The other failure, lookup, operational and reporting requirements remain open.

**2026-10-28 update:** item 2 now has implemented name-to-taxid handling, explicit
lookup failures, strict reference output checks and targeted retry rules, with
17 passing local workflow cases. Native lookup/reference acceptance is prepared
as a separate, small Crest check. See
[identity/reference contracts](identity-reference-contracts-2026-10-28.md).

1. Verify the full short-read nuclear route and the enabled short-read organelle
   branch, including successful circular/linear recovery, genuine no-recovery,
   execution failure, database setup failure and annotation failure. Current
   parser/no-recovery checks are not acceptance of those complete native routes.
2. Exercise taxonomy/reference lookup failure and missing species-specific
   references, including species resolved by name versus taxid. Confirm that a
   missing dependency cannot silently drop samples from downstream joins.
3. Finish deterministic-versus-transient retry/resource/dependency/publication
   review outside the accepted graph path. Preserve ordinary Nextflow caching;
   do not introduce cache reconstruction or dataset-specific path rescue.
4. Validate existing optional comparison plots against finalized assemblies and
   species-aware pair selection before main report integration. The branch
   now uses isolated pinned Python environments; native acceptance remains open.
5. Consolidate run-scoped Markdown reporting (category 6), including residual graph
   retention flags, excluded source intervals, cut dispositions and HiFi alternative
   placements. Disabled sections must be absent; no stale directory-wide discovery.

These are the next bounded tests, not a request for another expensive graph rebuild.
PSMC remains the first extension after core report/input/failure/cleanup acceptance.

## 2026-10-28 finalized comparison repair

The optional final comparison branch now uses metadata/species pair selection,
final FASTA/FAI dictionaries, strict synchronous PAF measurement and separate pinned
alignment/plot environments. Runtime SETUP_PAFR installation and retired R plot
wrappers were removed; accepted pangenome R figures retain their environment.
Eight unit checks, twelve workflow cases, normal resume, full-main input compilation
and the native-check wrapper with controlled external tools passed locally.
[Native comparison acceptance](final-comparison-contracts-2026-10-28.md) is next.
Category 7 remains open for that run, remaining dependency/provenance cleanup and
combined core integration. Category 6 remains the unified Markdown reporting pass.

### Native comparison execution result

`core-final-comparisons-01` completed 13/13 tasks with no failures or retries.
Final input fingerprints and five within-sample pairs reconcile. Sparse retained
alignment evidence and presentation limitations prevent complete visualization
acceptance; see [review and bounded follow-up](core-final-comparisons-01-assessment.md).
No assembly or graph rebuild is indicated.

### R synteny restoration (2026-10-28)

The optional pairwise figures again use pafr/gggenomes with ggplot2 and the
established minimap2 mapping-chain mode. The Matplotlib renderer is retired.
Raw/filtered mapping retention and interval-union coverage are now measured;
R rendering, workflow cases and cache reuse passed locally. Native recovery of
expected mapping coverage remains open. See [restoration and one-pair Crest command](final-comparisons-r-restoration-2026-10-28.md).
The user prefers R/ggplot2 for plotting throughout the pipeline. Remaining active
Python figure renderers are migration items for the category 6 reporting pass.

### Native R comparison pilot accepted

`core-final-comparisons-r-01` passed all five tasks and restored CMat chromosome
synteny: 95.88%/96.44% retained chain-span coverage. Final input dictionaries are
unchanged from the sparse run. Both R figures were inspected. No additional code
fix is indicated; the five-pair cohort check is next. See
[pilot assessment and full-cohort command](core-final-comparisons-r-01-assessment.md).
