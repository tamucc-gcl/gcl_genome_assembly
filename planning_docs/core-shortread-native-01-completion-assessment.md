# Native short-read completion accepted — 2026-10-27

Reviewed the replacement `core-shortread-native-01.assembly.review.tar.gz`,
run 1553560, after the final-inventory wiring corrections.

## Acceptance

Exit 0. Nextflow reports three completed tasks, 15 cached, zero failures, ignored
failures or retries. No WARN/ERROR entries appear in the controller log.
The new tasks are eligibility (23.7 seconds), shared reference (5 seconds) and
run summary (5 seconds). All successful native assembly and organelle tasks stayed
cached. The earlier launcher exit-1/missing-summary defect is resolved on Crest.

The native shotgun route's execution/completion gate is accepted. This does not
declare biological completeness or close every contract in categories 2/3/7.
No further repeat of this same pilot is needed without a concrete new defect.

## Output checks

- The summary records one expected and one finalized `spades` primary assembly,
  with zero Hi-C read sets. No artificial haplotype pair was created.
- Eligibility correctly labels it `collapsed_or_primary`, finalized,
  non-chromosome-scale, and not a graph member. Its species cohort is withheld
  for insufficient eligible individuals and no suitable chromosome-scale reference.
  Graph construction was explicitly disabled.
- A final linear reference is nevertheless selected from the actual finalized
  assembly, independently of graph eligibility. Its reference dictionary matches
  every published FAI entry: 1,973 sequences and 6,358,369 bp. The manifest records
  final-stage SHA256
  `2de5cf00bf4f2fcae5b21265b35a6fb0b892aa5943438a871aaef399d10a6412`.
  This review verifies recorded identity/dictionary consistency, not an independent
  rehash of the full FASTA, which is intentionally omitted from the packet.
- All 26 common result files in the preceding local inspection are byte-identical.
  There is no change to assembly names or the recovered organelle output.
- The previously omitted FCS action and taxonomy reports are now included. The
  action report has zero action rows: no removal, trimming or splitting was
  requested. Taxonomic similarities alone are not decontamination actions.
- Six annotated linear mitochondrial fragments still total 4,945 bp. This remains
  partial recovery. Native annotation and the empty-target retry correction have
  completed; no complete mitogenome is claimed.

## Remaining limitations and follow-up

The 6.36 Mb nuclear test assembly remains far below the genome-size estimate.
This intentionally reduced-work pilot validates execution and interfaces, not
production short-read assembly performance or biological completeness. Do not use
its selected linear reference as evidence that it is biologically adequate for a
later analysis; consumers still need their own suitability checks. Purging
optimization remains deferred.

The summary contains the standalone sentinel text `No reference alignment.` when
chimera detection is disabled. Remove this presentation artifact during the agreed
Markdown reporting pass; no workflow rerun is justified solely for this wording.

Remaining core gates are the recorded species-name/taxid and missing-reference
failure checks, the deterministic/transient retry and dependency audit, and
run-scoped report/operational acceptance. The native short-read pilot is no longer
an outstanding gate. Source plans remain the canonical 13-category checklist and
the input/failure/cleanup checkpoint.
