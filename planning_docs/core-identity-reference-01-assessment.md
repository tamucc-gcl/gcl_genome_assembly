# Native identity/reference acceptance — 2026-10-28

Packet: `core-identity-reference-01.review.tar.gz`.

## Outcome

Accepted. Exit 0; four completed tasks, zero failed/ignored/retried tasks. The local
taxonomy snapshot was reused through its normal stored-output contract. No assembly
or graph tasks ran. The acceptance-report task's upstream assertions confirm both
synthetic sample identities survived and name/taxid routes agree.

| Check | Observed result |
|---|---|
| Exact scientific-name lookup | Spratelloides delicatulus → 373251 |
| Installed taxonkit | Complete ten-column taxonomy record, with species/genus and full lineage |
| Installed reference retrieval | NC_009588.1, Spratelloides delicatulus, 16,617 bp |
| Output integrity guards | Reference task succeeded after checking a single nonempty matching pair, minimum length and organism annotation |
| Combined evidence | identity_reference_acceptance.md plus taxonomy/reference TSVs present |

Task durations including scheduler overhead: name resolution 24.6 seconds,
reference retrieval 4.9 seconds, taxonomy 9.9 seconds, evidence 4.9 seconds. Peak
RSS was 15.4 MB for streaming name resolution, 40.9 MB for retrieval and 525.6 MB
for taxonkit. The reference query used the reduced one-CPU/4-GB request.

The first-launch `-resume` warning is expected in this new isolated launch
directory. Unmatched full-pipeline config selectors are also expected because this
diagnostic intentionally imports only the relevant processes. Neither is an
execution failure. There are no ERROR entries in the controller log.

The raw superkingdom field reads `unclassified cellular organisms superkingdom`.
The full lineage explicitly contains Eukaryota, Metazoa and the vertebrate/fish
clades. The production kingdom, BUSCO, genetic-code and telomere helpers use
kingdom/full-lineage matches, so this raw rank label does not invalidate those
classifications. Preserve the raw record; no special-case species patch is needed.

The reference FASTA/GenBank contents are intentionally not in this compact packet.
Their existence/integrity was checked by the native task; this review does not
independently reparse their sequences or certify the best biological reference.

## Checklist reconciliation

Category 2 (supported general-purpose input handling) is accepted from the parser/
mixed-species/skip matrix, native HiFi and shotgun runs, the empty-reference final
inventory repair, 17 identity/reference scenarios and this native name/ID check.

Category 3 (organelle recovery and required failure/fallback contracts) is accepted
from the controlled positive/negative workflow matrix plus native HiFi and
fragmented short-read recovery/annotation and this reference test. Biological
completeness is not required for a correctly recorded partial/no-recovery outcome.

These are component acceptance decisions, not a claim that every possible taxon,
input biology or external network failure has been tested. The broader retry,
resource, dependency pinning and optional-comparison checks remain category 7.
Unified current-run reporting remains category 6. The final combined cache-based
core integration check remains outstanding before freezing the core interface.

No repeat of this identity/reference pilot, shotgun pilot or graph build is
indicated by the supplied results. Next work is the remaining operational audit
and bounded optional-comparison validation, followed by Markdown reporting and
the combined integration check. PSMC remains the first post-core extension.
