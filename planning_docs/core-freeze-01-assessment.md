# Core-freeze-01 — Crest acceptance

Reviewed 2026-11-05 from `core-freeze-01.review.tar.gz`.

## Acceptance

The final native core-freeze confirmation passes. The interface/documentation and
version-record cleanup gate is complete; PSMC is the next planned implementation.
This accepts the tested core routes, not unfinished optional analyses or biological
completeness of the graph.

- Wrapper exit: 0; Nextflow SUCCESS, with 64 completed and 875 cached tasks.
- No failed, ignored, retried or aborted tasks; no WARN/ERROR entries in the
  supplied Nextflow log. Approximate launch-to-finish time: 70 minutes.
- Core integration audit: PASS, ten finalized assemblies, five within-individual
  comparisons, one cached Cactus construction and zero rebuilt constructions.
- Published report validation: PASS, 543 local links checked on Crest.
- Software version issue table: header only, zero issues. All ten final Cooler
  observations are labelled and retain version 0.10.4. All ten TIDK_EXPLORE
  observations are labelled single-line records retaining tidk 0.2.65.
- The native tidk database warning is retained within its version observation;
  it no longer creates stray malformed records. This is not a new failed task.
- Four manual cuts remain exactly reconstructed/verified: three hap1 cuts and
  one hap2 cut. Both chr7/chr12 transitions remain retained by review.
- Sequence-graph, clipping, regional and variant audits pass. Source accounting
  remains REVIEW for the same five chromosome-scale representation flags; these
  are retained limitations, not newly detected assembly misjoins.

## Evidence comparison and execution scope

Compared 84 selected files directly with `detailed-report-01.review.tar.gz`:
final name maps, painting tables, event evidence packets, QC metrics and selected
pangenome JSON/TSV evidence. All 84 are byte-identical; none differ or are absent
from the baseline. This is a selected evidence comparison, not hashing complete
FASTA/BAM/graph products. The on-Crest integration audit separately verifies
baseline identity/name-map/reference consistency and cache reuse.

Assembly, scaffolding, read mappings, manual cut products and Cactus construction
stayed cached. Reruns included 41 tidk tasks, ten final contact-map tasks, eight
pangenome summary/plot tasks and five planning/report/provenance tasks. The
pangenome summaries therefore reran more broadly than the minimal version-record
forecast, but their selected retained evidence is unchanged. The supplied trace
shows the reruns, not a definitive cache-miss cause. No biological changes or
required additional biological run follow from this observation.

## Figure review and remaining presentation notes

Visually inspected the Crest-produced main chromosome structure figure,
haplotype/individual growth figure and CTlk hap2 final 1-Mb contact map. The
structure figure retains ordered assembly panels, a separate unplaced-total scale
and clear naming-versus-painting provenance. Growth shows both biological units.
Its quorum curves need not be monotone when an integer presence threshold changes;
they must not be mistaken for the union accumulation curve used for growth fitting.

The whole-scaffold contact map still has crowded tick labels, particularly small
and unplaced pieces. Record selective ticks/label spacing for its future scheduled
figure revision; the map renders and this is not a core-freeze blocker. This
review does not claim a new R migration of the legacy contact-map renderer.

The report records a working-tree revision rather than a release commit. Preserve
this archive and the synced source snapshot; a future release tag should identify
the exact released source. Do not invent a commit identity from the run name.

## Next work

Proceed to both agreed PSMC input routes and their callable-mask, fit, bootstrap,
scaling and report integration. The next pangenome visit separately includes the
[agreed graph-QC/MultiQC, growth-fit/classification and confidence-band exploration](pangenome-qc-growth-followup.md).
