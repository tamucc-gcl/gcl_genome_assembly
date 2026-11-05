# Detailed-report-01 — Crest acceptance

Reviewed 2026-11-05 from `detailed-report-01.review.tar.gz`.

## Execution and publication

- Launcher exit status: 0.
- Nextflow: 931 cached, eight completed, zero failed/ignored/retried tasks.
- Core integration audit: PASS; ten finalized assemblies, five within-individual
  comparison pairs, Cactus construction cached and not rebuilt.
- Published-report validation: PASS; 543 local links checked, including nested
  detailed reports. All 133 expected-output checks are PRESENT.
- Newly executed work is limited to the comparison plan, chimera review, run
  summary, sharing/regional plot renderers, version collection and the two report
  renderers. No assembly, mapping, graph construction or evidence assay reran.
- Chimera review took 47.6 seconds and about 696 MB peak RSS. Sharing/regional
  figures took 3.0/4.4 seconds. No new performance issue is indicated.

## Evidence preservation

Compared with main-report-02, all 29 checked baseline files are byte-identical:
ten final name maps, ten painting tables, six graph/regional audit JSONs, QC
metrics, assembly eligibility and the version-issues table. All six CTlk event
evidence.json packets are also byte-identical.

Each event's R render audit reconciles every assigned peer interval, consensus
interval, confidence-tier spanning-profile row, depth row, contact-profile row and
telomere row against its evidence packet. All 24 junction contact maps have PNG
and PDF outputs. Their sparse exports preserve matrix dimensions, raw totals and
zero/nonzero cell counts; original NPZ matrices are byte-identical to baseline.

Four cuts remain applied and verified, and the two chr7/chr12 events remain
retained by review. The five existing graph representation flags are unchanged;
they remain graph-review findings, not proposed new assembly cuts.

## Report and figure review

The combined evidence figure and the placement-confidence figure both remain in
the detailed reports. Regional and close-up PDFs are linked. Matching placement
site IDs and leader lines connect the top-panel markers to the lower bar groups.
Submitted cuts remain distinct green markers. Per-library junction maps carry
the resolution, original coordinate key and shared count scale; pooled whole-
scaffold images remain explicitly distinguished.

The revised combined close-up, placement-confidence figure, per-library junction
map and sharing partition were visually inspected from the Crest-produced files.
Focal chromosome colours, primary-versus-all-placement axes, readable contact
ticks, measured bar labels and PDF outputs are present. The run summary's chimera
link is restored. Twenty GenomeScope PNGs are now included in the review archive.
The pangenome hub includes the measured individual-sharing table while retaining
separate biological and technical evidence links.

No new blocking report, rendering or evidence regression was identified.

## Existing housekeeping and remaining gate

The provenance collector still lists 30 issues from cached malformed/unlabelled
version records (bare `0.10.4` and `0.2.65` lines, including associated no-labelled-
record notices). The entire issues table is byte-identical to main-report-02;
these are not new failures or evidence changes. Do not rerun expensive upstream
tasks merely to relabel old version files or invent tool labels from bare values.
Carry the producer/version-format cleanup into the interface/documentation freeze.

The native reporting integration subgate is accepted. Category 6 remains open for
the final interface/documentation freeze and disposition of the existing metadata
housekeeping. No further biological investigation or repeated report run is
required by this assessment. PSMC remains the first extension after the core gate.
