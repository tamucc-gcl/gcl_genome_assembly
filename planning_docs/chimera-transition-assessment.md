# Cohort chimera review and manual cutting

The active workflow is `workflows/chimera.nf`, implemented by `modules/misassembly.nf` and the `misassembly_*.py` scripts. The previous discovery, recovered-join, adjudication and review processes have been removed from the active workflow. Earlier investigation scripts and fixtures remain available for reproducing historical measurements; their output formats are not inputs to this workflow.

## First pass

1. Align every unordered pair of same-species assemblies, with CIGAR retained. These jobs run independently through the normal SLURM executor. Ten assemblies produce 45 pair jobs, each requesting 8 CPUs and 32 GB.
2. Build chromosome labels from the current harmonized name maps. Resolve blocks inside composites against plainly labelled cohort chromosomes using these new alignments. Both haplotypes of one individual supply one vote. Composite labels never label other composites recursively. Conflicting and unmapped bases remain unassigned. A focal individual's propagated label cannot count as independent confirmation of itself.
3. Screen every chromosome-scale assembly, including the naming reference, for transitions between chromosome-sized arms. The default minimum is 1 Mb of assigned sequence per arm (`chimera_min_arm_bp`). Small labelled islands are audited rather than promoted into chromosome-fusion decisions. Nearby measurements of the same ordered transition are grouped into one event.
4. Map the individual's HiFi reads once for each assembly with an event. Assemblies without events do not get a redundant mapping. Read continuity tolerates indels up to 50 bp, counts each primary molecule once, and is measured at 1 kb probes with 1 kb flanks. MAPQ 20 is the primary assay; MAPQ 30 and all primary placements provide sensitivity and ambiguity checks. Depth is measured separately from continuity.
5. Test verified current AGP gaps in and within 250 kb of each transition as explicit alternatives. Preserve the distinction between a detected transition and a nearby possible cut. Inspect intervening chromosome assignments before associating a gap with the transition. Assess direct gap-spanning HiFi molecules, immediately adjacent flank observability, per-library Hi-C counts, and native primary-contig placement.
6. Measure Hi-C using the pairs' actual source FASTA, matching last-scaffolding AGP and read-library identities. Compare gap assays with length-matched gaps and sequence assays with controls of the same separation. Controls require independent same-chromosome context; sequence controls also need measured HiFi continuity. Coverage matching and control counts are retained. Measured zeros are not missing measurements.
7. Publish a Markdown review and one editable decision file. Finishing may run concurrently in evidence-only mode because it uses the unchanged source assemblies. This report consistently uses pre-finishing coordinates.

A species without usable chromosome assignments or informative independent comparisons is explicitly marked unassessed. Absence of candidates is never substituted for a failed assay. Each assembly's report includes a collapsed screening coverage table.

## Suggested actions

| Action | Meaning |
| --- | --- |
| CUT | Recommend one exact verified gap for human approval. At least two independent individuals separate its flanks; there are no qualified same-chromosome or conflicting votes; both immediate HiFi flanks are observed but fewer than two molecules span the gap; no continuous native path opposes unjoining; intervening chromosome context is consistent; two libraries show calibrated contact loss. |
| RETAIN | Local HiFi supports continuity, or at least two independent individuals support same-chromosome context without contrary local votes. Broad intervals with at least 95% supported probes and no weak pocket longer than 3 kb can support retention; an independently supported exact gap-cut recommendation takes precedence. |
| SUSPECT | Chromosome separation has independent support, but the exact unsupported adjacency has not been established. The report gives the transition range, weak read-support intervals and any exact gap alternatives. |
| REVIEW | Evidence is incomplete or conflicting. The report identifies which measurements are missing or contradictory. |

These are explicit review heuristics, not claims that a biological fusion is impossible or statistically proven. Native assembled-path continuity is context, not independent molecule support. Same-chromosome context does not prove exact adjacency. Hi-C participated in assembly and is corroborating evidence. Low-confidence read placements are diagnostic only. Contact-loss calibration currently requires at least 100 within-side pairs on both flanks and five controls within a factor of ten in each flank's counts. A conservative upper ratio `(cross+3)/sqrt(left_within*right_within)` must fall below one quarter of the matched controls' tenth percentile. The exact controls and thresholds are recorded.

No recommendation selects a cut. Automatic cutting is deferred.

## Published output contract

Start at `assembly/chimeras/README.md`. It contains a compact assembly overview, suggested actions, and evidence cards separating evidence for cutting from evidence for retaining. Each candidate has a readable assembly-scoped ID such as C01; gap alternatives use C01G01, C01G02, etc. The cards show chromosomes, original scaffold coordinates, exact available gaps, direct read counts, contact counts by library and calibration status. They link to focused chromosome/continuity/depth plots and peer comparisons.

| File | Purpose |
| --- | --- |
| README.md | Cohort review, action summary and evidence cards |
| review.tsv | Editable decisions, initially all selected=NO |
| cut-instructions.md | Exact editing instructions, including adding an undetected cut |
| assembly_registry.tsv | Every assessed source assembly, FASTA path, checksum and status |
| mapping_manifest.json | Retained pre-finishing BAM/index/provenance paths for IGV or explicit reuse |
| coordinate_status.md | Coordinate contract |
| assemblies/ID/report.md | Assembly evidence, screening coverage, candidate plots and scaffold plots |
| assemblies/ID/evidence.json | Used measurements, screening audit, controls, source checksums and tool versions |
| assemblies/ID/*.tsv | Read support, chromosome blocks, Hi-C and telomere measurements, only when measured |
| assemblies/ID/plots/*.png | White-background candidate plots; whole-scaffold contact maps, telomere signal and binned contact continuity |

The binned pooled contact overview has 300 kb flanks. Candidate decisions use the separate exact 250 kb per-library assays, not this pooled overview. Temporary query FASTAs, SAMs, reference/peer PAF duplicates, unused IGV files, historical-join recovery and old diagnostic HTML are not published. Shared pair PAFs and BAMs remain in the normal Nextflow work cache. Source paths in the mapping manifest require retaining that cache.

## Manual decisions and a second run

Copy review.tsv outside generated output before editing. For a chosen gap option set `selected=YES`, `review_disposition=CUT`, reviewer and reason. Keep its supplied coordinates, source checksum and `action=UNJOIN_UNSUPPORTED`. Keep or defer an event with `selected=NO` and disposition RETAIN or DEFER. All completed reviews require a reviewer and reason.

For an internal sequence cut use `action=BREAK_PROBABLE_MISJOIN`, `localization_status=localized`, an exact reviewed cut_bp and blank gap fields. Add an undetected cut as a new row with a unique ID, copying assembly, assessment_sha256 and coordinate_stage from the registry. Coordinates are zero-based on the original pre-finishing scaffold; cut_bp divides [0,cut_bp) from [cut_bp,length). Do not substitute a broad range's midpoint for localization.

The existing source-bound cutter checks the actual FASTA checksum, scaffold identity, literal N gap, conflicting cuts and minimum piece size. It preserves every base, emits neutral fragment names and verifies exact reconstruction. Corrected sequences undergo chromosome reassignment and then finishing. The report retains original-source evidence and records submitted decisions and verified applied actions. Manual-only edits invalidate cutting, reassignment, finishing and downstream outputs; cohort discovery and read mapping remain cached.

## Full Crest run

Sync this checkout to the Crest repository first. From the project root containing data/, work/ and gcl_genome_assembly/:

```bash
mkdir -p logs
sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_pipeline.sbatch \
    data/assembly_samplesheet.csv data/hic_readsets.csv
```

This launches the real main.nf through run_assembly_checkpoint.sbatch, with the production `.nf/assembly` launch directory, work cache, `-resume` and normal per-process SLURM resources. Outputs go to `comparisons/chimera-pipeline-JOBID/results`; the controller packages the Markdown reports, tables, evidence and plots into `comparisons/chimera-pipeline-JOBID.review.tar.gz`. Existing assembly/scaffolding tasks can remain cached. The new cohort alignments and corrected evidence calculations must run on the first replacement pass.

To apply reviewed cuts in another full run:

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_pipeline.sbatch \
    data/assembly_samplesheet.csv data/hic_readsets.csv \
    "$PWD/comparisons/chimera-reviewed" \
    --chimera_break "$PWD/data/reviewed-chimera-cuts.tsv"
```

## Validation before Crest

The Python tests cover ordinary-indel continuity, rejection of large deletions/skips, separate treatment of known N gaps, MAPQ sensitivity, deleted-base depth, minor-island suppression, independent-individual voting, cohort-derived composite labels, propagation independence, PAF inversion, distance-matched controls, missing contacts, conflicting hypotheses and cut/retain recommendations. A real report-to-cutter test applies both an approved gap cut and a manually added internal cut and verifies sequence reconstruction and second-run reporting.

`tests/run_chimera_channel_reproduction.sh chimera_cohort_wiring_check.nf` exercises every replacement process with Nextflow 23.10.1, including publication and the naming-reference case without a self-PAF. This local integration uses explicitly labelled external-tool doubles; it verifies channels and script interfaces rather than real mapping accuracy. Actual minimap2/samtools/tidk execution, full-data resource use and the CTlk biological conclusions must be checked with the full Crest run above. The historical junction labels are not hard-coded into discovery or recommendations.

The local integration also passes with the `evidence` and `manual` scenarios, including per-library contact counts and delivery of corrected assemblies to finishing. A separate local resume check did not recover cached tasks, including the unchanged cutter; cache reuse has not been empirically verified in this local test environment. The Crest launch retains the production cache paths and `-resume`; check its trace for upstream cached tasks on the full run.

## Discovery performance update (2026-10-09)

Discovery now uses four CPU workers per assembly, parallelized across peer comparisons (64 GB requested). Catalog construction retains its existing two-CPU allocation. Binary-search indexes replace repeated scans of all chromosome-label intervals. Adjacent catalog labels merge only when chromosome identity and supporting individuals are equivalent; the projection index can merge equivalent assignments after independence filtering. Unassigned gaps and conflicting labels remain unassigned. Consecutive CIGAR match/mismatch operations merge, while insertions and deletions remain explicit boundaries. Published projected tracks merge equivalent adjacent intervals without bridging gaps. Out-of-scope PAF records are discarded before reverse CIGAR conversion.

Each peer and scaffold logs its start, alignment-record count, compact output-interval count and elapsed time to the task `.command.err`. Serial and parallel outputs are deterministic. Six added tests compare indexed projection with the prior algorithm, including randomized conflicting labels, both strands, indels, provenance changes and process workers. The full Python suite passes 272 tests. The fragmented-label synthetic benchmark (`python tests/benchmark_misassembly_projection.py`) produced identical assignments and approximately 460-fold faster projection on the local test; this is not a full-genome runtime estimate.

Use the same Crest launch command above after syncing the updated checkout. The alignment process and command are unchanged, so the 45 completed pairwise alignments remain eligible for cache reuse. Updated catalog/discovery helpers intentionally invalidate those stages and their dependent evidence/report tasks. Actual speed and memory use must be checked on Crest.
