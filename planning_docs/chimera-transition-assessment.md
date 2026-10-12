# Cohort chimera review and manual cutting

The active workflow is `workflows/chimera.nf`, implemented by `modules/misassembly.nf` and the `misassembly_*.py` scripts. The previous discovery, recovered-join, adjudication and review processes have been removed from the active workflow. Earlier investigation scripts and fixtures remain available for reproducing historical measurements; their output formats are not inputs to this workflow.

## First pass

1. Align every unordered pair of same-species assemblies, with CIGAR retained. These jobs run independently through the normal SLURM executor. Ten assemblies produce 45 pair jobs, each requesting 8 CPUs and 32 GB.
2. Build chromosome labels from the current harmonized name maps. Resolve blocks inside composites against plainly labelled cohort chromosomes using these new alignments. Both haplotypes of one individual supply one vote. Composite labels never label other composites recursively. Conflicting and unmapped bases remain unassigned. A focal individual's propagated label cannot count as independent confirmation of itself.
3. Screen every chromosome-scale assembly, including the naming reference, for transitions between chromosome-sized arms. The default minimum is 1 Mb of assigned sequence per arm (`chimera_min_arm_bp`). Small labelled islands are audited rather than promoted into chromosome-fusion decisions. Chromosome-arm discovery uses regional dominance to tolerate small conflicting label islands. Independent-individual consensus refines the transition using 100 kb regions followed by 1 kb regions. Conflicting haplotypes cannot supply two votes. Peer-only transitions contradicted by continuous local cohort identity remain in the audit rather than the review list. Nearby measurements of the same ordered transition are grouped into one event; their union is a screening range, not the final localized interval.
4. Map the individual's HiFi reads once for each assembly with an event. Assemblies without events do not get a redundant mapping. Read continuity tolerates indels up to 50 bp, counts each primary molecule once, and is measured at 1 kb probes with 1 kb flanks. MAPQ 20 is the primary assay; MAPQ 30 and all primary placements provide sensitivity and ambiguity checks. Depth is measured separately from continuity.
5. Scan literal N runs in the assessment FASTA, recovering both current AGP gaps and gaps inherited inside scaffold components. Test verified gaps in and within 250 kb of each transition as explicit alternatives. Preserve the distinction between a detected transition and a nearby possible cut. Inspect intervening chromosome assignments before associating a gap with the transition. Assess direct gap-spanning HiFi molecules, immediately adjacent flank observability, per-library Hi-C counts, and native primary-contig placement.
6. Measure Hi-C using the pairs' actual source FASTA, matching last-scaffolding AGP and read-library identities. Compare gap assays with length-matched gaps and sequence assays with controls of the same separation. Controls require independent same-chromosome context; sequence controls also need measured HiFi continuity. Coverage matching and control counts are retained. Measured zeros are not missing measurements.
7. Publish a Markdown review and one editable decision file. Finishing may run concurrently in evidence-only mode because it uses the unchanged source assemblies. This report consistently uses pre-finishing coordinates.

A species without usable chromosome assignments or informative independent comparisons is explicitly marked unassessed. Absence of candidates is never substituted for a failed assay. Each assembly's report includes a collapsed screening coverage table.

## Suggested actions

| Action | Meaning |
| --- | --- |
| CUT | Recommend one exact verified gap for human approval. At least two independent individuals separate its flanks; there are no qualified same-chromosome or conflicting votes; both adjacent HiFi neighborhoods are observed but fewer than two molecules span the gap; no continuous native path opposes unjoining; intervening chromosome context is consistent; at least one available Hi-C library shows contact loss against suitable controls, with no calibrated library disagreeing. |
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

## Assessment rendering update (2026-10-10)

Candidate chromosome plots now batch all intervals of each chromosome/peer into one drawing collection. Every assigned interval is retained; blank gaps remain blank, and evidence JSON/tables are unchanged. This replaces hundreds of thousands of individual Matplotlib collections in CTlk. On the local synthetic test, 30,000 intervals dropped from 52 seconds to about one second; 228,704 intervals across nine peer rows rendered in 1.5 seconds. These timings measure rendering only and are not a total assessment runtime prediction. Two added tests verify exact interval geometry, preserved gaps and collection counts independent of interval counts; 274 Python tests pass.

Assessment now logs elapsed stage times to `.command.err`, including external commands, per-event read assays, Hi-C scanning, candidate/scaffold rendering and evidence writing. The previous CTlk jobs cannot identify elapsed times for individual internal stages. The broad discovered ranges (7–43 Mb) need biological/localization review separately; drawing efficiency does not validate them.

That rendering-only revision changed the assessment helper and left discovery and mapping cache inputs unchanged. The newer boundary-refinement revision below changes discovery inputs too; follow its explicit retained-BAM instructions.

## Candidate figures and refined-boundary rerun (2026-10-12)

The main candidate figure has five aligned panels: peer chromosome assignments,
HiFi spanning molecules, HiFi depth, separate per-library Hi-C contact-support
profiles, and telomere motif counts. Pale chromosome-colored backgrounds on the
four measurement panels use local independent-individual agreement. Unassigned
and conflicting regions stay neutral. Orange marks the assessed transition;
dashed red lines mark available gap cuts. Chromosome colors remain consistent
between figures.

The separate `plots/scaffold_*.contacts.png` figures contain only the pooled
contact matrix. The Hi-C profiles use 100 kb bins and 300 kb flanks. They provide
spatial context; cut recommendations use the separately recorded exact-gap
assays and matched controls. Zero contact support with observable flanks is
plotted as zero; unobservable support is a gap in the curve. Plot profiles and
regional shading records are retained in `evidence.json`; exact measured gap
contacts remain in `contact_support.tsv`.

Local replay of the archived 1517368 chromosome tracks narrowed the CTlk hap1
chr9→chr7 bracket to 63,044,000–63,051,000 bp, without supplying a known cut
coordinate. The literal gap must still be recovered from the actual assessment
FASTA and assayed on Crest. No new biological cut recommendation is established
by this local replay.

After syncing the checkout, launch from the same production directory. This
runs the full architecture and uses its existing cache. Pass the previous
published mapping manifest explicitly because discovery and its helper inputs
have changed; `-resume` alone cannot guarantee that the mapping process is
reused. The mapper checks source checksum, sample, stage and BAM dictionary
before accepting a retained BAM.

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_pipeline.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv \
  comparisons/chimera-refined-review \
  --chimera_hifi_bam_manifest "$PWD/comparisons/chimera-pipeline-1517368/results/assembly/chimeras/mapping_manifest.json"
```

The launcher uses `-resume` and `--publish_dir_mode symlink`, publishing cached
files and entire output directories through symbolic links. The review archive
dereferences links and remains self-contained. Read the new
`comparisons/chimera-refined-review/results/assembly/chimeras/README.md` and
per-assembly reports. Existing assemblies and pair alignments are cache-eligible;
discovery, local assessment, plotting and review rerun. Retained BAM validation
runs without remapping HiFi reads. Every cut remains unselected until manual
approval through `review.tsv`.


## Evidence repair validation (2026-10-12)

Policy `cohort-review-v2` fixes the issues found in the refined-review packet.
A single-individual fusion screen is rejected when at least two independent
individuals have qualified same-chromosome flanks and outnumber the separate
and conflicting local votes. Missing consensus bins cannot prevent this
classification. A real cohort transition or independently conflicting evidence
remains reviewable. Rejected screens retain original observations and local
votes in the audit, with a collapsed human-readable table.

Gap association has three states: confirmed, unresolved, contradicted. A
cohort-localized gap within the recorded spatial resolution can be associated
without inventing sequence homology through an N gap. Actual contradictory
intervening orders are retained separately. Association does not authorize a
cut. Every gap lists all unmet cut prerequisites.

Immediate read probes remain recorded. Adjacent flank observability now uses
10 kb neighborhoods, sampled at 1 kb with 1 kb anchors, requiring at least five
probes, median support of two molecules and at least 80% supported probes.
This is an explicit review heuristic; it does not convert missing read support
into proof of a misassembly. Boundary-only weak observations are exported as
sampled positions, rather than zero-width weak intervals.

Control pools include up to 96 spatially separated sites, selected across
scaffolds before filling remaining slots. All selected controls are counted in
the same per-library Hi-C scan, including controls rejected for chromosome
context or read continuity. `control_registry`, `control_links`, and per-target
`control_audit` retain geometry, context, counts and rejection reasons.
`control_support.tsv` exposes these measurements for review. The calibration
requirements remain 100 within-side pairs per flank, five controls and a
factor-of-ten coverage match on each flank. No threshold was relaxed to force
an outcome.

Two calibrated Hi-C libraries opposing measured read continuity produce an
explicit REVIEW conflict. A single-library disagreement is recorded even when
read continuity supports retention. Retained boundaries are collapsed and
excluded from the actionable candidate count. No absent evidence file gets a
report link. Figures name their assembly and display ambiguous read/depth
measurements on separate, labelled right axes.

To validate the repaired evidence on Crest after syncing code:

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_pipeline.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv \
  comparisons/chimera-evidence-repaired \
  --chimera_hifi_bam_manifest "$PWD/comparisons/chimera-refined-review/results/assembly/chimeras/mapping_manifest.json"
```

The existing launcher resumes the production cache, links published outputs,
and keeps cutting disabled. Neighborhood read assays and expanded control
measurements must run on Crest; replaying an archived packet cannot create
those measurements. Review its new evidence before choosing any cuts or
continuing the wider pipeline improvement pass.

## Publication correction (2026-10-12)

Crest run 1532235 still launched with a later `--publish_dir_mode copy`,
overriding the checkpoint launcher's earlier `link` setting. Its 962 queued
publication operations therefore did not test the local link-mode repair.
The comparison launcher now explicitly requests `symlink`, which publishes
both files and directories without copying their contents. Nextflow 23.10.1
verification: `tests/publication_mode_check.nf` publishes an individual file
and a directory containing 1,000 files, asserting both are symbolic links and
the nested contents are readable. The review archive still dereferences links.
Verify the final publication argument in the Crest launch log after syncing;
an explicit trailing `--publish_dir_mode symlink` also overrides older launchers.

## Single-library decisions and visible transitions (2026-10-12)

Policy `cohort-review-v3` works with one Hi-C library. A cut recommendation requires
at least one available library with contact loss against suitable controls, plus
the existing independent chromosome, exact gap, adjacent read-observability and
native-path requirements. A calibrated library without demonstrated loss prevents
an unqualified cut recommendation; an uncalibrated extra library is unavailable
evidence rather than a veto. Manual selection remains mandatory for applying cuts.

Any calibrated Hi-C loss conflicting with HiFi continuity triggers REVIEW, including
loss in just one library. MAPQ-20 continuity that is not confirmed at MAPQ 30 also
triggers REVIEW. Conflicting calibrated library results trigger REVIEW regardless
of the read-continuity verdict. No thresholds for contact/control qualification
were relaxed.

All detected chromosome transitions now appear in the main action table and both
main/per-assembly evidence sections. Counts distinguish detection from provisional
retention. Same-chromosome false screens remain in the collapsed rejection audit.
Reports explain why contact loss cannot be judged when controls are unavailable,
and explicitly list weak intervals from the stricter read check.

Validation: 297 tests passed; the Nextflow 23.10.1 cohort integration passed with
one Hi-C library using external-tool doubles. Replaying the actual repaired-symlink
measurements gives hap1 C01 SUSPECT, C02/C03 RETAIN, C04 REVIEW, and hap2 C01/C02
REVIEW. All regenerated report file links resolve. This replay changes policy and
presentation, not measured evidence, and applies no cuts.

After syncing, validate the updated production workflow with:

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_pipeline.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv \
  comparisons/chimera-review-policy-v3 \
  --chimera_hifi_bam_manifest "$PWD/comparisons/chimera-evidence-repaired-symlink/results/assembly/chimeras/mapping_manifest.json" \
  --publish_dir_mode symlink
```

## Final report and manual-interface corrections

The final reporting pass preserves the v3 biological classifications and uses
the existing evidence. Sister haplotypes are compared both through independently
detected chromosome-pair transitions and through local flank measurements. A shared
chromosome pair is explicitly not treated as proof of identical junctions.

The review file now retains `event_id`, `event_assessment`, `event_reason` and the
parent recommendation separately from `option_advice` and `option_association`.
Unassociated nearby gaps are labelled as not suggested cuts in the report, overview
and review file. Exact gap options are visible alongside the transition interval.
Neighboring read support is distinguished from absent support across an actual gap.
The false library-disagreement message is removed when all measured libraries show
no demonstrated contact loss.

Validation: 301 tests pass, including generated manual gap cuts, manually added
internal cuts, and sequence reconstruction. The actual v3 packet regenerates with
six unselected decision rows, unchanged classifications and no broken links or
candidate anchors. This completes the identified report/interface fixes. Cut
selection and verification are the next phase; automation remains deferred.

Production validation now runs the same CHIMERA, ASSEMBLY_FINISHING,
FINALIZE_ASSEMBLY, COLLECT_NAME_MAPS, ASSEMBLY_RUN_SUMMARY and FINAL_READ_EVIDENCE
modules used by main.nf. Small synthetic fixtures and explicit external-tool
doubles verify software interfaces; they do not validate mapping accuracy or
biological recommendations. Finishing deliberately changes sequence lengths so
stale final name maps fail the test. Manual selections must survive reassignment,
finishing, final naming and the final mapping reference digest.

The production run summary links the generated chimera review when detection is
enabled. Disabled detection omits that section. The normal pipeline launcher
packages this summary with the review archive and checks that it was published.

After syncing, run from the production project root:

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_pipeline.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv \
  comparisons/chimera-production-final \
  --chimera_hifi_bam_manifest "$PWD/comparisons/chimera-review-policy-v3/results/assembly/chimeras/mapping_manifest.json"
```

This launches main.nf with production SLURM resources, resume enabled, the existing
production work directory, retained source-validated chimera HiFi BAMs, a separate
output directory and symlink publication. It generates evidence and review files;
cuts remain unselected. Review results/assembly/chimeras/README.md and review.tsv.
Keep edited decisions outside generated outputs. A subsequent production run may
supply --chimera_break /absolute/path/to/reviewed.tsv to apply explicitly selected
cuts and rerun affected reassignment, finishing, naming and final mappings.
Automation remains deferred.

Validation of the production integration (2026-10-12): all 301 Python tests pass.
Nextflow 23.10.1 integration passes for one-library evidence and two-library manual
cut application through finishing, final naming, name-map collection, linked run
summary and final HiFi reference digests. A second manual invocation reuses cached
alignment, discovery, mapping, assessment, cutting, reassignment, finishing,
finalization and final read mapping. The run summary regenerates because it records
the invocation name. The direct-JAR local harness now passes official Java module
options into its child JVM so its asynchronous cache writer persists entries.
Full biological validation still requires the production Crest run above.

Regional evidence improvements (2026-10-12)
-----------------------------------------

Policy cohort-regional-review-v4 retains the original arm-transition bracket if
fine localization lacks two independent, qualified chromosome-flank observations.
Clipped projected blocks retain their original alignment density. The chr7/chr12
screen is supported by approximately 13–15 Mb of chr7 sequence in each independent
haplotype comparison; the old narrow placement relied on a small assignment island.

Assessment searches the entire transition bracket plus 500 kb on each side for
MAPQ20 losses and sustained MAPQ30 losses. Overlapping regions are merged, up to
24 are directly assayed, and omitted regions and assay-edge censoring are reported.
Each has independent chromosome context, all three HiFi placement tiers, actual
bridging-read identities with observed competing placements, and separately
calibrated per-library Hi-C counts. Alternatives have ranked review rows and blank
cut coordinates; their assay probes are never selected cuts. Nearby sustained
losses prevent local continuity from being interpreted as validation of the join.

Focal chromosomes use contrasting blue/orange consistently across haplotypes;
other chromosomes use muted colors and a separate legend. The report and plots
distinguish the original transition bracket, ranked alternative regions and exact
gap options. Manual cuts remain explicitly reviewed, source-bound selections.

Validation: 311 general Python tests pass, including regional losses, placement
ambiguity, missing measurements, fine-assignment islands, read-level bridge audits,
alignment-density clipping and manual decision safety. Nextflow integration tests
cover one-library evidence and two-library manual cutting through finishing,
renaming and final read mapping, including cache resume. Tool doubles test wiring;
biological accuracy is assessed with the full Crest rerun below.

The saved CTlk profiles nominate the chr4/14 low-confidence region around
35.217–35.496 Mb, the 42.831 Mb weak point, and the hap2 support loss around
39.864–39.894 Mb. All-primary placements restore support in the first and last
regions. Fresh alternative-site contact measurements are required before any
new interpretation; old nominal-boundary measurements are not transferred.

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_pipeline.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv \
  comparisons/chimera-regional-review \
  --chimera_hifi_bam_manifest "$PWD/comparisons/chimera-production-final/results/assembly/chimeras/mapping_manifest.json"
```

The manifest reuses source-validated BAMs from the last real run; it is optional
in a fresh production run. Pairwise alignment and unchanged upstream assembly
processes resume from the normal cache. Discovery, assessment and reporting
regenerate. Outputs use symbolic links, and the launcher writes a self-contained
comparisons/chimera-regional-review.review.tar.gz for reassessment.
