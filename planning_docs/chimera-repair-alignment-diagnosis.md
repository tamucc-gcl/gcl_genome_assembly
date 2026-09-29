# Chimera repair: alignment diagnosis and completion plan

2026-09-29. Static inspection of uploaded PAF/AGP/report text and pipeline source only. No aligner, assembly, analysis script or tests executed.

## Evidence and diagnosis

Inputs: chimera-comparison-20260929-110257-reports.tar.gz and chimera-comparison-20260929-110257-alignments.tar.gz. Current single-library run: high_euclid, job 1500578; success, 103 completed and 186 cached. Historical reference was CMat hap2; current reference is CMat hap1. Reference chromosome numbering and assembly sequences differ.

The original harmonization candidates comprise two CTlk BREAK_CANDIDATE records plus excluded fragments. The current table comprises one CTlk hap2 BREAK_CANDIDATE plus excluded CPla fragments. CBau, CLim and CMat have no composite candidates in either table. The increase from 2 to 67 evidence tables comes from the added broad diagnostic/recovery branch, not 65 new automatic cut candidates. Restore candidate-driven evidence scope.

### Concrete unsafe-location case (regression fixture, not a sample-specific rule)

Current CTlk hap2 scaffold_3 is 76,699,408 bp. Last-round AGP:
- object 1..30,202,138: scaffold_14, reverse orientation;
- gap 30,202,139..30,202,238;
- object 30,202,239..73,725,257: scaffold_7 source 4,306,001..47,829,019, forward;
- another gap, then scaffold_39.

The caller proposes 30,202,188. PAF records (zero-based, half-open query intervals) show:
- 29,987,160..30,623,478 aligns to reference scaffold_7 (consensus chr7), spanning the proposed cut;
- 30,623,615..31,397,238 also aligns to chr7;
- 32,161,237..33,665,002 aligns to chr7;
- 33,826,222..35,788,875 aligns to reference scaffold_12 (chr12).

These are substantial alignments with MAPQ 60. Their endpoints are evidence anchors, not an exact breakpoint estimate; repetitive/minor competing alignments and alignment operations must be considered.

Recovered older join J00000685 occupies current interval 33,004,681..33,004,781. Both local flanks were assigned chr7. Its presence inside the coarse 33..34.5 Mb transition interval does not justify cutting it. The source AGP places the following component h2tg000137l_1 at source 7,108,544..11,956,579. The observed switch is therefore consistent with lying within that component's current placement, rather than at either proposed N-gap. Corrections between stages prevent treating source arithmetic alone as verified provenance.

The main failure is in chimera_joins.py:
1. assign() reduces each last-round AGP component to one chromosome label, although a component can itself be an older scaffold containing multiple chromosome regions.
2. It credits the full PAF matching-base count to every overlapping component, rather than support clipped to that component. Repeated/overlapping alignments can inflate counts.
3. Unassigned components inherit labels from nearest neighbours by component index.
4. transitions() then places the transition at the component boundary. A nearby AGP gap is declared callable.
Thus the reported 51 bp distance is agreement with a constructed component-boundary transition, not agreement with the actual alignment transition. Coordinate integrity alone does not certify biological breakpoint location. Keep auto off until this is repaired.

## Minimal implementation sequence

### 1. One chromosome scope
Export the actual selected scaffold IDs already held in harmonize_names.py: qset_by_id for queries and chrom for the reference. Include selection method, flags, original ID, length, reference ID and coordinate-stage provenance. Do not reconstruct membership from chromosome names, lengths or voter status.
Pass this manifest through detection, evidence, automatic cuts and supplied-file cuts.
Remove chimera_min_span as a scope gate. Audit the separate 5 Mb member-footprint gate so it cannot silently exclude genuine small chromosomes; use alignment coverage relative to inferred chromosome lengths where a support criterion is necessary. Keep scope separate from evidence confidence. Ambiguous chromosome inference must be reported, not replaced by a new absolute threshold.

### 2. Restore selective evidence
Primary evidence is requested only for in-scope composite candidates and their supported junctions. No default all-scaffold sliding-window profile request. Put optional broad audits in supporting output, off by default for expensive profiles.
One diagnostic panel per candidate junction interval, with exact proposed gaps overlaid. Do not duplicate nearly identical whole-scaffold figures for a midpoint and several neighbouring gaps.
Use tidk alone; missing or malformed results are explicit unavailable/error states.

### 3. Derive transitions from alignment support
Use existing minimap2 alignments in the verified pre-finishing frame. Compute non-double-counted support over actual query intervals; use alignment operations where base-level placement is required and available. Never credit an entire record to a merely overlapping interval.
Retain ambiguity and mixed components; do not fill unknown labels to manufacture a transition.
Construct intervals between supported chromosome regions, then intersect with verified current or recovered AGP gaps. The existing coarse windows may help display evidence but cannot establish cut coordinates.
Zero compatible gaps -> REVIEW/no safe gap. Multiple compatible gaps -> REVIEW/ambiguous. A unique gap still needs direct flanking support consistent with the transition. A continuous confident same-chromosome alignment spanning a proposed gap is a conflict requiring review, not a majority-label override.
Reference self-assessment without alignments is unavailable, not clean.

### 4. One auditable decision contract
Each candidate junction records chromosome membership, distinct-individual support, sister-haplotype support, direct alignment anchors, competing support, interval bounds, matching gap IDs, assembly SHA256, evidence availability and decision reason.
Do not copy a scaffold-wide verdict onto every junction. Sister haplotypes are not independent biological corroboration. Recompute the existing conservative concordance policy per junction.
AUTO requires chromosome membership, sufficient unambiguous concordance, a uniquely justified gap, current provenance and valid resulting pieces. Diagnostic midpoints and in-component transitions are never automatic safe cuts.
Supplied files can explicitly approve reviewed biological decisions, but cannot bypass membership, sequence identity, interval validity or diagnostic-only restrictions. Reject unsupported instructions visibly.
Do not use a weak Hi-C ratio or telomere enrichment alone as permission to cut. Distinguish candidate, proposed action and actually applied action.

### 5. Validate before changing assemblies
Add synthetic fixtures for mixed last-round components, same-chromosome alignment across an upstream gap, overlapping alignment records, reversed components, edits/removal, multiple plausible gaps, chromosome membership, small true chromosomes, missing evidence, and per-junction sister/individual support.
Use sample-neutral fixture names and shuffled IDs to verify no name-specific rules. The observed 30.20 Mb proposal must be rejected without forcing a particular alternative cut or historical outcome.
Evaluate corrected decisions on the CURRENT cached assembly inputs with chimera_break=false. Confirm the current and historical labels can differ without forcing old calls.
Acceptance: primary report contains only scoped, justified candidates; every proposed automatic cut has explicit local support and a valid current gap; no excluded scaffold changes in either cutting mode; uncertainty remains visible.
Only after review should an auto/file-cut checkpoint run. Cutting changes downstream finishing inputs and may legitimately rerun finishing. Restoring the second Hi-C library is a later controlled input comparison.

## Cache and evidence collection

No new assembly is required to repair or test decision logic. Exporting the membership manifest may invalidate the combined harmonization task (which currently also computes PAFs); preserve or separate existing alignment inputs before promising those alignments remain cached. Assembly, mapping and scaffolding processes should not change.

Keep comparison helpers under scripts/comparisons and generated artifacts under comparisons. After the next evidence-only checkpoint, from the cluster project root:

~~~bash
bash gcl_genome_assembly/scripts/comparisons/collect_chimera_comparison.sh \
    tst/genome_assembly_store \
    genome_assembly
~~~

Upload comparisons/chimera-comparison-<timestamp>-reports.tar.gz and retain the alignment archive. Preserve baseline directories and work/cache. Check the completed log/trace against the report invocation rather than relying on the latest submitted job name.

## Remaining limits

This inspection establishes why the 30.20 Mb cut is unjustified. It does not establish an exact safe alternative, prove a biological translocation/misassembly, or prove where the old hap1 chimera went. Resolving cross-run sequence correspondence requires matching the actual assemblies, not comparing chromosome numbers. That is separate from correcting the demonstrable caller defect.


## Implementation checkpoint

The candidate-scope and location repair is implemented for a detection-only rerun. See [the rerun handoff](chimera-repair-rerun-handoff.md) for exact changes, remote checks, cache expectations and collection commands. Distinct-individual per-junction concordance remains a stated limit before enabling automatic cutting; this document describes the full repair target, not a claim that every step has been validated.
