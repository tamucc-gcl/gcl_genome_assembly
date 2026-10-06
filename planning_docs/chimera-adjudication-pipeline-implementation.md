# Chimera evidence, adjudication and correction implementation

## CTlk result with H01 only

Use the isolated comparisons/ctlk-boundaries-1511285/H01/unjoined.fa as the one-break candidate. Hap1 scaffold_1 becomes 63,051,325 bp and 75,395,654 bp records; all original sequence, including 100 terminal gap Ns, is preserved. Hap1 has 363 rather than 362 records, 953,570,870 total bp, and N50 69,062,049 rather than 71,708,076. Hap2 is unchanged.

The left piece is chr9-associated. The right piece still has chr7/chr12 correspondence. Hap1 scaffold_5 remains the 74,765,583 bp chr4/14 composite; hap2 scaffold_5 remains the 83,430,767 bp chr7/12 composite; hap2 scaffold_7 remains the 82,879,935 bp chr4/14 composite. Separate chr4/chr7 pieces remain intact. Output scaffold ranks may change after finalization. Do not call this a fully resolved assembly or convert local-path support into biological-fusion confirmation.

## Integrate into the existing core, before sequence finishing

Extend workflows/chimera.nf and the existing CHIMERA_JOINS/CHIMERA_EVIDENCE/CHIMERA_SEQUENCE_CONTEXT/BREAK_CHIMERAS machinery. This is a permanent core workflow, not a dependency on investigation job directories, CTlk IDs, historical assemblies or user-written shell batches.

Stable input channels: pre-finishing FASTA + AGP + source-stage provenance; all current sister and peer assemblies with quality/scope metadata; HiFi reads or matching retained BAM/placements; per-library Hi-C pairs and manifests; retained native graphs when available. Evidence collection must not silently rebuild assembly inputs. Validate hashes, dictionaries and coordinate lifts before using measurements.

### 1. Detect candidate chromosome-scale discordance

Use long collinear chromosome-block correspondence across usable current assemblies, including sister comparison. Rank multi-block chromosome-sized scaffolds and proposed boundaries. Six haplotypes from three individuals count as three independent individuals. Pair labels in two sisters do not establish homologous junctions; compare actual reference endpoints and orientation. Preserve competing/repeat placements and fragmented pieces rather than forcing a single chromosome identity.

Detection must NOT be restricted to an existing gap or a join rare across haplotypes. The chr4/14 and chr7/12 composites in both CTlk sisters demonstrate why a haplotype-frequency gate can hide relevant cases. Presence of a discordant chromosome-scale structure creates an evidence candidate, not a cut authorization. Existing vote rules can rank candidates but cannot authorize edits.

Output candidates.tsv with stable IDs, literal assembly/scaffold interval, FASTA checksum, stage, source component/graph provenance, correspondence blocks and competing explanations.

### 2. Produce a bounded evidence packet per candidate and comparison controls

Collect lanes concurrently while batching expensive work per assembly: one Hi-C pairs scan; one reusable HiFi mapping/retrieval; batched peer-window mapping; graph neighborhoods only at implicated regions. Preserve native graph files at initial assembly, not via later historical hunts.

The human packet should include one HTML page per candidate with:

- Chromosome correspondence on both sisters and independent peers, with long collinear blocks, competing short anchors and uncertainty ranges.
- Each Hi-C library's contact map, cross/within-flank distance profiles, competing chromosome partners and matched controls. Explicit mappability and informative/non-informative labels.
- Local HiFi coverage, narrow spanning molecules and representative read alignments across a proposed boundary, plus competing placements. Broad-core zero bridges are not evidence of absence.
- Native graph neighborhood, overlap/orientation, alternate paths and read-placement support where available.
- An action proposal with exact coordinates, reasons for and against, remaining uncertainty, and projected resulting pieces.

Provide TSV/PAF/BED alongside HTML, an IGV session tied to exact reference names/checksum and BAMs, and a machine-readable provenance/status record. Figures should be generated from the measurements and show controls; an unavailable lane is not a zero-support lane. IGV files alone do not replace a readable decision summary.

### 3. Adjudicate; keep local-path and chromosome-scale decisions separate

Output decisions.tsv with separate candidate-level chromosome-scale status and boundary-level action. Use RETAIN, BREAK_PROBABLE_MISJOIN, UNJOIN_UNSUPPORTED and UNRESOLVED, plus local_path_support, localization_status, assay_informative and auto_eligible fields.

A locally supported sequence path can remain inside a chromosome-scale candidate requiring review. For CTlk, J02's tested seam is locally retained because three MAPQ60 molecules span it; the entire chr7/12 placement is not automatically approved. H01 is a reviewed conservative gap-unjoin proposal. UNRESOLVED must remain visible in final reports and must not become an implicit fusion call.

Suggest cuts only at literal verified gaps or defensibly localized internal boundaries. Never snap a broad transition midpoint to whichever gap is nearest without testing that gap. Internally supported boundaries must not be cut merely to reproduce peer chromosome labels.

### 4. Use one source-bound action file for both manual and automatic application

Proposed normalized schema (not the current legacy called-joins schema):

id, assembly, coordinate_stage, assessment_sha256, scaffold, action, cut_bp, gap_start, gap_end, decision_source, evidence_packet_id, reason.

Coordinates are 0-based; cut_bp splits sequence into [0,cut_bp) and [cut_bp,length). Gap unjoin policy should preserve every base. For the validated H01 candidate, gap_start=63051225, gap_end=63051325 and cut_bp=63051325; terminal gap Ns remain on the left piece. A gap-edge split must be validated by its actual all-N interval, not by requiring NN immediately across the selected cut position.

Manual mode uses only the rows explicitly supplied in the reviewed action file; it never applies all suggestions. Automatic mode produces the same action schema from eligible decisions, then uses the same guarded applicator. Evidence-only mode writes candidate/evidence/decision outputs and leaves all assemblies intact.

Required application guards: exact FASTA checksum and stage; unique scaffold names; integer coordinates; actual gap validation for unjoins; boundary localization/evidence identity for internal cuts; duplicate/conflicting-action rejection; whole-set minimum piece lengths for multiple cuts; no out-of-range edits. Never silently skip a malformed selected row. Report every applied and rejected action.

Manual internal cuts should be permitted through the explicit action schema, with mandatory review identity and localization rationale. They must not bypass coordinate/base-conservation guards. They need not masquerade as gap calls to pass the current breaker.

### 5. Conservative automatic policy

Ship evidence-only as default. Enable automatic gap unjoining separately from internal cutting. Initial auto eligibility requires verified gap localization, robust discordant chromosome blocks from multiple independent individuals, informative support measurements on both flanks and controls, an opposing/support-loss signal beyond peer-label disagreement, no informative molecule/graph support that contradicts the proposed action, and sequence/minimum-piece guards. Missing or non-informative evidence makes a candidate review-required, not auto-eligible.

Automatic internal cutting should remain disabled until a localized classifier is validated against intact junctions and known/accepted errors. This is a policy to implement and calibrate, not an already validated algorithm. Current nextflow.config and workflows/chimera.nf deliberately reject chimera_break=auto; do not remove that rejection without a decision-stage eligibility gate.

H01 should initially be a reviewed-file action, not automatically declared eligible from zero bridge/contact counts: immediate flank mappability and gap controls were imperfect. Human acceptance can authorize the conservative unjoin while the automatic rule remains stricter.

### 6. Apply once, verify, then continue the core

Write corrected FASTA, rewritten AGP/source-coordinate lift, per-piece sequence hashes, action audit and input/output accounting. Unchanged assemblies pass through. Reconstruct each original scaffold from its output pieces and compare checksum, including Ns; total-base equality alone does not detect substitutions, reorderings or balanced loss/duplication.

Use neutral piece identifiers, then rerun chromosome assignment/harmonization on the resulting pieces. Do not inherit one chromosome label or chromosome class from the parent. Preserve composite/unresolved annotations on pieces that still contain multiple blocks.

Run final QC and downstream finishing/graph construction from corrected outputs. Existing BAM/pairs coordinates are either validly lifted or explicitly labeled as pre-correction; never attach them to renamed/split FASTAs unchanged. Apply correction after the last scaffolding pass by default. If another scaffolding pass is requested, carry excluded-adjacency constraints and test for recreation of the disputed join.

## Existing-code gaps that must be fixed

py_scripts/break_chimeras.py currently accepts reviewed gap rows only, rejects auto, requires a cut strictly inside Ns, and derives chromosome names/classes from side labels. Therefore the validated H01 gap-edge split and its still-composite right piece require a schema/application update. Do not simply feed the new review worksheet to the legacy breaker or turn on its auto branch.

workflows/chimera.nf already separates detection/evidence/application and routes unselected assemblies through. Reuse that framework. Extend evidence collection from the tested helpers, removing sample/job hardcoding and adding candidate/control manifests. Replace the obsolete vote-to-cut concept with evidence-backed eligibility. Keep manual and auto application mechanically identical after selection.

## Implementation sequence

1. Upgrade the source-bound applicator, neutral naming, AGP/lift/audit outputs and H01 end-to-end regression, including partial unresolved composites.
2. Promote bounded evidence lanes into Nextflow modules with cached shared measurements and human HTML/IGV reports; preserve graph inputs and library identity.
3. Implement the decision/suggestion layer and reviewed-file round trip. Integrate before finishing and rerun assignment after edits.
4. Validate conservative automatic gap eligibility against controls and accepted corrections, then expose it as opt-in. Keep automatic internal cuts disabled until their own validation is adequate.

No pipeline application code or production settings were changed by this design document. The investigation scripts demonstrate measurements; they are not yet a general-purpose integrated adjudication process.
