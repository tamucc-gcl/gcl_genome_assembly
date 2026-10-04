# Complete assessment of suspected chromosome misassemblies and fusions

Date: 2026-10-03. Status: proposal for discussion; no new adjudication policy,
assembly edits, or analysis jobs are authorized by this document alone.

## Objective and completion criteria

For every questionable chromosome connection, establish (1) what sequence is
connected, (2) when that connection appeared, (3) whether the individual's raw
data support it over plausible alternatives, and (4) what assembly action is
justified. Assess an adjacency, not an entire named scaffold as one indivisible
case. A scaffold can contain both correct and incorrect connections.

Distinguish three decisions: structural interpretation, confidence, and safe
correction location. A convincing misassembly diagnosis does not automatically
identify an exact cut. Conversely, a precisely located AGP gap is not evidence
that the join is wrong. Supported adjacency establishes neither an ancestral
fusion mechanism nor a chromosome count by itself.

The final deliverable is a complete decision table, an evidence packet per
adjacency, and a proposed correction manifest. Unresolved cases must state the
specific observation that would resolve them and the interim assembly treatment.
No sample names or historical coordinates become production decision rules.

## Current evidence and why it has not settled the question

The current chromosome assignments identify substantial arms, but their PAF
records lack base-level CIGAR information. Overlapping chain spans make the
diagnostic midpoints poor proxies for actual junctions. HiFi tiling or a read
across one midpoint may lie entirely within shared repeat sequence; it does not
necessarily link the two distinctive arms. Several wide assignment intervals
are longer than typical reads, so failure to span their full width is not itself
evidence of error.

The existing assignment review also shows that the largest composite contains
parts of chromosomes, rather than three complete chromosomes. Some associations
recur across both haplotypes at coarse resolution. These observations change the
hypotheses worth testing, but establish neither correctness nor a biological
fusion. The same reads, repeats, reference and assembly algorithms can generate
correlated outcomes in both haplotypes.

The investigation must distinguish:

- A false scaffolder connection between otherwise correct contigs.
- A connection already present in a source contig: repeat collapse, an incorrect
  graph traversal, or a haplotype switch/other contig error.
- A connection introduced or altered by correction or sequence replacement.
- A correct adjacency obscured by chromosome naming, ambiguous alignment, or a
  problem in a comparison assembly.
- A genuine structural arrangement in this individual, potentially a fusion,
  translocation, or other rearrangement.

## 1. Freeze inputs and retain the evidence

Inventory existing published files and retained task outputs before rerunning
anything. Record checksums, sample/read-library identity, tool versions, effective
commands, assembly stage and coordinate conventions. Test links for readability;
a filename alone is not proof that the target survived scratch cleanup.

Required retained products:

| Stage | Products to locate and retain |
|---|---|
| Hifiasm | Haplotype contig GFAs/FASTAs, raw and processed unitig GFAs, available read-placement annotations, low-quality intervals and logs |
| Contig processing | FASTAs before/after purge, filtering, Inspector correction and decontamination; edit, deletion and coordinate records |
| Each YaHS invocation | Exact input FASTA/index, final AGP/FASTA, all generated intermediate AGPs, BIN, effective parameters and logs |
| Between scaffolding rounds | Corrected/decontaminated FASTAs and edit records linking the two rounds |
| Assessment/finishing | Assessment FASTA and mappings; proposed/applied cuts; gap-fill and extension edits; final sequence/name maps |
| Read evidence | HiFi BAM/index and its reference identity; Hi-C pair/BAM evidence with library identity, filtering and deduplication provenance |

The current module source explicitly emits haplotype contig GFAs and final YaHS
AGPs. It does not explicitly emit all raw/processed unitig graphs or intermediate
YaHS AGPs. This is an output-contract observation, not proof that those files are
absent on the cluster. First inventory them; then propose explicit output and
publication declarations for any needed products not reliably retained.

Hifiasm's unitig graphs and read annotations can expose information lost in a
contig FASTA ([hifiasm output documentation](https://hifiasm.readthedocs.io/en/latest/interpreting-output.html)).
YaHS produces initial-correction and intermediate-round AGPs as well as its
final AGP ([YaHS documentation](https://github.com/c-zhou/yahs)).

If required historical artifacts are missing, assess whether final AGPs and
retained sequences already answer the origin question. Reproduce only the
necessary task if they do not. Changing output declarations cannot recover files
already removed from an old task. Keep reproduction outputs separate and verify
their relationship to the original run before treating them as historical proof.

**Deliverable:** an availability/provenance manifest and an explicit list of
genuinely missing evidence. No blanket assembly rerun.

## 2. Reconstruct every potentially relevant connection

Audit all joins on the suspect composite scaffolds, including joins outside the
current seven diagnostic intervals. Also screen all chromosome-scale scaffolds
for previously unflagged structural discontinuities using the existing inferred
chromosome set, not a new absolute size threshold. Smaller pieces remain useful
for locating missing/duplicated sequence and alternative partners.

Trace the sequence through this chain:

`hifiasm -> contig filtering/purging/correction/decontamination -> YaHS round 1
-> scaffold correction/decontamination -> YaHS round 2 -> assessment/cuts
-> gap filling -> telomere extension -> final naming`

For each connection, identify its first observed stage and whether it is an AGP
join, an internal source-contig connection, an edit boundary, or still unknown.
Cover hifiasm's actual Hi-C/phasing/scaffolding settings as well as YaHS. Trace
source orientation and component subdivisions. Validate coordinate transfer by
sequence alignment where correction changed bases, lengths or contig structure;
do not propagate coordinates arithmetically through an unverified edit.

Inspect every original gap, later filled gap, source-component boundary,
correction/replacement interval and relevant graph branch within the suspect
region. A midpoint falling inside a W component does not clear that component.
Post-assessment finishing cannot have caused a structure already present at
assessment, but its effects must be followed into the final assembly.

**Deliverable:** an adjacency ledger with first-appearance stage, source pieces,
orientation, current coordinates, uncertainty interval and lineage confidence.

## 3. Verify chromosome-arm identity and account for all sequence

Produce base-resolved alignments through the uncertain regions, retaining
competing placements. Minimap2 requires an alignment-producing option such as
`-c`/`--cs` or `-a`; approximate mapping spans alone are insufficient for this
purpose ([minimap2 FAQ](https://github.com/lh3/minimap2/blob/master/FAQ.md)).

Compare against multiple suitable peer assemblies, both haplotypes, and all
other pieces containing the implicated arms. Select peers by assembly/evidence
quality, not hard-coded sample identity. Repeat assignments with an alternative
reference and with the original reference excluded. Fragmented assemblies can
contribute local sequence evidence even if they cannot vote on chromosome order.

For each arm, report orientation, distinctively aligned sequence, ambiguous
sequence, missing intervals, duplicated intervals and all alternative locations.
Use union coverage rather than summing overlapping alignments. Determine whether
the proposed structure combines whole chromosomes, partial arms, or a short
repeat-related excursion. Check reciprocal products and copy balance across the
individual; absence of a reciprocal product is contextual evidence, not a veto
on nonreciprocal biology.

**Deliverable:** an arm-content map and competing structural hypotheses. Close
assignment-only false alarms before attempting assembly cuts.

## 4. Establish whether the available molecules can resolve each junction

Annotate candidate regions with existing repeat annotations, or targeted
[RepeatMasker](https://www.repeatmasker.org/RepeatMasker/) and
[Tandem Repeats Finder](https://tandem.bu.edu/trf/trf.html) analyses if needed.
Retain tidk telomere evidence. Add self/alternative alignments and k-mer copy
information to identify distinctive anchors. Unmasked sequence is not
automatically unique; a telomeric motif is not a fusion diagnosis.

Locate the nearest reliable anchors on both arms. Report the distance between
them and compare it with actual HiFi read lengths and local, haplotype-aware
coverage. Distinguish a unique locus shared by both homologues from a sequence
repeated at unrelated loci. Parental haplotype markers must not be assumed when
parental data are unavailable.

Classify evidence opportunity before interpreting missing support:

- **Directly testable:** available reads should connect the anchors.
- **Potentially graph-resolvable:** no single read spans, but informative
  variants/overlaps may distinguish competing paths through the region.
- **Not resolvable with current molecules:** the repeat exceeds informative
  reach and alternative paths remain indistinguishable.

Calibrate spanning opportunity with matched well-supported regions of similar
repeat content, coverage and anchor separation. Report expected opportunity and
uncertainty rather than treating zero observed spans as equally meaningful
everywhere.

**Deliverable:** anchor locations, competing copies, repeat structure and an
explicit statement of what existing data can and cannot test.

## 5. Evaluate raw-molecule and assembly-graph support

Reuse retained HiFi alignments when their reference matches the assessment
sequence. Gather all reported placements of candidate read names, including
secondary/supplementary alignments and clipping. Existing aligner output can
omit alternative placements; absence of an emitted alternative is not proof
of uniqueness.

Where needed, align recruited reads competitively against both haplotypes and
whole-assembly competing loci. Never estimate uniqueness using only an isolated
candidate FASTA. Recruitment must include reads anchored on either side and
relevant clipped/unmapped or competing-copy reads, not only reads already
spanning the disputed join. Document residual recruitment bias.

Count distinct molecules that connect the distinctive anchors in a coherent
alignment, with error, clipping and alternative-placement information. Separate
these from repeat-only spans and local tiling. In diploid competitive mappings,
low MAPQ caused by homologous copies must not automatically be interpreted as
an unrelated repeat or lack of support. Evaluate allele consistency where
informative variants exist, and avoid double-counting supplementary records.

Use existing Inspector results and propose CRAQ in diagnostic-only mode on
compatible full BAMs. CRAQ supports supplied mappings and identifies clipping/
coverage patterns; its default classifications need checking for these phased
assemblies. Agreement between tools using the same reads is corroboration, not
independent experimental replication ([Inspector](https://github.com/ChongLab/Inspector),
[CRAQ](https://github.com/JiaoLaboratory/CRAQ),
[CRAQ paper](https://www.nature.com/articles/s41467-023-42336-w)).

Inspect the corresponding raw and processed hifiasm graph neighborhoods, using
[Bandage](https://rrwick.github.io/Bandage/) for review. Identify branches,
repeat nodes, coverage changes, read-supported paths and possible phase switches.
An assembler-selected path or a chain of repeat-overlapping reads is not itself
proof of a unique chromosome connection.

**Deliverable:** evidence for the current adjacency and for each plausible
alternative, plus a clear explanation of any unresolved repeat ambiguity.

## 6. Test Hi-C support at the arm and chromosome scales

Assess each Hi-C library separately before combining them. Use appropriate
enzyme/mappability information, duplicate handling and comparable distance
scales. Inspect both raw and normalized contact maps; low coverage can distort
normalization. Audit the filters actually applied upstream of YaHS, including
MAPQ, pair selection and secondary/supplementary handling, rather than inferring
them solely from requested YaHS parameters.

For every proposed connection, assess contacts between distinctive arm sequence,
distance-dependent continuity across the boundary, abrupt block separation,
orientation consistency and stronger alternative chromosome partners. Compare
with matched ordinary joins and repeat-rich controls. Examine the whole contact
neighborhood, not only the lowest window or the plotted midpoint.

Use a standard viewer such as Juicebox; YaHS supplies a supported route for
generating these maps ([YaHS contact-map documentation](https://github.com/c-zhou/yahs#generate-hic-contact-maps)).
Preserve library identity and exact coordinate provenance. Lift positions only
through verified mappings; changed sequence may require remapping. A library
already used for hifiasm/YaHS is not held-out validation. Agreement across two
libraries helps assess robustness but can share repeat-mapping biases.

**Deliverable:** a consistent set of arm-scale and local maps, library-specific
support summaries, and explicit support for alternative partners where present.

## 7. Resolve contradictions with one targeted follow-up batch

Only unresolved cases proceed. Each job must state which competing hypotheses
it can distinguish and what result would change the decision.

| Uncertainty | Targeted follow-up |
|---|---|
| Suspect YaHS join between correct contigs | Replay the relevant scaffolding stage from verified inputs, comparing the current join with a split representation and evidence-supported alternatives; assess libraries separately if informative |
| Connection already inside a source contig | Inspect graph/read paths and, if informative, locally reassemble reads recruited from both anchors and alternative copies |
| Apparent connection created during correction | Compare input/output sequence and raw support for the original and replacement paths |
| Ambiguous chromosome assignment | Obtain better base-resolved multi-peer comparisons; avoid changing assembly sequence to satisfy a label |
| Repeats longer than informative molecules | Specify additional data that would resolve the ambiguity; do not keep tuning mapping thresholds |

Score alternatives with the same data and comparable representation. Do not
choose by N50, chromosome-size appearance, or ability to force agreement with
the reference. Local reassembly must include competing reads/copies and account
for diploidy; otherwise it can reproduce the same biased path. A second
assembler using the same reads adds algorithmic corroboration, not new molecules.

If required, longer reads spanning distinctive anchors, optical maps with
sufficiently informative labels, or chromosome-scale FISH/karyotyping can
address different remaining questions. Targeted PCR is useful only for suitably
sized, uniquely anchored junctions. None is an automatic requirement for every
case. Distinguish establishing adjacency from establishing chromosome count or
the evolutionary direction of a rearrangement.

## 8. Adjudicate, propose edits and validate the result

| Outcome | Required reasoning | Proposed action |
|---|---|---|
| Supported misassembly | Discriminating evidence contradicts the current adjacency, with adequate evidence opportunity and a coherent alternative/explanation | Break or repair at a separately verified location |
| Supported adjacency / rearrangement | Distinctive sequence linkage and compatible wider structure survive repeat, placement and provenance checks | Retain; reserve "fusion" for stronger mechanistic/chromosomal evidence |
| Assignment/reference artifact | Structural discrepancy disappears after correcting the comparison or membership interpretation | Correct assignment/reporting; do not cut sequence |
| Unresolved | Available data cannot distinguish repeat/phase/structural alternatives | Record uncertainty and choose a transparent interim assembly treatment |

No universal number of reads or fixed MAPQ cutoff establishes these outcomes.
Require multiple independent molecules when making a molecule-based claim, and
calibrate sufficiency against coverage, ambiguity, read length and matched
controls. Weigh evidence by what it measures; do not take a simple majority vote
of correlated tools or plots.

For unresolved scaffolder-created joins, conservative de-scaffolding can be
considered as an assembly-delivery choice, explicitly labelled provisional.
It is not proof of a misassembly. Internal contig cuts need stronger localization
and a specific rationale. Never manufacture a gap or use an assignment midpoint
just to make a cut executable. Preserve uncertain sequence rather than silently
deleting it.

The edit manifest must bind every action to an input checksum, coordinates,
orientation, coordinate convention, reason, evidence and confidence. Rebuild
from the earliest necessary stage and rerun affected downstream stages only.
Confirm expected sequence retention, changed junctions, chromosome-arm accounting,
read support and Hi-C structure after editing. A successful pipeline exit alone
does not validate the correction.

## Bounded implementation and review schedule

**Batch A: one complete assessment from existing evidence.** Inventory/retention,
lineage reconstruction, base-resolved arm assignments, repeat/anchor assessment,
reused-read diagnostics and library-specific Hi-C review. Produce the decision
table for all cases together. This replaces serial midpoint investigations.

**Batch B: one targeted ambiguity-resolution batch.** Run only the selected
competitive mapping, local graph/reassembly or scaffolding comparisons that can
change an unresolved decision. Review all resulting cases together.

**Decision checkpoint.** Close each case with an adjudication and proposed action,
or a stated evidence limitation and specific new-data requirement. Reopen only
for new discriminating evidence or a demonstrated analysis defect, not merely
because the preferred answer has not emerged.

Most inventory, lineage and retained-alignment summaries should be relatively
lightweight. Base-resolved comparisons and competitive remapping are moderate
tasks. Whole-sample remapping, assembly reproduction, repeat-library construction
and new sequencing are escalation tasks, not unconditional first steps. Estimate
resources from a representative task before launching expensive batches.

## Eventual general-purpose pipeline implementation

Separate immutable assembly/read provenance and reusable mappings from evidence
extraction, adjudication and editing. Expanding a report should not invalidate
the mapping merely because a candidate table gained a column.

Build custom code for coordinate lineage, arm accounting, molecule/anchor
summaries, report integration and the edit audit. Use published aligners,
assembly assessors, repeat tools and graph/contact viewers for their established
roles. Initially emit recommendations for review; expand automatic breaking
only after calibration.

Before enabling automation, test planted misjoins, repeat-derived false spans,
supported heterozygous arrangements, phase switches, low coverage, reversed
components, filled gaps and changed coordinates. Include matched real supported
joins and independently assessed historical examples. Assess false breaks as
well as missed errors, stratified by repeat content and evidence opportunity.
Simulations test mechanics and expected behavior; they do not make ambiguous
real events ground truth.

One concise evidence packet per adjacency should contain a lineage/arm diagram,
the decisive read/graph evidence, local plus arm-scale Hi-C views, and the decision
with alternatives rejected or unresolved. Keep detailed tables downloadable.
Use ggplot2/tidyverse/patchwork for new R figures. Comparison scripts belong under
`scripts/comparisons/`; outputs under `comparisons/{TESTID}/`; preserve baselines
under `comparisons/baselines/`. No further assembly changes are part of this plan.
