# Assessment of the revised chimera/fusion plan

Date: 2026-10-03.
Reviewed: `chimera-fusion-assessment-plan_claudeAdditions.md`, revision 2.
Scope: assessment and recommendations only. No pipeline settings, cut policies,
or analysis jobs changed. The proposed teleost chromosome-number prior is
excluded from this review at the user's request. The user's 15-chromosome
species context must not become a target count that forces individual assemblies
to conform.

## Overall judgment

The revision substantially improves the plan's protection of genuine structural
variation. Its best additions are explicit heterozygous arrangements, junction-
level recurrence, retention of inter-scaffold Hi-C evidence, contig-end provenance,
unknown physical lengths of scaffold gaps, and calibration against real-data
controls. Adopt these.

The revision is not yet a specification for an automatic classifier. Several
quantitative predictions are idealized but written as biological expectations;
some claims of independence or exclusivity are too strong; and breakpoint
localization remains less explicit than event classification. Correct those
issues before implementing the scoring. Keep detection, structural assessment,
phasing assessment, breakpoint localization and edit authorization separate.

## Verified source findings

These are static findings in the current checkout, not runtime validation of the
cluster installation or proof that any real event has been cut.

| Revision claim | Assessment |
|---|---|
| The concordance heuristic can endanger private heterozygous arrangements | Confirmed as a policy risk: `harmonize_names.py` can assign BREAK_CANDIDATE when the sister is separate and other voters separate the pair. `break_chimeras.py` uses that verdict for auto selection. Coordinate, chromosome-scope and supported-gap checks also apply, so a private event is not invariably cut. Those checks do not establish that an adjacency is biologically wrong. |
| Recurrence can suppress detailed evidence | Confirmed: another individual's chromosome-pair recurrence can yield NOT_A_CANDIDATE; the reference emits an empty called table without a reference PAF. Add an audit route that does not depend on the vote or reference identity. |
| Recurrence is too coarse | Confirmed: concordance uses unordered reference-chromosome pairs, and scaffold-level aggregate verdicts can mix different junctions. Use oriented, anchor-defined adjacencies with explicit coordinate uncertainty. |
| Reference PAFs lack detailed alignment and alternatives | Confirmed in the harmonization commands: `--secondary=no`, with no `-c`, `-a` or `--cs`. `--eqx` alone does not request base alignment. Keep the coarse naming pass if useful, and add a separate detailed assessment pass. |
| Alternative Hi-C partners are discarded | Confirmed: `chimera_hic_pairs.py` skips cross-scaffold pairs in candidate evidence extraction. Full retained input pairs may still contain them; inspect before remapping. |
| Called-table headers disagree | Confirmed: the module's empty/stub header has 30 columns, while the Python schema has 39 before provenance stamping. Aggregate collection assumes a shared header. Centralize the schema and validate named fields for empty and populated inputs. |
| Report overstates cut validation | Confirmed in `r_scripts/generate_summary_report.R`: text claims Hi-C/telomere confirmation before cuts, while those diagnostics are not the auto-cut gate. Replace fixed prose with the actual decision/evidence record. |
| YaHS MAPQ settings imply weak filtering | This would be an incorrect inference. The current upstream pairtools filter already applies `hic_min_mapq` (default 30) and UU selection. YaHS's round-specific values are not the full filtering history. |

Stale cohort-specific comments and the unsupported assumption that a real fusion
must occur in both haplotypes should be removed. Line-ending normalization is
ordinary maintenance, not a biological prerequisite for the investigation.

## Scientific changes required

### 1. Replace fixed heterozygote predictions with conditional models

The table's lambda/2 and half-strength Hi-C values are useful illustrations only
under restrictive assumptions: balanced allelic sampling, comparable recovery,
the stated mapping representation, and negligible alternative/background signal.
They must not become cutoffs.

The claim that the remaining half of reads end at a real heterozygous junction
is especially misleading. Sequencing read ends occur throughout the region;
clipping depends on molecular geometry, mapping and the alternative structure.
Neither all reads clipping at a misjoin nor half clipping at a heterozygous
event is a universal prediction. Apparent alignments can also cross artificial
Ns; count aligned bases, deletions and supplementary placements explicitly.

For a read-opportunity model, define S as the total physical span needed to cover
both required anchors, including the anchor bases. A useful starting approximation
is `lambda = C_carrier / E[L] * E[(L-S)+] * q`, where q represents recovery through
mapping/filtering and C_carrier is coverage from the chromosome copies carrying
the adjacency. If C_carrier is already single-homologue coverage, do not halve it
again. This is a proposed calibrated model, not a validated likelihood for these
data. Estimate uncertainty and use local controls; repetitive-region coverage is
not uniform. At an unknown-length AGP gap, do not substitute the placeholder
length for S. Report a range or unassessable opportunity.

Likewise, a simple raw-contact mixture is `alpha * mu_joined + (1-alpha) *
mu_alternative`, with mapping and coverage effects included. It is not generally
half the normalized homozygous signal. Use empirical controls and simulations
to estimate its range.

The Hardy-Weinberg ratio in the revision is algebraically correct under its
assumptions. It illustrates why heterozygotes matter; unknown allele frequency,
population structure and selection prevent using it as this cohort's event prior.

### 2. Keep the improved Hi-C analysis, qualify its conclusions

Distance bands, repeat buffers, multiple orientations and alternative partners
are valuable. Remove the assertion that repeats can affect only junction-proximal
bins. Distant repeat copies, collapsed representations, incorrect placements and
contact normalization can influence wider regions; UU/high MAPQ is conditional
on the alignment reference. Conversely, do not dismiss strong distal contact
as repeat artifact without demonstrating that explanation.

Known Hi-C biases include sequence uniqueness, GC content and restriction-site
geometry ([Yaffe and Tanay](https://www.nature.com/articles/ng.947)). An all-end-pair
background helps, but does not automatically absorb those effects or telomere/
centromere clustering. Match or stratify controls by end type, available length,
repeat/mappability, coverage, library and, when reliably known, chromosome context.
Unannotated centromeres must not be treated as known control coordinates.

Adapt band sizes to informative arm length and coverage. Treat the suggested
0.2-1/1-5/5-20 Mb bands as examples, not general-purpose constants. Use block-aware
resampling and library-level sensitivity checks; individual Hi-C pixels are not
independent biological replicates. Distinguish sampling uncertainty from model
bias, and account for scanning many positions or choosing the strongest end pair.

A different statistic from the same scaffolding data is not held-out evidence.
For an upstream sensitivity experiment, building with one library and assessing
with the other is informative, but shared mapping biases remain. Do not imply a
newly computed distal statistic eliminates selection bias from the original join.

### 3. Add an explicit breakpoint-localization stage

The plan currently explains how to assess an event better than how to locate
its erroneous connection. Add this bounded procedure:

1. Define an uncertainty interval from distinctive chromosome-arm anchors and
   assembly provenance. Permit multiple disjoint intervals when placements are
   ambiguous; do not force an overlapping PAF midpoint into one breakpoint.
2. Enumerate actual source-component boundaries, gaps, graph branch/overlap
   boundaries, correction edits and clusters of informative read discontinuities
   within the interval. Inspect nearby source boundaries even when a coarse
   diagnostic position lies inside a component.
3. Evaluate both sides of each proposed breakend: orientation, unique-anchor
   support, compatible clipped/split reads, coherent alternative placements,
   copy context and distinct molecules supporting the current or alternative
   arrangement. Correct N-gap length uncertainty and repeat equivalence explicitly.
4. Report a breakpoint interval or a set of equivalent locations when the
   sequence cannot distinguish a single base. A lower Hi-C window identifies a
   structural neighborhood, not a nucleotide cut coordinate.
5. Authorize an edit only after a separate structural decision. For an artificial
   scaffolding gap, component boundaries often define the operational split. For
   contig errors, require a sequence-level location or preserve an unresolved
   bridge as a separately tracked record when that treatment is justified.
6. Validate all edits on a multi-junction scaffold jointly. Fixing one incorrect
   connection must not remove an independently supported connection or create
   a haplotype-inconsistent sequence path.

Outputs should separately record adjacency confidence, left/right breakend
uncertainty, phase confidence, physical-gap uncertainty and action eligibility.

### 4. Separate biological adjacency from correct haplotype reconstruction

Protecting a real rearrangement does not mean accepting every assembled version
of it. An adjacency can exist in the individual while the scaffold incorrectly
combines sequence from opposite homologues. Assess physical adjacency and local
phase consistency separately; a phase error can need correction without claiming
that the biological chromosomes are separate.

Rename the proposed deliverable to an **adjacency evidence/genotype table**.
Unphased diploid Hi-C cannot generally identify which named assembly haplotype
carries an event. Haplotype labels are not coordinated across chromosomes in
Hi-C phasing ([hifiasm README](https://github.com/chhylp123/hifiasm)). Report sample-
level support and phase-unresolved states; assign 0/1/2 copies or per-haplotype
states only when supported. Missing support is not a 0/0 genotype when opportunity
is low. A local adjacency table is not a complete karyotype.

The new genotyping section is less dependent on scaffolding, not reference-free.
Sniffles-style calling against a common reference requires valid read alignments
to that reference; existing per-assembly BAMs cannot simply be supplied unchanged.
Treat this as targeted corroboration where breakend geometry is callable, not a
mandatory whole-cohort remapping project. Reference orientation, gaps and repeat
resolution limit what any long-read SV caller can detect.

Recurrent supported linkage across the sampled individuals does not establish
population fixation. Matching repeat-mediated errors can recur at matching
locations. Evaluate each individual's own molecule evidence and distinguish
whole-chromosome events from local breakend associations.

### 5. Qualify read, telomere and orthologue evidence

Replace "only a single molecule" with multiple distinct, well-anchored molecules
or another resolving evidence source. One read can be artifactual; a uniquely
resolved graph path or appropriately informative longer-range data can also help.
Telomeric sequence does not invariably identify a former chromosome end, and
absence of an interstitial array does not exclude an old or non-telomeric fusion.

GFA A-lines are excellent for locating relevant reads, but are assembler-derived
placements, not an independent raw-read proof of every inferred overlap. Retrieve
and verify the molecule sequence where it is decisive.

BUSCO orthologue order is useful complementary arm evidence, not a guarantee of
local unique mappability or breakpoint identity. Restrict to reliable orthologue
assignments, check duplications, and use consistent lineage definitions and
coordinate stages. Reuse existing tables; do not make a new distant-outgroup
annotation exercise a prerequisite for resolving current joins.

### 6. Use added tools selectively

HMM-Flagger is a reasonable addition for diploid coverage/copy anomalies. Its
intended dual-assembly input and satellite-related coverage biases need attention
([HMM-Flagger documentation](https://github.com/mobinasri/flagger)). Independently
mapping the complete diploid reads to each haplotype and concatenating those BAMs
is not equivalent to competitive diploid mapping. Inspect input compatibility
before claiming this can reuse existing BAMs without remapping. A coverage-normal
misjoin can escape it; a coverage anomaly does not localize a cut by itself.

Pretext with the proposed tracks is useful for consistent manual review
([curationpretext](https://pipelines.tol.sanger.ac.uk/curationpretext)). Inspector,
CRAQ and optional repeat-aware remapping with Winnowmap can provide corroboration;
avoid making every tool mandatory for every case. Reanalyses of the same reads
may provide different information but must not be counted as independent votes.

### 7. Revise the calibration and edit policy

Retain-and-flag is a reasonable default for unresolved connections, provided the
uncertainty reaches the FASTA/AGP companion metadata and downstream analyses.
Neither retaining an unresolved scaffold nor provisionally separating it should
be described as proving a fusion or fission. Exclude claims that depend on the
uncertain adjacency, rather than discarding all valid local information from it.

Thinning half the positive evidence does not recreate a heterozygote: it lacks
the alternative arrangement's reads, background contacts and mapping competition.
Use mixtures of both arrangements and remap a selected subset end-to-end. Keep
simple thinning as a coverage stress test. Pure AGP liftover tests are useful but
miss mapping errors induced by a changed sequence representation. Hold out entire
loci/repeat families or individuals for validation rather than splitting highly
related windows between training and testing.

No false cuts among finite controls is a necessary target, not proof of zero
false-break risk. Report control counts and uncertainty alongside misjoin recall
and unresolved rate; an algorithm that never cuts must not pass on specificity
alone. Under an ideal independent-binomial model, zero errors among N controls
has an approximate 95% upper error bound of 3/N; correlated controls require
more careful treatment.

Remove only verified artificial AGP gap sequence during de-scaffolding. Preserve
internal sequence and ambiguous bases unless an explicitly audited edit changes
them; do not strip every N or require every internal cut to produce a third
record. Check interval identity/orientation and one-to-one sequence provenance,
not merely total non-N length.

Freeze cut-derived ends during adjudication and do not use subsequent extension
as independent validation. A permanent ban on extending all cut-derived ends is
too broad: a correctly restored chromosome end may later be extendable with
adequate unique support. Preserve pre-extension evidence and report that decision
separately.

Keep graph-specific input splitting as a separate, optional design decision.
The chromosome-splitting mode of Minigraph-Cactus can lose interchromosomal
relationships ([Cactus documentation](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md)).
That does not justify automatically splitting every retained arrangement at an
uncertain coordinate. Require a reversible coordinate manifest and a sequence-
retention audit; assess the exact configured Cactus version/mode before choosing
this workaround. It is not a prerequisite for the chimera assessment.

## Upstream prevention experiments

Freeze software versions and effective commands first. Hifiasm and YaHS conda
specifications are currently unpinned. Older hifiasm documentation describes
`-u` differently from current source help, which documents a 0/1 argument.
Verify the installed binary before specifying a rerun. Sources:
[current hifiasm option parser/help](https://raw.githubusercontent.com/chhylp123/hifiasm/master/CommandLines.cpp),
[older parameter reference](https://hifiasm.readthedocs.io/en/latest/parameter-reference.html).

| Suspected cause and trigger | Controlled comparison | Protection against losing real structure |
|---|---|---|
| Contig connection first appears during hifiasm post-joining | Compare current `hifiasm_post_join = 1` with post-joining disabled; for a compatible installed version this is `-u 0` | Check whether a questionable connection disappears while distinctive-anchor-supported connections remain. Fragmentation alone is not improvement. Keep dual scaffolding off in the initial comparison; it is already false. |
| Purging removes an alternative path or changes local copy representation | Compare current internal purge level 3 with a less aggressive setting such as level 1; separately test the external purge_dups stage when provenance implicates it | Assess both haplotypes, missing/duplicated sequence, reliable gene/k-mer completeness and phase. Review explicitly supplied similarity thresholds alongside purge level; changing level alone may not reproduce that level's documented default behavior. |
| Raw unitig ambiguity is driven by repeat-overlap limits | Only with graph evidence, test repeat-overlap retention options such as `-D`, `-N` or `--max-kocc` | These are not generic anti-chimera switches. They can alter correction/overlap computation, resource needs and cache validity. Preserve genuine alternate paths; defer unless simpler stage tests fail. |
| Hifiasm's Hi-C misjoin correction or telomere logic creates/removes the connection | Change one implicated option, such as `--l-msjoin` or telomere behavior, using a version-appropriate comparison | The misjoin option's length is a unitig-size eligibility threshold, not a universal confidence threshold. A true repeat-rich rearrangement can also be broken by aggressive correction. |
| YaHS creates joins using ambiguous or inconsistent links | Compare audited valid-pair inputs, retained MAPQ/alternative-placement strata, library-specific and combined scaffolding; use sample-matched library metadata | Current upstream MAPQ filtering is already 30 by default. Raising YaHS's own -q may have no effect with coordinate-sorted BAM; repair/standardize input semantics first. Blanket repeat masking or discarding all homolog-ambiguous pairs can destroy genuine linkage. |
| False joins first appear in coarse YaHS rounds | Replay the implicated invocation with its final problematic resolution(s) omitted or capped, keeping other inputs fixed | Large configured resolutions are not themselves a defect or a chromosome-length limit. YaHS adapts its usable range. Use actual logs/AGPs to associate round with resolution; score true connections lost as well as false joins removed. |
| Round 2 creates/reintroduces errors | Compare round-2 scaffold error correction enabled versus current disabled; compare with skipping round 2 only if source/correction history makes that a valid comparison | Do not disable round 2 unconditionally when upstream correction/decontamination needs re-scaffolding. Track both error correction and new joins. |
| Terminal repeat-driven YaHS joins | Test terminal-telomere protection if available in the installed version, first on annotated, confidently terminal ends | YaHS supports a telomere-motif option, but a real interstitial fusion array exposed by earlier fragmentation could be mistaken for a terminal end. Use as a controlled comparison, not a universal veto. |
| Inspector/correction creates a false replacement | Compare the exact replacement with its input; test less permissive correction or review-only treatment of unsupported complex edits | Do not disable all correction based on one case. Require raw support for whichever sequence is retained. |
| Finishing fabricates support or reconnects rejected pieces | Freeze implicated ends during review, increase uniquely anchored molecule requirements for extension after calibration, and retain gap-fill/extension edits | Finishing cannot explain pre-finishing errors, but can hide or modify them. Rejoining an explicitly rejected adjacency should require explicit authorization and new evidence. |

YaHS input-sort/MAPQ behavior, adaptive resolutions, correction switches and
telomere protection are documented in its
[official README](https://github.com/c-zhou/yahs). These are candidate experiments,
not recommendations to change all production defaults at once.

Use retained hifiasm corrected-read/overlap caches only with compatible versions,
identical input identity and settings for the cached computation. Graph-only
experiments may reuse them; changes to correction/overlap inputs require an
appropriate rebuild. Copy or stage immutable baseline caches into separate
experiment locations so one experiment cannot modify another's evidence.

For every comparison record: verified false joins removed, supported connections
lost, new questionable joins, copy/sequence completeness, haplotype consistency,
runtime and memory. Match outputs by sequence/anchors because scaffold names and
haplotype labels may change. Accept a new default only when it improves the
defined error outcomes without an unacceptable loss of real linkage on held-out
controls and unaffected assemblies. Do not optimize for 15 scaffolds, N50, or
one sample's known phenotype.

## Recommended implementation order

1. Correct the documented policy/report/schema defects; keep current production
   cutting disabled. Propose deactivating vote-only auto explicitly while a new
   authorization policy is developed. Retain reviewed, checksum-bound actions.
2. Complete one baseline assessment using provenance, oriented anchors, existing
   reads and library-specific inter/intra-scaffold contacts. Add the explicit
   breakpoint-localization and phase-confidence outputs.
3. Select a small first batch of upstream experiments from the first-appearance
   results: post-join/purge comparisons for hifiasm-origin cases; filtering,
   resolution or correction comparisons for YaHS-origin cases. Do not perform a
   full parameter grid or tune every stage before knowing the origin.
4. Resolve remaining contradictions with targeted competitive mappings or local
   reconstruction. Add cohort genotype calls only where the data can support
   them; do not make complete population genotyping a gate for every assembly edit.
5. Calibrate on realistic mixed-arrangement and mapping-aware controls, then
   integrate evidence-driven recommendations and finally conservative automation.

This preserves the revision's strongest biological protections while making the
investigation capable of identifying an actual erroneous connection, explaining
its origin, and preventing recurrence upstream.
