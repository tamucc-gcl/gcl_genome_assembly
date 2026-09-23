# Assembly and pangenome refocus — agreed implementation plan

**Status:** decisions agreed with the user; implementation and experiments not yet performed.  
**Date:** 23 September 2026.  
**Supersedes:** the recommendations in the [earlier revised plan][previous] where they conflict with the choices recorded here. The [static review][review] remains the detailed evidence ledger.

## 1. Purpose and scope

Produce reliable assemblies and a meaningful report for a single short-read sample, multiple HiFi/Hi-C samples of one species, or mixed input types across several species. Automatically run applicable analyses, announce skipped samples, and avoid silently changing the intended biological cohort.

Streamline the maintained workflow. Prefer published tools over bespoke analysis; replace overlapping implementations rather than keeping both indefinitely. Preserve specifically requested results and important structural-variant capabilities.

**Outside scope:**

- Annotation, including repeat/gene annotation: handled by the user's separate pipeline.
- Individual-level population analysis and variant filtering for that purpose.
- Mapping/genotyping new samples against the pangenome.
- TellSeq and ONT implementation: future development only. Assembly-based eligibility should accommodate future technologies without technology-specific exceptions.
- Repairing the nonfunctional optional BlobTools/decontamination-evidence branch: keep unsupported/disabled and reserve for a likely separate rebuild. Existing assembly decontamination is not thereby disabled.
- A standalone pangenome/analysis-only entry point: future enhancement. Use normal full-pipeline resume now.

This planning work used static source inspection and documentation. No pipeline commands, analysis scripts, tests or experiments were executed, and no pipeline source was changed.

## 2. Input handling and assembly routing

### Sample-sheet behavior

Process eligible rows and announce every skipped sample with its reason in the startup summary and final report. A Hi-C-only row may remain in the sheet: it is an expected assembly skip because it lacks supported assembly reads.

Distinguish intentional exclusions, unsupported/malformed rows, missing inputs and actual execution failures. Ambiguous identities such as duplicate sample IDs must never silently assign files to the wrong sample; the implementation must handle them explicitly.

Use proper delimited-file parsing; preserve taxid, individual, assembly and haplotype identities independently of filenames. Validate normalized graph names for collisions. Propagate supported per-row overrides, including genome-size and deduplication choices.

### Effective assembly method

Select the method automatically with an explicit override. Record the effective method and representation in metadata and the report.

| Inputs | Route |
|---|---|
| HiFi + Hi-C | HiFi assembly with applicable Hi-C phasing/scaffolding |
| HiFi alone | HiFi partially phased assembly |
| HiFi + short reads, optionally Hi-C | HiFi route; supplementary reads remain separate available inputs |
| Short reads alone | Short-read assembly route |
| Hi-C alone | Announced assembly skip |

Downstream routing follows what was assembled, not simply whether short reads are present. Unsupported combinations must not silently take a different route.

## 3. Assembly processing policies

### Organelle recovery and screening

Unsuccessful organelle recovery must not block the main assembly. Report recovery separately from screening. Use an appropriate available screening reference where supported; otherwise state that screening is incomplete. Never claim guaranteed nuclear-only output from an incomplete screen.

Automatically separate confidently identified organelle contigs, preserving them and their removal evidence. Retain and flag ambiguous/partial matches. Review is optional and never a prerequisite to completing the assembly. Do not automatically excise a partial organelle-like interval from a nuclear contig. Distinguish genuine tool errors from a valid no-recovery outcome.

Prefer supported tool evidence and transparent rules over building a new organelle classifier. Nuclear mitochondrial insertions are a real reason to avoid blanket match-based removal; MitoHiFi explicitly addresses their distinction from mitochondrial contigs. [MitoHiFi paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10354987/)

### Purging and short-read conditioning

Preserve current defaults:

- Hifiasm internal purging remains enabled with the current configuration.
- Additional purge_dups remains off by default.
- Redundans reduction, scaffolding and gap closing remain enabled on the appropriate short-read route.
- Repair per-sample dedup controls so requested choices actually govern execution.
- dedup=none means no additional reduction; it does not disable hifiasm's internal assembly behavior.
- Keep Redundans scaffolding/gap closing independently controllable.
- Preserve removed sequences and record effective processing.

No new biological purging strategy or threshold change is agreed here. Targeted evaluation can inform a future change.

### Correction and finishing

Retain current defaults: Inspector contig/scaffold correction enabled where applicable, gap filling on the scaffolded HiFi route, Teloclip extension enabled where applicable, and Pilon disabled. Make the steps explicitly controllable; fix concrete routing and coordinate defects without bundling in a change of scientific defaults.

### QC

Retain none, final_only and all_stages, with final_only the default. Do not add new routine QC checkpoints.

This does not waive output-integrity requirements: missing/corrupt required files, broken name maps or lost samples must fail explicitly rather than masquerade as low biological quality.

Report biological-quality concerns without introducing new automatic exclusions. Explicit sample exclusions remain available. With QC disabled, say unassessed rather than passed.

### Harmonization and chimera handling

Assess chromosome-scale suitability from the assemblies themselves. Use the existing assembly-based definition as the starting point and verify its implementation.

- Suitable assemblies can vote in reference/consensus selection.
- Fragmented assemblies can be nonvoting participants and still enter the pangenome if otherwise eligible.
- HiFi-only assemblies will commonly be nonvoters because of fragmentation, not because of a blanket technology exclusion.
- Comparative harmonization bypasses cleanly when suitable comparison inputs are unavailable; finalization and identity tracking still run.
- Chimera detection remains enabled where applicable.
- Cutting remains opt-in, using the existing explicit automatic-call or supplied-breakpoint modes.
- No mandatory review checkpoint is introduced.

Before any cutting, correct species pairing, short-read bypass loss and absent/invalid name-map handling. Breakpoints must refer to the exact assembly being modified. Old AGPs do not automatically remain valid after correction, decontamination, gap filling or telomere extension. Maintain a validated coordinate transformation or remap the necessary evidence; choose the smallest correct implementation, not the cheapest incorrect one.

## 4. Pangenome admission and failure policy

Pangenomes run automatically per species unless globally disabled.

Requirements:

1. At least two eligible biological individuals of the same species. One diploid's two assemblies do not satisfy this threshold.
2. At least one suitable chromosome-scale reference.
3. Eligible biological representations and valid input files/identities.

Eligible representations include Hi-C/trio-phased haplotypes, HiFi-only partially phased haplotype pairs, and genuine haploid assemblies subject to suitability. Exclude collapsed diploid assemblies from automatic admission. Never relabel primary/alternate as two equivalent haplotypes. Short reads alone are not a disqualifier for a genuine haploid assembly.

Hifiasm documents its HiFi-only default as two partially phased assemblies, distinct from primary/alternate output. Preserve that distinction in metadata. [Hifiasm documentation](https://hifiasm.readthedocs.io/en/latest/pa-assembly.html)

If no suitable reference exists, skip that species' graph and report why. Cactus explicitly emphasizes a chromosome-scale reference; this does not require all nonreference inputs to be chromosome-scale. [Cactus reference guidance](https://github.com/ComparativeGenomicsToolkit/cactus/blob/v3.2.1/doc/pangenome.md#reference-sample)

Do not add new automatic exclusions for unusual size, completeness or duplication warnings. Report limitations and honor explicit exclusions.

**Actual execution failure:** continue independent samples, but withhold the affected species' pangenome by default when a required assembly step fails. Unaffected species may finish. Planned eligibility skips are different from an unplanned failed member of the intended cohort. Implement explicit failure/completeness accounting rather than merely ignoring process errors.

## 5. Graph roles and variant products

| Product | Agreed role |
|---|---|
| CLIP | Main biological graph and main reporting view |
| GREF(CLIP) | Coordinate extension and sole default variant catalog |
| FULL | Compact supporting QC and clipping-sensitivity role |
| Conventional assembly-reference VCF | Opt-in downstream deliverable, not a parallel default analysis |

Avoid a parallel FULL report suite. Produce corresponding machine-readable sharing summaries and compact clipping comparisons; do not duplicate every figure.

FULL means the unclipped construction for the selected inputs/settings, not proof that every assembly base is represented. Distinguish removed sequence from changes in sharing-category totals.

Gref's synthetic paths are coordinate objects, never extra biological samples. Validate the installed build, biological path allowlist, matching coordinate FASTA and fragment provenance. Current Cactus documentation describes gref as experimental and as adding paths without new nodes or edges. It also describes parent/nonreference coordinate relationships; counts and allele-length sums must not be mistaken for independent events or unique affected bp. [Cactus gref guidance](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md#graph-reference-paths---gref)

**Default variant representation:** supported standard vcfbub-processed gref catalog. Fine decomposition is requested explicitly or triggered by an enabled consumer that needs it. Follow the pinned tool's supported raw-to-wave route rather than assuming decomposition of the already-filtered parent catalog is equivalent. Preserve source IDs and invalidate incompatible graph fields. [Cactus VCF guidance](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md#vcf-output)

If Cactus necessarily produces conventional VCFs, they can remain construction artifacts without receiving duplicate analysis/reporting.

## 6. Sharing, growth and similarity

### Panacus as the preferred engine

Replace overlapping custom counting/growth/report calculations with Panacus where it meets the required outputs. Retain thin integration/aggregation only where necessary; do not retire required capabilities before feature coverage is demonstrated. [Panacus paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11665632/)

Required outputs:

- Sharing categories per chromosome.
- Sharing categories per haplotype.
- A chromosome × haplotype table with bp and percentages.
- The full exact sharing-count distribution.
- Descriptive growth/core curves.
- Main CLIP figures plus compact FULL sensitivity results.

Panacus supports regional subsets and sample/haplotype grouping, but exact cross-tabulation/export support must be verified in the pinned release. [Panacus documentation](https://github.com/codialab/panacus/wiki/YAML-report-file)

Use biological haplotypes as the default sharing unit; individual grouping is available on request. Define sharing across the species cohort, then summarize by chromosome/haplotype. Do not subset to one haplotype and thereby classify everything as private.

### Categories and units

Use expanded categories with configurable thresholds. Proposed defaults presented during the decisions:

| Category | Default |
|---|---|
| Private | Exactly one eligible haplotype |
| Core | 100% |
| Soft-core | ≥95% and <100% |
| Shell | ≥15% and <95%, excluding private |
| Cloud/rare | <15%, excluding private |

Assign private first to keep bins disjoint. Display actual count boundaries, retain empty categories when warranted, and validate threshold ordering. Label core as strict only when its threshold is 100%. These are reporting conventions, not universal biological boundaries.

Retain both count and denominator. Distinguish graph sequence bp from haplotype-spelled bp and repeated traversals. Define chromosome assignment, unplaced sequence, and missing/sex-specific chromosome handling explicitly during implementation. Do not silently discard those cases or equate missing representation with confirmed absence.

### Statistical interpretation

Use published expected growth calculations. Retire the custom base-independence uncertainty band and unsupported open/closed threshold as part of replacing overlapping bespoke analysis; do not introduce a new extrapolation suite by default.

HWE does permit independent allele draws at a locus under the relevant random-mating model. The earlier blanket statement that homologues must always be treated as dependent was too strong. Single-locus HWE alone does not establish genome-wide independence across linked sequence or eliminate shared technical errors. Descriptive observed-cohort sharing/expected curves do not require a claim of population sampling independence. Any future uncertainty or extrapolation must state and validate its sampling assumptions.

### Haplotype similarity views

Retain heatmap, ordination and tree views, preferably from one Panacus similarity matrix. Verify the matrix export, distance definition, missingness handling and resource cost before replacing ODGI similarity.

Use established implementations for PCoA of a distance matrix and a neighbor-joining similarity tree. Label them as descriptive sequence-content comparisons, not population-genetic or evolutionary conclusions. Do not silently turn missing distances into zero or alter eigenvalues/branches solely to obtain a plot.

If Panacus cannot provide a suitable matrix, evaluate a published alternative before retaining bespoke computation. Individual-level variant PCA is out of scope.

## 7. Required structural-variant capability — experimental milestone

**Inversions, duplications and translocations are required general capabilities.** Tool selection remains unresolved.

Do not remove existing bespoke analysis until targeted experiments establish an adequate replacement, or support retaining a clearly bounded part. This is not a commitment to maintain all candidate tools permanently.

At the relevant rewrite stage:

1. Inventory the current classifier, untangle/rearrangement and inversion-rescue outputs and their intended biological questions.
2. Select bounded known-answer examples plus representative real loci, including nonreference insertions, nested variation, inverted duplications and ambiguous repeats.
3. Compare sensitivity, false positives, class definitions and genotype/event correspondence.
4. Include both chromosome-scale and fragmented/partially phased inputs. State each method's eligibility.
5. Measure alignment plus analysis cost, not only the classifier's runtime.
6. Determine how results relate to gref coordinates and prevent double-counting across callsets.
7. Choose a replacement/retention strategy and document capabilities lost or unresolved.

Candidates investigated:

- bcftools/vcflib: standard sequence-allele statistics, not a complete mechanism classifier.
- SVIM-asm: assembly-alignment SV calling; candidate for broader input representations. [Documentation](https://github.com/eldariont/svim-asm)
- SyRI: chromosome-scale synteny/rearrangement analysis; orientation and input requirements matter. [Documentation](https://github.com/schneebergerlab/syri)
- Pantree: graph-native alternative variant representation; not a drop-in classifier for the existing gref catalog. [Documentation](https://github.com/oclb/pantree)
- Truvari: useful for compatible-callset comparison, not a universal replacement caller.

No tool benchmark or selection has been completed. Avoid an expensive full-cohort experiment before bounded tests justify it. If none provides sufficient coverage, retaining bespoke analysis requires explicit tests, honest class definitions and documented limits.

## 8. Optional and retired analyses

| Component | Decision |
|---|---|
| Private-sequence mapping/k-mer evidence | Retain temporarily, explicitly enabled, for the specifically requested run; outside the default report. Reassess afterward |
| Self-mapping QC | Opt-in, bounded reproducible read sampling; assess value/cost before any default promotion |
| New-sample mapping/genotyping | Outside scope |
| Mapping indexes | Generate as required by enabled QC or explicit export; identify compatible build/settings |
| Progressive-growth branch | Remove from maintained pipeline |
| Whole-genome 2D visualization | Remove from maintained pipeline; retain useful graph artifacts for external tools |
| Annotation/population analysis | Outside scope |

Before using private evidence, repair base-alignment/identity handling, search sensitivity interpretation, pooled versus per-assembly coverage definitions, no-read-support classification, k-mer schemas, private/control keys and calibration units. Outputs provide evidence, not automatic proof of novelty/collapse/expansion. Repeat annotation belongs in the separate annotation workflow.

Self-mapping evaluates consistency with contributing reads, not independent graph accuracy. Select a published workflow appropriate to supported input types; record read subsets, seeds, mapper/index versions and metrics. Avoid creating a new bespoke scoring system or automatic biological rejection rule.

## 9. Report and publication

**Primary deliverable: one detailed GitHub-ready Markdown report.**

- Analysis sections are collapsible and collapsed by default using details/summary markup.
- Organize species → analysis, with nested detail as useful.
- Keep a compact visible run-status summary so failed/skipped samples are apparent.
- Include required chromosome/haplotype figures and tables.
- Use repository-relative links to figures and downloadable artifacts.
- No JavaScript expand-all controls and no standalone HTML report requirement.
- Published-tool interactive reports may remain supporting artifacts if useful, but essential results must be readable from Markdown.

GitHub supports Markdown content within collapsed sections. [GitHub documentation](https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/organizing-information-with-collapsed-sections)

**Publish directly into the user-specified results directory.** The user manually changes it between runs. Do not create automatic run directories, latest links, directory reconciliation or deletion.

Retain a manifest of outputs from the current run. Generate the report from that manifest, not an indiscriminate scan of old files. Identify graph flavor, coordinate space, variant tier, cohort, units and status without repeating lengthy caveats in every section. No parallel FULL report.

## 10. Correctness work that remains mandatory

The [35-finding review][review] supplies source locations and detailed rationale. Recheck current source before patching.

| Review findings | Required disposition |
|---|---|
| 1–3, 15, 34–35 | Canonical taxid/assembly keys; remove Cartesian species pairing and per-species first-item broadcasts; fix optional outputs and multi-species aggregation; correct DSL2 comments |
| 4–8, 22, 26 | Repair private evidence before the requested run; reassess permanent support afterward |
| 9–10 | Replace bespoke growth logic; explicit biological units and descriptive scope |
| 11–14, 21 | Assembly representation, input parsing/routing, collision-safe IDs and harmonization bypass |
| 16–18 | Exact cut coordinates, retained short-read bypass and valid name maps |
| 19 | BlobTools/decontamination evidence unsupported; future rebuild |
| 20, 24 | Graph inclusion accounting and consistent variant/coordinate contracts |
| 23 | Rearrangement correctness included in structural-variant experiments |
| 25 | Replace similarity branch where possible; validate matrix/ordination/tree handling |
| 27 | Nonblocking recovery and conservative organelle separation |
| 28 | Preserve finishing defaults and QC modes; fix coordinate/routing defects; no new biological QC admission gate |
| 29 | Honor effective sample controls; TellSeq/ONT deferred |
| 30–32 | Enforce resource budgets; sensible retries/recovery; pin tools/databases |
| 33 | User-managed publication directory; unique artifact identities and manifest-driven report |

Pin software environments/container identities and database releases/checksums. Share immutable databases across compatible runs; do not update them in place. Validate upgrades on representative small inputs.

Budget concurrent threads, nested workflow memory and disk explicitly. Avoid globally retrying deterministic failures. Preserve expensive recoverable state until acceptance; do not alter cache records to force reuse.

## 11. Implementation order and rerun strategy

Use ordinary Nextflow resume. A future analysis-only entry point is desirable but is not part of this rewrite.

| Stage | Scope | Validation / likely rerun boundary |
|---|---|---|
| 1 | Freeze baseline inputs, configs, versions, outputs and cache locations | Static inventory; preserve prior results |
| 2 | Repair input identities, routing, skips/failures, optional channels and species joins | Tiny mixed-input/two-species fixtures |
| 3 | Repair harmonization bypass, voting eligibility, cut coordinates and finalization | Small coordinate-changing examples; no production cutting until valid |
| 4 | Implement agreed graph admission and CLIP/gref/FULL roles | Existing valid artifacts where possible; small graph fixtures |
| 5 | Replace sharing/growth/similarity components with published-tool outputs | Verify required chromosome × haplotype results before retiring custom code |
| 6 | Conduct structural-variant tool experiments | Bounded known-answer and representative cases; choose implementation |
| 7 | Repair temporary private evidence and add optional self-mapping QC | Small evidence fixtures and reproducible read subsets |
| 8 | Build Markdown report and direct-publication manifest | Singleton, mixed species, skipped/failed sample and missing-optional cases |
| 9 | Production validation | Resume valid upstream work; rerun only affected stages where scientifically valid |

Changes to assembly sequence, cohort, graph reference or construction model can require long rebuilds. Allow them. Downstream/report changes should not gratuitously invalidate construction.

Process scripts, names, environments, inputs and calling workflow names can affect cache reuse. Avoid unnecessary relocation/renaming of expensive processes early; do not preserve incorrect code merely to keep a cache. [Nextflow cache guidance](https://docs.seqera.io/nextflow/cache-and-resume)

### Acceptance criteria

- Every input row has a recorded outcome and effective processing route.
- No cross-species file pairing or unannounced cohort reduction.
- All agreed singleton/mixed-input cases produce meaningful reports.
- Failed required sample processing withholds its species' graph while independent work can finish.
- No primary/alternate pair is mislabeled as equivalent haplotypes.
- Optional analyses do not leave references to uncalled processes.
- Graph sharing uses biological path identities; synthetic gref paths never inflate denominators.
- Required chromosome/haplotype tables reconcile to their documented units.
- Variant FASTAs, genotypes, coordinate maps and representation tiers agree.
- Chimera cuts use current coordinates and preserve nonparticipating assemblies.
- Reports show current-run artifacts and explicit unavailable/failed states.
- No new manual review gate, annotation workflow, population workflow or automatic output-directory management is introduced.

## 12. Historical baseline and remaining questions

Treat previous-run numbers as reported observations until checked, not validated acceptance targets. In particular:

- The reported fall in AC=AN records after gref does not alone prove removal of false variants. Compare biological sample sets, ploidy, missingness, representation and corresponding events.
- INS/DEL asymmetry is not required to approach one.
- A chromosome returning to the cohort-median size is useful evidence but does not alone validate its repaired structure.
- Reconcile the original plan's differing individual/sample/haplotype counts using a canonical ledger.
- The configured gref image tag and current/tagged documentation differed in feature description; pin and verify the actual working image and bundled versions.
- Dataset-specific size/coverage outliers remain case investigations, not hard-coded general rules.

The outstanding design milestone is **structural-variant method selection**. Other open details are implementation contracts: Panacus export coverage, chromosome accounting, self-mapping methods and the repaired temporary evidence workflow. Resolve these with bounded tests rather than adding parallel permanent analyses by default.

## 13. Decision register

| # | Agreed decision |
|---|---|
| 1 | Automatic per-species pangenomes, globally disableable |
| 2 | Admit QC-suitable HiFi-only partially phased pairs |
| 3 | Assembly-based chromosome-scale voting suitability |
| 4 | At least two eligible individuals per species |
| 5 | Exclude collapsed diploid assemblies; genuine haploids remain eligible |
| 6 | Skip graph without suitable chromosome-scale reference |
| 7 | GREF(CLIP) sole default variant catalog; conventional export opt-in |
| 8 | Fine decomposition on demand or required by enabled consumer |
| 9 | INV/DUP/translocation required; targeted method experiments during rewrite |
| 10 | Panacus replacement preferred; chromosome/haplotype summaries mandatory |
| 11 | Haplotype sharing default; individual grouping on request |
| 12 | Expanded tiers with configurable defaults and exact boundaries |
| 13 | Temporary private evidence for requested run; reassess afterward |
| 14 | CLIP primary; FULL compact QC support |
| 15 | Haplotype heatmap, ordination and tree, preferably one Panacus matrix |
| 16 | Individual population PCA outside scope |
| 17 | Bounded self-mapping QC opt-in; new-sample mapping outside scope |
| 18 | Remove progressive-growth and whole-genome 2D branches |
| 19 | Normal resume now; analysis-only entry point future |
| 20 | Detailed GitHub Markdown, collapsed sections |
| 21 | Process eligible rows; announce skips, including Hi-C-only |
| 22 | Automatic effective assembly route with override |
| 23 | Organelle recovery not required for main assembly completion |
| 24 | Conservative separation; ambiguous matches retained; review optional |
| 25 | Keep purging defaults; repair per-sample controls |
| 26 | Keep none/final_only/all_stages QC; final_only default |
| 27 | Conditional harmonization with explicit bypass |
| 28 | Chimera detection default; cuts opt-in |
| 29 | Preserve correction/finishing defaults |
| 30 | Biological-quality warnings do not add automatic exclusions |
| 31 | Continue independent work; withhold affected species' graph after required failure |
| 32 | Pin software and databases; deliberate validated upgrades |
| 33 | Publish directly to user-specified directory; no automatic management |

## 14. Input update: multiple independent Hi-C libraries

The current dataset has one sample with one HiFi file and two independently prepared Hi-C libraries, using the same protocol/kit and originating from the same individual as the HiFi reads. Both libraries contribute to one assembly and one set of haplotype outputs; they do not create additional biological samples or pangenome members.

Retain the ordinary one-row-per-sample sheet. Add an optional Hi-C read-set table with `sample_id`, `library_id`, `readset_id`, `hic_r1`, and `hic_r2`. For this dataset it contains only the two library rows for the exceptional sample. For samples represented in that table, require the main sheet's Hi-C fields to be empty rather than silently combining the two input mechanisms. Other samples retain their existing single-pair fields.

- Validate sample references, complete mate pairs, unique read-set identities and duplicate file assignments before scheduling assembly. Repeated sample IDs in the read-set table are intentional; repeated biological-sample rows remain invalid. Resolve all read sets before deciding whether a sample has Hi-C data.
- Retain separate library identities through preprocessing and QC. The two current libraries can use the same protocol settings, but their quality need not be identical.
- Pass both paired file lists, in matching order, to hifiasm for joint phasing. Its documented `--h1` and `--h2` options accept multiple files; avoid an unnecessary large concatenated copy for this consumer.
- Wherever the selected downstream workflow removes PCR duplicates, retain library-aware handling: deduplicate independent libraries separately before combining usable contacts. Runs/lanes of the same library instead need duplicate handling across those runs. Do not introduce an extra deduplication pass merely to support multiple inputs.
- Combine compatible contacts for the sample's scaffolding workflow using the selected tool's supported input contract. Preserve library provenance and prevent read-name collisions from causing cross-library mate matching.
- Report one assembly summary with compact library contributions, QC and explicit failures in collapsed details. Do not silently drop a failed library or label a partial-input assembly as complete; apply the agreed failure policy unless the user explicitly excludes that library.
- Keep per-library preprocessing independently resumable where practical. A change to contributing Hi-C reads must invalidate affected phasing/scaffolding results; reuse upstream work only where supported and safe.

This is a planned input-handling change, not a capability of the existing single-path parser. Validate the rewrite with a single-library case, two independent libraries, two runs from one library, and invalid or duplicated mate assignments. No pipeline execution is authorized by this planning update.

Tool reference: [hifiasm Hi-C input options](https://hifiasm.readthedocs.io/en/latest/parameter-reference.html#hi-c-integration-options); [pairtools duplicate-removal and merge operations](https://pairtools.readthedocs.io/en/latest/).

[review]: <C:/Users/jselwyn/Documents/Codex/2026-09-23/referenced-chatgpt-conversation-this-is-an/outputs/pipeline-static-review.md>
[previous]: <C:/Users/jselwyn/Documents/Codex/2026-09-23/referenced-chatgpt-conversation-this-is-an/outputs/pangenome-refocus-plan-revised.md>
