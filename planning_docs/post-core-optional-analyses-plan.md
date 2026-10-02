# Post-core optional analyses: implementation and report integration

Date: 2026-10-02. Status: agreed scope and proposed implementation; not a claim of runtime validation.

## 1. Scheduling decision and precedence

Finish the current assembly, QC, core pangenome, consolidated Markdown report and cleanup before implementing these extensions. This document moves the former step 6 **after step 7**. It also brings resurrection of the nonfunctional BlobTools/decontamination-evidence workflow into the future extension phase. It supersedes older planning documents where they require SV classification before core completion or exclude that resurrection from the rework.

The revised sequence is:

1. Production HiFi assembly with both Hi-C libraries, with QC/pangenome/optional analyses off and chimera evidence on; cuts off for review.
2. Review all assemblies, chromosome assignments, reference selection and chimera evidence; preserve a baseline and choose automatic or explicit manual breaks.
3. Prepare the integrated run, including the short-read sample with simplified settings supplied by the caller script, not production defaults.
4. Resume with chosen breaks, QC and core pangenome enabled; produce final-coordinate read mappings where required.
5. Accept assemblies and core graph outputs, documenting unresolved limitations and exact input/reference identities.
6. Complete the former step 7: one functional GitHub Markdown report, publication/status handling, core cleanup and general-purpose integration checks.
7. Implement the extensions below incrementally against that accepted baseline. Each extension contributes to the existing report as part of its own acceptance, not in a later reporting project.

The current SyRI pilot may finish and its results should be retained, but it does not block core release and does not substitute for graph-wide SV analysis. Do not expand experiments solely to complete the old scope.

## 2. Core release gate

The baseline must finish without any extension enabled. A single short-read assembly, same-species phased assemblies and mixed-species/input cohorts must have meaningful outcomes and reports; unsupported analyses are explicit skips. Representative routing tests need not be expensive full biological assemblies.

Core includes final assembly/QC status, reference and cohort accounting, compact FULL-versus-CLIP support, CLIP sharing/growth, required per-chromosome/per-haplotype/chromosome-by-haplotype partitions, descriptive similarity views and existing graph variant exports. Missing required core summaries cannot be relabeled optional to close the release.

The report uses collapsed analysis sections, a visible run-status summary and relative links. It consumes current-run records, not all files found in the results directory. Disabled, skipped, unavailable and failed are distinct. Core cleanup removes obsolete dependencies and active dead branches after replacement checks, while retaining useful comparison experiments separately.

## 3. Execution and publication framework

Production is launched by `/work/birdlab/GCL/spratelloides_delicatulus_genome/scripts/run_assembly.sbatch`, publishing to `genome_assembly`. The shared work and launch directories remain `work` and `.nf/assembly`. The user owns launcher overrides and chooses the output directory explicitly.

Experiments use `scripts/comparisons` for helpers and `comparisons/{TESTID}` for results. Baselines live under `comparisons/baselines/{BASELINEID}`. No automatic directory rotation, deletion, latest links or promotion logic. A different output directory isolates published files, not the shared execution history: serialize runs using the shared session and record the explicit Nextflow resume target where necessary.

Retain work/cache, source revision, launch arguments, sample/readset sheets, tool/database versions and trace. Standalone experimental results do not automatically become Nextflow-cached production tasks. Integrate accepted implementations into stable modules/workflows; do not permanently run production analysis from comparison scripts.

### Shared inputs

- Assembly manifest: species, individual, assembly, representation/haplotype, processing stage, FASTA checksum and eligibility reasons.
- Reference manifest: selected linear assembly, checksum, selection evidence and override provenance, independent of graph execution. Singleton selection must be added before PSMC needs it without relaxing graph admission rules.
- Graph manifest: construction method/version, cohort, path identities, graph role, coordinate conventions and input assembly identities.
- Read/mapping manifest: individual/library/readset, reference identity, mapper/settings, BAM/index and completeness. Matching file names alone do not establish reuse.
- Artifact/status records: owner, analysis, inputs, outputs, completion/skip/failure reason and report links.

Each extension is independently opt-in and has narrowly scoped parameters. Any new flag names in this plan are design intent, not runnable options until implemented. Keep `main.nf` small; named workflows orchestrate tool processes. Default-off extensions must not instantiate expensive tasks or request their databases. No production sample-specific choices.

### Common report contribution contract

Every extension emits structured summary/status tables, provenance and a collapsed Markdown section with relative links to figures and downloadable data. The report collector receives explicit current-run artifacts keyed by species/individual/analysis. Required fields include analysis version, assessed scope, method, input identities, missingness/limitations and result status. Link large native files rather than embedding them.

Use tidyverse/ggplot2 for maintained R plots and patchwork for combined figures. Native published-tool viewer datasets can be supplementary; essential interpretation must remain available in Markdown. Separate expensive computation from rendering so cosmetic changes do not rerun mappings, fits or graph construction. Report inclusion and link checks are acceptance requirements for every extension.

## 4. Analysis inventory and current state

| Extension | Current state | Scope |
|---|---|---|
| PSMC | Supplied scripts reviewed; production workflow planned | Per-diploid-individual history, independent of graph enablement |
| Graph SV interpretation | Graph VCF groundwork exists; classification unresolved | Multi-sample graph alleles and justified INV/DUP/translocation interpretation |
| Optional SyRI | Branch implemented; full chromosome-set pilot running at plan date | Supplementary assembly-pair evidence, disabled by default |
| Graph self-mapping QC | Agreed opt-in; method to select | Bounded consistency checks using contributing reads |
| Private/accessory sequence evidence | Legacy branch requires repair/replacement | Temporarily retained for requested analysis; reassess long-term maintenance |
| Additional graph summaries/exports | Core figures already developed; audit actual coverage first | Only useful additions beyond required core summaries |
| BlobTools/decontamination resurrection | Nonfunctional optional evidence branch; rebuild planned | Published-tool evidence workflow, followed by separately tested cleaning integration |

## 5. PSMC

### Inputs and eligibility

Use a linear final assembly reference shared with the species graph when that reference exists, never a synthetic GREF path. With one diploid individual, select between its phased assemblies using documented available quality evidence and deterministic tie-breaking; allow an explicit override. A single FASTA from a diploid individual is distinct from a biologically haploid individual, which is unsupported for standard PSMC.

Two phased assemblies from the same individual use the Dipcall route and its confident-region BED. One assembly plus the same individual's complete HiFi reads uses diploid genotype calling. Do not use haplotype-assigned reads alone to infer diploid heterozygosity. Short-read PSMC is not part of the initial addition. Record partial-phasing limitations.

The read-based caller/model and thresholds remain a method decision, especially for nonhuman species. An existing self-mapping BAM is reusable only if reference checksum, read set and mapping requirements match the intended analysis; most other individuals need mapping to the shared reference.

### Processes and outputs

Reference preparation -> Dipcall or mapping/calling -> callable mask and diploid consensus -> PSMC input -> observed fit -> separately scheduled bootstrap replicates -> tabular export and plotting -> report contribution.

Callability must distinguish confidently homozygous sequence from missing or uncertain sequence. Record depth, mapping/genotype quality, gaps and exclusions; a variant-only VCF is insufficient. Allow explicit autosomal/interval scope and external masks; chromosome naming alone does not establish autosomal status. No annotation pipeline is added.

Retain callable BED, per-contig callable bp and heterozygosity summaries, input missingness, PSMC input, observed/bootstrap fits, parameter manifest, curves and logs. Require species-specific mutation rate and generation time for scaled plots; otherwise publish explicitly unscaled results. Expose bin size and fitting intervals. A small bootstrap pilot precedes a full requested set; verify randomization behavior of the selected build and retain generated replicate inputs/provenance.

### Acceptance and cost

Test known homozygous/heterozygous/missing intervals, masks, singleton selection, graph on/off reference consistency, true-haploid skip, poor callability, truncated fit detection and plotting scales. Compare R exports with upstream output. Validate both assembly-based and read-based input modes before claiming both supported.

Mapping/calling may be substantial; bootstraps multiply fitting work but can be scheduled independently. Scaling/style changes rerun plots only; more replicates should reuse completed fits. No assembly or graph rebuild is required.

Report: individual, reference, input mode, callable fraction, heterozygosity, assumptions, observed curve and bootstrap spread. Bootstrap spread is not uncertainty in mutation rate, generation time or genotype error.

## 6. Graph SV interpretation and supplementary SyRI

### Intended outcome

The primary target is a graph-derived multi-sample variant product with useful structural-event interpretation, not a chromosome-complete-only pairwise catalogue. Preserve original graph alleles and genotypes. Avoid treating a valid VCF as proof that every downstream tool interprets its event representation correctly.

### Implementation gates

1. Audit existing undecomposed VCFs, graph paths, GREF metadata and sample/haplotype genotypes. Document missingness, clipping and coordinate tiers.
2. Build small known-answer examples covering insertions/deletions, balanced inversions, tandem/dispersed/inverted duplication, translocations, nested variation, repeats and fragmented inputs.
3. Evaluate published extraction/classification methods against those examples and bounded real regions. Record missed classes and false positives. Do not filter SVs solely by ALT-minus-REF length: balanced events can have zero length difference.
4. Decide whether the existing graph suffices. A graph-method change is permitted and rebuilds the pangenome onward, not the assemblies.
5. Emit original calls plus interpreted event relationships and explicit QC. Preserve missing genotypes; no record is not automatically homozygous reference. Synthetic paths do not count as samples.
6. Evaluate independent read support only with a defined question and compatible mappings. A read-based cohort caller is a possible later method choice, not an already-authorized replacement for graph analysis.

SyRI remains optional supplementary evidence with explicit chromosome-set restrictions. Its skips are not absence of SVs. Preserve its current pilot, native calls and advisory QC; do not require expanding it before the core release. Existing self-BAM provenance checks do not validate SVs.

### Acceptance, outputs and cost

Preserve allele/event parentage to avoid double counting; distinguish event counts from affected bp and net sequence gain. Demonstrate balanced-event detection, appropriate fragmented-input missingness and traceable coordinates. Produce class/size/chromosome summaries only for supported interpretations; selected diagnostic plots remain optional, not a plot for every event.

Cost ranges from export/classification to a full graph rebuild; establish bounded resource measurements before expanding. New read mapping can also be expensive. Separate graph construction, extraction, interpretation, evidence and rendering.

Report: represented cohort, assessed classes/regions, classification method, original/interpreted downloads, QC, missingness and limitations. Do not label advisory-no-flags calls biologically validated.

## 7. Optional graph self-mapping QC

Inputs: CLIP graph, compatible indexes, contributing individuals' reads and explicit read types. This is consistency QC, not new-sample graph genotyping or independent truth validation.

Select a published supported mapper workflow; record versions and exact index settings. Reproducibly sample reads per individual, preserving pairs where applicable. Keep sampling/indexing/mapping/metrics/rendering separate. Linear self-mapping BAMs are not graph mappings.

Assess standard available mapping rate, identity, quality and ambiguous alignment measures with documented definitions. Avoid a custom composite score or automatic assembly rejection. Compare unlike read types only with appropriate qualifications. Generate/export indexes only for enabled consumers or explicit requests.

Acceptance: deterministic subsets, correct sample/index joins, paired-read integrity, appropriate empty/unsupported handling and bounded resource use. Report subset sizes/seeds, metrics and limitations. Sampling bounds mapping cost; index construction can remain substantial. Rebuild indexes if the graph or compatible mapper settings change.

## 8. Temporary private/accessory sequence evidence

Inputs: accepted sharing categories, biological units, traceable sequences/coordinates, source assemblies, available reads or read k-mers, and matched control regions. Optional repeat evidence comes from externally supplied annotation, not new repeat annotation here.

First specify private-to-haplotype versus private-to-individual and the observed cohort denominator. Extract intervals without conflating graph distinct sequence with assembled copy length. Prefer published alignment/k-mer tools plus minimal aggregation.

Repair or replace the legacy alignment identity calculations, search-sensitivity interpretation, private/control keys, pooled versus individual coverage, k-mer schema and units. Distinguish absent reads, low detection power and ambiguous mappings from true absence. Never automatically label unsupported sequence novel, collapsed or expanded.

Outputs: region FASTA/coordinate ledger, private/control evidence tables, support availability, compact comparisons and parameter provenance. Tests must include repeats, missing reads, known shared/private controls, multi-copy sequence and coordinate boundaries. Read k-mer/database generation may be expensive; pilot selected regions and reuse compatible stores.

Report: explicitly enabled temporary analysis, evidence and limitations. After the requested run, decide whether its value justifies permanent maintenance; do not retain duplicate bespoke and published implementations indefinitely.

## 9. Additional pangenome summaries, exports and figures

Audit actual current products before adding anything: individual/haplotype grouping may already exist. Required growth/sharing, per-chromosome/haplotype partitions, similarity heatmap/PCoA/NJ and compact clipping QC belong in the core release.

Optional additions include alternative biological grouping, focused regional views, conventional linear-reference VCF export alongside GREF, fine decomposition required by a consumer, and local graph diagnostics. Sharing segment-length distributions need a defined question; they are not SV size spectra.

Keep cohort-wide sharing definitions fixed when selecting chromosomes/haplotypes for attribution. State distinct graph bp versus assembled copy bp and denominator; keep unplaced/composite sequence explicit. Regional totals may overlap and must not be silently summed as a partition. Synthetic paths are excluded from biological counts.

Use published table exports before custom calculations. Validate on small known graphs before large node-by-group exports. Render PNG for Markdown and a suitable export format from existing tables. No unsupported openness classification, custom base-independent confidence bands, progressive-construction growth branch, individual population PCA or whole-genome 2D default.

Outputs/report: only figures answering a distinct question, their tables and a methods/denominator note. Plot-only changes reuse calculations. Expensive table exports require a bounded resource pilot.

## 10. BlobTools/decontamination workflow resurrection

### Scope and separation

Rebuild the currently nonfunctional optional BlobTools/decontamination-evidence branch, likely from scratch. This does not mean the existing production contamination-removal tools are all nonfunctional. Inventory their roles and outputs before deciding what to reuse, replace or retire.

Develop in two stages: **evidence/reporting first**, then **explicit sequence-cleaning integration**. Evidence collection on accepted final assemblies must not silently edit them. If new evidence identifies real contamination, a deliberate cleaning run creates a new assembly version and legitimately invalidates affected downstream mappings, QC, graphs and analyses. This is the exception to analysis-only reuse, not a reason to hide a required correction.

### Inputs and tool choice

Inputs include stage-identified FASTAs, target taxon, sequence lengths/GC, compatible coverage evidence, taxonomic hits and database provenance. Existing BUSCO and contamination-screening outputs can be reused only where their assembly coordinates and formats match. Missing reads or BUSCO must be visible rather than manufactured values.

Evaluate the actively supported BlobToolKit toolchain and its recommended published workflow before rebuilding custom wrappers. Check compatibility with Nextflow 23.10.1 and cluster containers/resources: do not introduce a silent runtime upgrade or nested pipeline execution. Decide between direct supported imports and a bounded integration of published components based on required features and maintenance cost.

Use existing FCS screening where appropriate; do not substitute a visualization for a contaminant classifier. Compare taxonomic, coverage and sequence evidence without assuming every non-target hit is contamination. Preserve ambiguous/cobiont/organelle interpretations and unclassified sequence.

### Processes and databases

Database preparation -> assembly metrics -> coverage preparation/import -> taxonomic searches/screen import -> BlobToolKit dataset -> summaries/static report plots -> optional cleaning decisions and separate application.

Retain shared versioned databases using the established storeDir strategy with release/checksum/completion records. Assess capacity before large downloads; avoid per-sample copies, in-place database upgrades, incomplete downloads appearing valid and concurrent writers to one store. Taxonomic search is likely among the more expensive extension tasks.

The cleaning stage must distinguish whole-sequence removal from interval trimming/splitting, emit exact decisions and coordinate maps where relevant, retain removed sequence and reconcile bp totals. No assembly- or sample-specific production exclusions. Parameters/explicit decision files belong to the run. Review is optional where a tested automatic policy is selected; evidence generation is not a mandatory manual approval gate for producing assemblies.

### Acceptance and report integration

Test clean inputs, known contamination, organelles/cobionts, ambiguous taxonomy, unclassified sequence, low coverage, missing read evidence, empty screening results and mixed species. Include a small seeded contamination case with known expected behavior, plus a bounded real assembly. Verify any removal policy separately from viewer generation and demonstrate appropriate downstream invalidation after sequence changes.

Report: assessed sequence stage/checksum, evidence completeness, GC/coverage/taxonomic overview, unclassified fraction, candidate removals or actual actions, retained/removed bp and tool/database versions. Native interactive BlobToolKit output is supplementary; static essentials and tables live in the main Markdown report. Maintained R adapters follow ggplot/tidyverse/patchwork; do not reimplement the toolkit's science simply for plot style.

## 11. Proposed order within the extension phase

All begin only after the core report/cleanup gate. Order below is a recommendation, not a reason to delay core release:

1. Finalize shared manifests/report contribution interfaces without gratuitously moving accepted expensive processes.
2. PSMC: one individual and both supported input modes, then requested cohort/bootstraps.
3. BlobToolKit evidence-only pilot and report integration; decide cleaning integration from results. Genuine contamination findings take precedence over cache preservation.
4. Bounded graph self-mapping QC.
5. Requested temporary private/accessory evidence, optionally reusing compatible evidence from prior extensions.
6. Focused graph SV method evaluation and integration; retain the SyRI pilot as supporting evidence.
7. Additional justified graph exports/figures and final extension cleanup.

No exact runtime promises before resource pilots. PSMC mapping/calling and bootstraps, taxonomic searches, graph indexes and any graph rebuild can be long; summary/plot/report tasks should generally be much smaller. Prioritize independent lightweight implementation while the user runs bounded expensive jobs, without concurrent resumes into the shared session.

## 12. Per-extension delivery and acceptance checklist

- Explicit scope, supported inputs, method decisions and parameter ownership documented.
- Published tool/version selected with a bounded scientific check; no silent fallback to a different method.
- Disabled/skipped/failed behavior, singleton and mixed-species cases tested as applicable.
- Computational inputs are stable and narrow; report timestamps and unrelated configuration do not become expensive-task inputs.
- Tests and biological jobs run by the user on the cluster; static review locally is not runtime verification.
- Handoff includes exact launch command, expected reused/new tasks, source revision and an explicit small-output collection command under scripts/comparisons. Collected archives omit large BAM/FASTA/graph files unless necessary.
- User-selected test output directory contains reports, metrics, tool versions, logs and manifests sufficient for review.
- Accepted extension runs through the normal pipeline and contributes a functioning section to the existing report; links resolve from both test and production roots.
- Repeat a resumed run to check the relevant cache boundaries; do not broaden expensive testing without a reason.
- Retire superseded implementation/dependencies only after acceptance. Preserve useful experiments and future limitations in the comparison/planning folders.

## 13. Deferred beyond this plan

Gene/repeat annotation remains external. General population analysis, routine new-sample graph genotyping, ONT/TellSeq support and a standalone analysis-only entry point are not added here. Read-based cohort SV calling is a future method option requiring an explicit scope decision. Deferred CTlk/chr12/reference-choice investigations belong to the new-Hi-C assembly review and must not become permanent sample-specific branches.

## References and related plans

- [PSMC tool and input/bootstrap/scaling guidance](https://github.com/lh3/psmc)
- [Dipcall](https://github.com/lh3/dipcall)
- [vg graph VCF export](https://github.com/vgteam/vg/wiki/VCF-export-with-vg-deconstruct)
- [Minigraph-Cactus documentation](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md)
- [BlobToolKit maintained repository and workflow recommendation](https://github.com/genomehubs/blobtoolkit)
- [NCBI Foreign Contamination Screening](https://github.com/ncbi/fcs)
- [PSMC detailed earlier design](completed-run-review-and-psmc-plan.md)
- [SV implementation checkpoint](assembly-sv-integration-checkpoint.md)
- [Core pangenome figure inventory](pangenome-figures-and-regional-plan.md)
- [Earlier agreed decisions](assembly-pangenome-agreed-plan.md)

Tool capabilities/compatibility must be rechecked at implementation against pinned releases; current upstream documentation is not proof that the installed version exposes every interface.
