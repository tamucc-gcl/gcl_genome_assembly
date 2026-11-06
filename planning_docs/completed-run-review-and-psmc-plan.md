# Completed-run review and optional PSMC integration
Prepared 2026-09-28. Planning and static script review only; no pipeline, PSMC, R, or tests executed.

## Immediate checkpoint: assess the completed run

The user restarted with the short-read sample excluded. Expect five HiFi individuals and ten final haplotype assemblies. The Hi-C-only sample is a reported skip if retained in that run's sheet; an absent row is not a skip. Short-read execution is deferred, not validated or removed from supported inputs.

Provide the final Nextflow log, Slurm launcher output, exact submitted samplesheet and Hi-C readset table, launch command/overrides, and execution trace if produced. Also provide assembly_run_summary.md from that invocation; a listing of assembly/final with file sizes; final .fai and .name_map.tsv files; harmonization reference/voter/candidate reports; and small chimera/finishing reports. Large FASTAs/BAMs are not needed for initial review. A directory listing or archive of these small reports is sufficient.

Review gates:
1. Match run name/revision, command, and input membership. Require actual successful workflow termination; a scheduler job ending and existing result files are insufficient.
2. Reconcile latest task attempts: completed/cached outcomes, resolved timeouts, ignored errors, failures, and unexpected missing downstream processes. Check the double-library sample's two haplotypes and both mapping rounds explicitly.
3. Reconcile ten final FASTAs, indexes, and maps against those exact input individuals. Detect stale short-read products from earlier runs without deleting them. Published file existence alone does not identify current-run output.
4. Compare FAI total length, largest scaffolds, sequence count and stage summaries per haplotype. Investigate large unexplained losses/gains, truncated/empty output, incorrect map lengths, duplicate names, and divergence between homologous assemblies. N50 alone is not assembly validation.
5. Review decontamination removals, correction/scaffolding changes, gap filling and telomere extension for anomalous changes. Review Hi-C library contribution and filtering retention, not just job completion.
6. Review harmonization reference, voters, unplaced/composite assignments and chimera flags. The completed run predates the coordinate fixes; its chimera positions are not approved for automatic cutting. Existing reports cannot establish correct liftover through sequence correction.
7. QC mode none means biological quality is unassessed beyond available evidence. Missing BUSCO/consensus/contamination QC is not a pass.

Do not require the new eligibility report or coordinate-audit files from this old run: those arrive after syncing the rewrite.

After review, preserve old logs/small reports and keep work/cache/history intact; sync all edited and new files; run the documented synthetic tests and input validation remotely; then resume the assembly checkpoint. Keep pangenome, PSMC and automatic cutting off. Expected upstream cache reuse is verified in the log, not guaranteed. The filter_hic_bam label now requests 96h. A reduced testing sheet should exclude the short-read sample rather than repeatedly launching/cancelling it. Retain a separate future short-read acceptance case.

## Scope and scheduling of PSMC

Add optional POST_ASSEMBLY_ANALYSES -> PSMC_ANALYSIS after final assembly acceptance. It is an explicit, bounded exception to the earlier exclusion of individual-level population analyses; broader population genetics remains outside scope. Default run_psmc=false. This original design is superseded by the implementation in docs/psmc.md.

The PSMC flag must work independently of run_pangenome, qc_mode and optional visualization switches. Disabled means no mappings, calling, bootstrap or plotting tasks registered. Enabled but unsupported inputs yield a reasoned skip; failed requested analysis is not represented as a successful empty plot.

### Individual reference selection (corrected 2026-11-06)

PSMC uses a reference from the individual being analyzed, independently of the species graph reference and cohort readiness. Record assembly ID, species, individual, final FASTA identity/checksum, sequence dictionary and selection reason.

Use the individual's own hap1 for a finalized pair, or its sole assembly otherwise. This deterministic choice is not a quality ranking. The previous shared-reference design departed from the user's intended example and was corrected after the first pilot. Do not restore it. Different individual references do not prevent combining demographic curves. Keep species and scaling groups separate, and retain per-individual masks and reference provenance.

### Input modes

- Two phased haplotypes from the same diploid individual: retain Dipcall as the assembly-based mode. It produces variants and a confident-region BED; preserve the BED, since PASS variants alone do not define callable sequence. Never pair haplotypes from different individuals. Partially phased assembly limitations must be reported. [Dipcall documentation](https://github.com/lh3/dipcall)
- One FASTA representing a diploid individual plus that same individual's HiFi reads: map the full diploid read set to the chosen reference and call diploid genotypes. Do not restrict to haplotype-assigned reads. Caller/version/model and filtering remain a targeted method decision, especially for nonhuman species; a human-trained model is not assumed appropriate.
- A biologically haploid individual is unsupported for standard PSMC. A single haplotype FASTA from a diploid individual is a different case. PSMC models a diploid sequence. [PSMC documentation](https://github.com/lh3/psmc)
- Missing HiFi for the single-FASTA mode: skip with reason. Short-read PSMC support is not part of this addition.

Use linear final-assembly coordinates and input callability, not graph VCFs. Separate modes in provenance; compare their heterozygosity/callability before treating estimates as interchangeable.

### Callability and sequence scope

The read route must distinguish confident homozygous sites from uncovered/uncertain sites: a variant-only VCF cannot do this. Combine coverage, mapping/base quality, genotype confidence, depth bounds, and explicit exclusion masks. Mask ambiguous mapping, indel/complex regions as appropriate, gaps, organelles and identified non-diploid regions. Preserve missing sequence as unknown, never turn it into homozygous reference by default.

Replace the hard-coded chr[0-9]+_1 selector with explicit contig/interval scope in reference coordinates. Names and chromosome-scale status do not establish autosomal identity. Allow a user autosome BED/list and exclusions; unknown sex-chromosome status is reported as unresolved. Do not require repeat annotation within this pipeline; optional external masks can be consumed.

Emit per-individual and per-contig callable bp/fraction, heterozygous SNP counts and rate per callable bp, depth summaries, exclusions, PSMC input missingness and segment counts. Insufficient callability/heterozygosity should produce a clear unsupported/unreliable outcome, not an apparently meaningful empty curve.

### Processes, resources and reproducibility

Separate reference preparation, Dipcall OR read mapping/calling, callability/consensus, PSMC input, observed fit, bootstrap replicates, and plotting/report aggregation. Use Nextflow scheduling with one bootstrap per task and a concurrency limit, rather than nested GNU parallel consuming the allocation. Give expensive mappings/calls stable input contracts; rescaling plots must not rerun them. More bootstraps should reuse previous fits where supported.

Pin tools and helper scripts; eliminate runtime wget from master. Record parameters and independently controlled seeds using the chosen PSMC build's supported interface after verification. Expose fitting parameters, bootstrap count, block/bin settings and filters. Treat ten replicates as a pilot; propose 100 for a full analysis, with the user's compute budget controlling the choice.

Require positive species-specific mutation rate and generation time for years/Ne plots; do not carry the example fish values into every dataset. Otherwise provide explicitly unscaled results. The example fitting intervals are not universally appropriate: the author identifies the example as human-oriented and advises inspecting interval information. [PSMC documentation](https://github.com/lh3/psmc)

### Assessment of supplied scripts and publication

The shell script's Dipcall-to-PSMC conversion follows the documented approach. Preserve this published-tool route, while replacing hard-coded paths, sample/reference IDs, scaffold-name selection, unchecked shell execution and runtime downloads. Its QC section is currently comments, not implemented measurements.

The R scaling matches the documented formula for a 100-base bin. Pass the actual bin size explicitly. Its directory-wide file discovery can include stale replicates; consume an explicit current-run manifest instead. Require exactly one observed fit and the requested bootstrap set, valid positive scaling values, finite estimates and complete iteration blocks; it currently accepts EOF without the closing marker.

Prefer the upstream plotting tool where it meets requirements; retain only a small validated export/plot adapter if needed for fractional generation times, CSVs or combined figures. Compare identical selected iterations, bin scaling and terminal-interval handling against upstream output. The upstream plotter has its own iteration selection and last-interval policy. [PSMC plotting source](https://github.com/lh3/psmc/blob/master/utils/psmc_plot.pl)

Keep raw observed/bootstrap .psmc files, .psmcfa, callable mask, small QC tables, reference/parameter manifest and logs in addition to CSV/plots. These support reproducibility and rescaling without remapping. Large mappings follow the pipeline's work retention policy; do not automatically discard the scientific inputs as the example comment suggests.

Publish under post_assembly/psmc/<species>/<individual> with one collapsed section in the primary Markdown report: input mode, reference, callability, assumptions, completion and links to curves/data. Combined plots preserve individual and species identities. Bootstrap spread is not uncertainty in mutation rate, generation time, model assumptions or genotype error.

### Acceptance before enabling

Remote bounded checks: known homozygous/heterozygous/missing intervals; mask boundary correctness; one diploid pair; one single-FASTA diploid with HiFi; true haploid rejection; single individual without graph eligibility; graph on/off gives the same eligible-cohort reference; mixed species; empty/poor callability; missing bootstrap/truncated fit; nondefault bin size and fractional generation time; changed plotting parameters reuse fits. A small representative biological comparison must validate the new read-based mode before production use.
