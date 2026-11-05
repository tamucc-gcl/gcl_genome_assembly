# Frozen core interface — 1.0

Frozen 2026-11-05 against Nextflow 23.10.1 and the accepted detailed-report-01 run.
The interface covers supported inputs, flags, output locations and coordinate/decision
semantics. It does not freeze biology, dependency image tags or optional future analyses.
Changing this contract requires documentation and relevant regression checks.

## Entry points and configuration

The production entry point is `main.nf`; `-entry VALIDATE_INPUTS` validates inputs
without running assembly. Use normal Nextflow `-resume` with the original launch
session and work directory. Nothing reconstructs deleted cache files.

```bash
# Run inside the repository with project-specific paths/settings.
nextflow run main.nf -profile slurm -resume \
  --sample_sheet /project/data/samples.csv --hic_readsets /project/data/hic.csv \
  -params-file /project/data/assembly-params.json --outdir /project/results
```

The Hi-C table is optional. Relative read paths resolve against the table containing
them, not the launch directory. Explicit CLI settings override the parameter file.
Default values live in [nextflow.config](../nextflow.config); execution profiles
contain site resource/storage settings and must be adapted for another cluster.
The [default inventory](parameter-defaults.tsv) records the core snapshot. It does
not guarantee that every advanced parameter is appropriate for another organism.

For the Crest sibling layout `project/gcl_genome_assembly`, the checkpoint launcher is:

```bash
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch \
  data/samples.csv data/hic.csv assembly -params-file "$PWD/data/assembly-params.json"
```

Modes: `validate`, `assembly`, `graph_validate`, `graph_build`. Use `NONE` as the
second positional argument when there is no Hi-C table. Assembly mode sets graph
construction off unless explicitly overridden; graph modes select validation/build.
All biological tuning comes from the sample table, parameter file or explicit flags.
The convenience `run_assembly.sbatch` delegates to assembly mode with the same first
two arguments and forwards subsequent Nextflow options. Launchers use symlink publication.
Direct Nextflow's configurable `publish_dir_mode=link` means hard links, not symlinks.

## Inputs and identity

The sample CSV requires `sample_id`. Supported columns include `hifi_bam`,
`hic_r1`, `hic_r2`, `sr_r1`, `sr_r2`, `taxid`, `species`, and `ploidy`.
Omit unused columns or leave them empty. Sample IDs must be unique, beginning with
a letter/number and containing only letters, numbers, underscores, periods or hyphens.
Mate files must occur together; the same physical read file cannot be assigned twice.

```csv
sample_id,taxid,ploidy,hifi_bam,hic_r1,hic_r2,sr_r1,sr_r2
individual_A,12345,2,reads/A.bam,reads/A_R1.fastq.gz,reads/A_R2.fastq.gz,,
individual_B,12345,2,,,,reads/B_R1.fastq.gz,reads/B_R2.fastq.gz
```

An additional Hi-C CSV uses exactly the identity columns
`sample_id,library_id,readset_id,hic_r1,hic_r2`. Read-set identity is unique within
a sample; multiple read sets may belong to one library. Independent library
identities survive scaffolding/evidence generation. Shotgun reads are not Hi-C.

Row taxid/name overrides the corresponding global context according to input
validation. Positive taxids are authoritative; scientific-name-only input resolves
through local NCBI taxonomy with ambiguity checks. Species remain isolated throughout.
Ploidy describes the biological individual; it is not inferred from other samples.
Hi-C-only/no-assembly-read rows produce an explicit skip, not an empty assembly.

HiFi with Hi-C uses phased hifiasm and YaHS; HiFi alone uses hifiasm without Hi-C
scaffolding. Paired shotgun input uses the supported SPAdes/conditioning route and
its organelle workflow. Linked-read columns are validated but no linked-read
scaffolder is supported. Shotgun + Hi-C scaffolding is not wired.

## Primary controls

| Control | Supported behavior |
| --- | --- |
| `qc_mode` | `none`, `final_only`, `all_stages` |
| `run_report` | Main Markdown reporting independent of optional post-assembly work |
| `report_chromosome_painting` | Coordinate-matched painting from the selected final linear reference |
| `run_final_contact_maps` | Final contact plots; QC still controls final Hi-C evidence/metrics |
| `run_final_hifi_mapping` | Final-coordinate read evidence; distinct from pre-finishing chimera mapping |
| `chimera_detect`, `chimera_plots`, `chimera_hifi_context` | Candidate/evidence production and presentation |
| `chimera_break` | `false` for evidence only, or an absolute reviewed TSV path; `auto` is rejected |
| `run_post_assembly`, `run_pairwise_alignments` | Optional finalized comparison branch |
| `pairwise_alignment_mode` | Within-sample, reference or all-pairs modes as validated by the comparison planner |
| `run_pangenome`, `pangenome_validate_only` | Optional eligible-species graph build or input-only validation |
| `pangenome_growth`, `pangenome_popstruct` | Accepted sharing/growth/regional and descriptive similarity/ordination/tree outputs |

Reference selection is independent of graph admission/construction. Fragmented phased
assemblies may contribute graph sequence without becoming chromosome voters or references.
Completeness and scope are explicit; missing chromosome assignments do not imply absence.
Species-specific parameters (enzymes, mutation rates, motifs, assembly tuning) are not
project-specific defaults embedded in general launch scripts.

The SyRI/reference-comparison branch remains a separately enabled pilot, not accepted
graph-wide SV interpretation. PSMC, BlobToolKit rebuilding, graph read-mapping QC,
private/accessory validation, topology-based SV interpretation and extra graph products
remain [planned extensions](../planning_docs/post-core-optional-analyses-plan.md).
Retired prototypes for these products are not working optional modules. Requests to
enable explicitly retired pangenome branches fail with an unsupported-feature message.

## Manual cuts and coordinates

The first pass writes `assembly/chimeras/README.md`, `review.tsv`,
`assembly_registry.tsv`, `cut-instructions.md`, `coordinate_status.md` and
`assemblies/<assembly>/report.md`, with measured evidence and plots.
Keep the combined evidence overview, confidence/competing-placement figure and
per-library junction maps together. Whole-scaffold contact maps pool libraries.

Copy the review TSV outside generated results before editing. Each event has an
unselected event-level row; gap and alternative-location rows are separate options.

| Decision | `review_disposition` | `selected` | Required changes |
| --- | --- | --- | --- |
| Keep | RETAIN | NO | Reviewer and reason |
| Defer | DEFER | NO | Reviewer and reason |
| Cut a supplied gap | CUT | YES | Keep gap coordinates/action; reviewer and reason |
| Internal/manual cut | CUT | YES | Exact `cut_bp`, `BREAK_PROBABLE_MISJOIN`, localized status; reviewer and reason |

For an undetected cut, add a unique row with the assembly/scaffold identity,
`assessment_sha256`, `coordinate_stage=pre_finishing`, exact coordinate, reviewer
and reason. Preserve all current schema columns; the emitted instructions define
the required fields. Do not reuse a cut sheet on a different assessment FASTA.

All cut coordinates are zero-based on original pre-finishing scaffolds:
`cut_bp` separates `[0,cut_bp)` from `[cut_bp,length)`. Gap cuts preserve gap Ns.
Multiple cuts use original coordinates. The cutter validates source fingerprints,
coordinates, literal gaps, conflicts and minimum pieces and verifies exact parent
reconstruction. It emits cut audits, coordinate lifts and verification JSON.
Chromosome reassignment, finishing, final naming, final mappings and graph inputs
consume corrected sequences. A pre-finishing plot is not a final-coordinate mapping.

An assay site, low-support region, plot marker or automated recommendation is not
an applied cut. Only submitted actions with verification are reported as applied.

## Stable published outputs

All paths are relative to `outdir`. Enabled analyses contribute their own outputs;
disabled sections are omitted rather than represented by empty success figures.

| Path | Contract |
| --- | --- |
| `assembly_report.md` | Concise GitHub Markdown main report and current-run links |
| `reports/main/` | Detailed evidence index, decision/scope tables, completeness and report manifest |
| `assembly_run_summary.md` | Declared/accepted/skipped/completed assembly inventory |
| `assembly/final/<id>.fasta`, `<id>.name_map.tsv` | Final sequence and corresponding naming map |
| `assembly_eligibility.json`, `.tsv`, `.md` | Cohort membership, scope and shared reference decisions |
| `references/<taxid>/reference.json` | Final linear reference identity, fingerprint and sequence dictionary |
| `assembly/chimeras/` | Original-coordinate evidence, review instructions, selected-action verification |
| `qc/assembly/assembly_qc_metrics.csv` | Stage-aware QC measurements; missing values remain missing |
| `reports/assembly_qc_report.html` | Retained interactive ranking report |
| `contact_maps/`, `qc/hic/` | Final Hi-C maps, retention and cis/trans metrics |
| `organelle/` | Recovery status, assemblies and available annotations/plots |
| `pangenome/<taxid>/` | Graphs, construction/representation audits and enabled summaries |
| `provenance/` | Tool observations, dependency declarations, resolved execution context and issues |
| `report_validation.json`, `.md` | Post-publication report/link checks; does not certify biology |

Name-map chromosome assignment percentages describe named scaffold lengths; painting
coverage is a separate measurement. Composite names do not themselves establish a
true fusion. Graph representation flags are not automatic assembly-edit requests.

## Failures, caches and version records

Requested task failures retry according to resource policy and then stop. Unsupported
input routes and inadequate scope must have explicit status/reason. Missing measurements
must not be fabricated as zero. Do not mark a failed analysis complete because a report exists.

Large mappings/graphs remain in Nextflow work storage. Published symlinks depend on that
storage. Do not delete the cache and expect `-resume` to reconstruct missing inputs.
Database store directories are snapshots; refresh by selecting a new snapshot directory.

Version records use one nonempty `tool<TAB>version` line, optionally a three-column
`process<TAB>tool<TAB>version` schema. Multiline tool output is normalized into the
single version field. Unknown/malformed observations remain explicit issues, never
guessed from bare numbers. Existing cached records are immutable; producer formatting
fixes take effect when those producers rerun naturally. This cleanup changes only
Cooler/tidk version emissions, not their scientific commands or evidence thresholds.
