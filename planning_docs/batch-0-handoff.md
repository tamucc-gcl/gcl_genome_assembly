# Batch 0: first structural change and execution baseline

Status: first implementation slice of batch 0, statically reviewed; not compiled or run. Remaining workflow extractions and the batch 1 input rewrite are not complete. No commit, push, job submission or pipeline execution has been performed.

## Execution baseline supplied by the user

- Nextflow module: `nextflow/23.10.1`. Target its DSL2 syntax; do not introduce a new runtime/parser requirement during structural extraction.
- Launcher inspected as text: `C:/Users/jselwyn/Downloads/run_assembly.sbatch`. It has not been changed or executed.
- Entry point on the execution computer: `${PROJECT_ROOT}/gcl_genome_assembly/main.nf`; profile `slurm`.
- Launch directory: `${PROJECT_ROOT}/.nf/assembly`. With the normal local cache configuration, its `.nextflow/cache` contains task-cache metadata. Preserve the launch directory and `${PROJECT_ROOT}/work`.
- `-resume` is already present. Record the intended session when distinguishing fixture runs from production runs; do not accidentally resume the last unrelated fixture run.
- Data sheet defaults to `${PROJECT_ROOT}/data/assembly_samplesheet.csv`. Current relative read paths resolve against the sheet's parent, so `raw_bam/...` implies `${PROJECT_ROOT}/data/raw_bam/...`. Preserve that convention for both proposed sheets, placed together. Remote input existence was not checked.
- Results publish directly to `${PROJECT_ROOT}/genome_assembly`; current mode is `link`. Retain user-directed output-directory management. Hard links require compatible filesystems; do not change this override silently.

Keep the user's existing command overrides distinct from pipeline defaults. In particular, this launcher requests purge_dups, disables QC and post-assembly analysis, selects all pairwise alignments, enables the pangenome, and uses `chimera_break=auto`. Structural extraction does not change these settings or their current interpretation. The implementation needs an explicit typed contract for the `auto` option rather than silently reinterpreting it as a Boolean.

The current main workflow gates REPORTING on both QC and post-assembly outputs. Consequently this launcher's `qc_mode=none` and `run_post_assembly=false` suppress it. The agreed rewrite must provide a core status report independently, with optional sections marked disabled or unavailable. Do not enable costly analyses just to obtain that report.

The Slurm log directory must exist before submission because Slurm opens `logs/assembly-%j.out` before the shell body runs. The script also returns to the project directory after Nextflow; a later launcher-hardening change should explicitly preserve Nextflow's failure exit status. These are recorded for follow-up, not silently changed in the supplied launcher.

## Implemented slice

Added `workflows/read_preparation.nf` with the named workflow `READ_PREPARATION`. Moved the existing conversion, trimming, read-bundle assembly and genome-size estimation orchestration into it. Updated `main.nf` consumers to use named workflow outputs. Tool process definitions, commands, parameters and biological routing were not edited.

This boundary owns these processes permanently unless a later correctness requirement justifies a change:

| Previous task name | New qualified task name |
|---|---|
| BAM_TO_FASTQ | READ_PREPARATION:BAM_TO_FASTQ |
| TRIM_HIC | READ_PREPARATION:TRIM_HIC |
| TRIM_SHORTREAD | READ_PREPARATION:TRIM_SHORTREAD |
| ESTIMATE_GENOME_SIZE | READ_PREPARATION:ESTIMATE_GENOME_SIZE |

Task names change even though tool scripts are unchanged. Expect cache misses for these moved tasks and potentially their descendants because output paths change. No selector targeting these process names was found in the repository configuration scan; external/site configuration remains outside this inspection.

### Current extraction contract

Inputs are the existing `(meta, reads)` sample channel and `(sample, ploidy)` channel. Outputs are:

| Output | Shape and intended cardinality |
|---|---|
| reads | `(meta, hifi, hic_r1, hic_r2, sr_r1, sr_r2)`, one per successfully prepared sample |
| hifi | `(meta, fastq)`, one per HiFi sample |
| hic | `(meta, r1, r2)`, one per Hi-C sample under the existing single-library contract |
| shortread | `(meta, r1, r2)`, one per shotgun sample; trimmed when enabled, raw otherwise |
| qc_reads | `(meta, reads)`, current assembly-QC read selection |
| genome_size | `(meta, size_file)`, one per successful estimation |
| genome_results | `(meta, summary_file)`, one per successful estimation |
| versions | Existing version artifacts from short-read trimming and genome-size estimation |

Optional short-read trimming has a pass-through output, so consumers do not directly reference an uncalled process. No additional wait for the whole cohort was added. Existing failure behavior, metadata-selection problems and version-coverage gaps are intentionally not represented as fixed by this extraction.

## Remaining batch 0 ownership work

Keep CONTIG_ASSEMBLY, ORGANELLE, QC_PHASE, HARMONIZE_SCAFFOLDS, CHIMERA, PANGENOME and REPORTING at stable named boundaries where feasible. Do not add a top-level wrapper that renames every otherwise unchanged task merely to reduce line count.

| Remaining boundary | Ownership and contract to finalize before costly validation |
|---|---|
| Input/taxonomy preparation | Unique sample rows, optional library table, taxid-keyed lookup and sample identity/traits; explicit expected cohort and skip records |
| Reference/database preparation | Shared versioned references and database readiness; species-keyed outputs, independent of report settings |
| Contig refinement | Organelle separation, purge/correction/decontamination and short-read conditioning; assembly-stage records plus retained QC checkpoints |
| Hi-C scaffolding | Mapping/filtering/scaffolding rounds with explicit assembly coordinate version and library identity; no stale BAM/AGP reuse after sequence changes |
| Finishing | Gap filling, telomere extension and bypasses; final pre-harmonization assembly plus provenance and checkpoints |
| Harmonization/chimera/finalization | Species-keyed votes and bypasses, coordinate/name maps, unchanged nonparticipants and finalized assemblies |
| Graph construction versus analysis | Stable cohort/reference manifest and graph-role outputs; downstream analyses must not force construction edits for presentation needs |
| QC/report assembly | Stage-labeled assembly inputs and independent artifact/status records; optional QC must not suppress core reporting |

Each new boundary must expose named channels rather than force unrelated consumers to reference internal process names. Separate sample/haplotype/library identity and sequence-coordinate versions. Keep scientific default changes out of extraction commits; repair known behavior in the subsequent ordered batches.

## Concrete supplied dataset

The supplied eight rows describe seven biological samples: five HiFi samples, one short-read sample and one Hi-C-only sample. Sde-CTlk_104 has two independently prepared libraries from the same individual, with the same protocol. Sde-CTlk_101 should be retained as an announced assembly skip. Expected assembly jobs: six, subject to actual input validation. Do not count the duplicate library row as a new individual or infer graph eligibility from these counts alone.

Prepared two future-format files without modifying any raw-read path strings:

- [assembly_samplesheet.proposed.csv](assembly_samplesheet.proposed.csv): seven unique sample rows; Sde-CTlk_104's Hi-C cells are empty because its two libraries are in the companion table.
- [hic_readsets.proposed.csv](hic_readsets.proposed.csv): two rows for Sde-CTlk_104. Ex2/Ex3 library labels and read-set labels were assigned from the supplied filenames; they are identifiers, not inferred protocol differences.

**These sheets are not ready for the current parser.** Multi-library parsing and consumers are batch 1/2 work. Do not switch production to the proposed main sheet alone: that would omit both Hi-C libraries for Sde-CTlk_104. Do not use the duplicate-row sheet to try to enable this feature either; current parsing can skip the Hi-C-only duplicate instead of incorporating its library.

The batch 1 contract will introduce an optional `hic_readsets` parameter. Both table paths resolve relative to their own containing directory; keeping both tables alongside the existing sheet preserves the supplied relative paths. A sample listed in the library table must have empty inline Hi-C cells. Validate unknown sample references, missing mates, duplicate file assignments and library/read-set identities before scheduling work. Unlisted ordinary samples retain the original inline-pair format.

## Static verification and next checkpoint

Reviewed the source diff, moved body, emit names and all previous process-output references in main. Inspected module output declarations and repository configuration selectors. No Nextflow compiler, stub run, assembler, analysis script or data validation was executed. These checks cannot establish runtime correctness.

Do not start an expensive biological run for this extraction alone. Finish the upstream structural and input-contract changes first. The first user-run checkpoint should then establish Nextflow 23.10.1 compilation and bounded routing cases, followed by representative data only after those pass. Stub support must itself be inspected before suggesting that every branch is safe to stub.
