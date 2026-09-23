# Assembly-only checkpoint: user-run validation

Status: source changes implemented and statically reviewed. No Nextflow command, compiler, stub, biological tool or analysis script has been run by the assistant. This is a first runtime-validation candidate, not a completed or validated pipeline rewrite. Target runtime: Nextflow 23.10.1 with the supplied Slurm profile.

This handoff supersedes the earlier batch-0 note that the companion Hi-C table was not implemented. The complete agreed design and remaining roadmap still apply.

## What to sync

Sync all changed tracked files and the new workflows, functions, modules, assets and checkpoint launcher. Copying only `main.nf` is insufficient. Changes are in the repository working tree; the assistant has not committed, pushed or submitted anything. The user's existing planning files and updated Ex2_run1/Ex3_run1 identifiers were preserved.

On the execution computer, retain the existing `.nf/assembly` launch directory and `work` directory. Existing cached jobs can only resume when their task identities and inputs still match. The new workflow nesting and read-set metadata will invalidate substantial previous work; do not expect a cost-free first baseline. Later downstream edits should preserve these new upstream boundaries.

## Prepare the two tables

Use the seven-row [main sheet](assembly_samplesheet.proposed.csv) and the two-row [Hi-C table](hic_readsets.proposed.csv). On the execution computer, place them alongside the existing data sheet, for example:

- `data/assembly_samplesheet.csv`
- `data/hic_readsets.csv`

Do not retain the duplicated Sde-CTlk_104 row in the main sheet. Its inline Hi-C cells must be empty, and both libraries must be in the companion table. Ordinary samples keep their inline Hi-C pair. The supplied relative paths are unchanged and resolve against each table's directory: with these locations, `raw_bam/...` means `data/raw_bam/...`.

Library IDs are Ex2 and Ex3; read-set IDs are the user's Ex2_run1 and Ex3_run1. No raw-read files were inspected or verified on the execution computer.

## First: input-only check

The new `run_assembly_checkpoint.sbatch` is based on the user's supplied launcher. From the same project root used for previous submissions, create the Slurm log directory before submission and run:

```bash
mkdir -p logs
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch data/assembly_samplesheet.csv data/hic_readsets.csv validate
```

This selects the `VALIDATE_INPUTS` entry point. It parses both tables, checks file existence and identities, reports row decisions, and runs a small report-writing task. It does not convert reads, download databases, assemble genomes or build graphs. The report task uses Slurm, so the handler and its small child task may both appear in the queue.

Expected outcome, assuming paths are valid:

| Sample | Decision | Expected final assembly count |
|---|---|---:|
| Sde-CBau_104 | HiFi + one Hi-C read set | 2 |
| Sde-CPla_115 | HiFi + one Hi-C read set | 2 |
| Sde-CLim_110 | HiFi + one Hi-C read set | 2 |
| Sde-CTlk_101 | Announced skip: no supported assembly reads | 0 |
| Sde-CMat_203 | HiFi + one Hi-C read set | 2 |
| Sde-CTlk_104 | HiFi + two independent Hi-C libraries/read sets | 2 |
| Sde-CMat_061 | Short-read assembly | 1 |

Inspect `genome_assembly/pipeline/input_validation.md` and `.tsv`, and the current job log. Expect six accepted samples, one deliberate skip, no invalid rows, and 11 expected assembly outputs. An input-check exit code alone is insufficient: row-level invalid inputs are reported while valid samples can still proceed. Stop and correct unexpected skips/invalid rows or counts.

Validation uses `.nf/assembly-input-check` so it does not replace the assembly launch directory's default resume history. This separates execution bookkeeping only; it does not create or manage alternate results directories.

## Then: first assembly runtime validation

After the input check agrees with the table above, the assembly-mode command is:

```bash
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch data/assembly_samplesheet.csv data/hic_readsets.csv assembly
```

This is the first execution test of the changed assembly path. It may expose compilation, tool/environment or biological-data problems that static review cannot detect. Preserve failed task directories and return the first relevant error rather than repeatedly launching the entire run unchanged.

The launcher preserves the supplied settings except for these deliberate changes:

- Pangenome is disabled (`--run_pangenome false`).
- Automatic chimera cutting is disabled (`--chimera_break false`). Its use of coordinates after correction/gap filling/telomere extension still needs the planned validation. Detection remains as configured, but cut coordinates and coordinate-dependent evidence are provisional until that work is complete.
- The optional Hi-C table is supplied with `--hic_readsets`.
- Strict shell error handling preserves a failing Nextflow exit status. Logs include mode and Slurm job ID.

The user's QC mode remains `none`, post-assembly plotting remains disabled, and the existing correction/scaffolding/finishing defaults and command overrides remain in effect. This runs the assembly path through finalization and harmonization where applicable. It does not claim biological QC passed.

Results still publish directly to the manually chosen `outdir` in the launcher. Change that line yourself between runs if desired. No automatic results-directory reconciliation or cleanup was introduced. Publishing mode remains the supplied hard-link mode; it requires the work and results locations to support hard links.

## What changed

### Stable workflow boundaries

Input preparation, read preparation, contig refinement, and Hi-C scaffolding/finishing have named workflows with explicit outputs. Existing assembly, harmonization, chimera, QC and reporting workflows remain. `main.nf` is approximately 41,000 source characters, down from its previous monolithic arrangement; this is a source measurement, not proof of compilation.

### Inputs and multiple libraries

Replaced the naive CSV parser with Nextflow's quoted CSV splitting and whole-table identity validation. The old parser was removed rather than retained as a second implementation. Duplicate sample IDs, duplicate file assignments, malformed tables, invalid identifiers and unknown sample references are explicit errors. Individual missing/invalid reads are reported with row status; the assembly run can use valid rows. Hi-C-only rows are explicitly skipped.

Each Hi-C read set is trimmed independently. A sample is released when all of its expected read sets finish, with mates ordered by read-set ID. Hifiasm receives paired comma-separated file lists, as its [documented Hi-C interface supports](https://hifiasm.readthedocs.io/en/latest/parameter-reference.html#hi-c-integration-options).

Mapping retains library read-group tags and prefixes alignment query names by read set to avoid cross-run name collisions. Read-set manifests and mapping statistics are published. Same-library runs retain the same RG ID during merging via samtools' supported `-c` behavior. [Samtools merge documentation](https://www.htslib.org/doc/1.10/samtools-merge.html)

Pairtools reads the RG tags into extra columns and only marks coordinate duplicates when both library-tag columns also match. This implements within-library duplicate handling in a combined stream rather than physically splitting the BAM into libraries. It avoids a bespoke duplicate-removal algorithm. The filter environment now pins pairtools 1.1.2, a published release; the new environment needs to resolve on the execution computer. [Pairtools options](https://pairtools.readthedocs.io/en/latest/cli_tools.html), [published package](https://pypi.org/project/pairtools/)

### Assembly-path fixes

- Downstream routing and QC read selection use the effective assembler rather than merely the presence of shotgun reads.
- Per-sample purge_dups selection is honored. `dedup=none` disables Redundans reduction while keeping independently configured scaffolding/gap closing.
- Unsupported assembler/dedup combinations and primary/alternate-to-haplotype relabeling are rejected. Assembly file counts are checked before haplotype assignment.
- Hifiasm's per-task option variables are local to avoid shared-variable races.
- YaHS now receives its configured `min_mapq` key; the old code read a different key and fell back to 1.
- Finalization preserves the staged `input/` path. Actual sentinel assets are supplied for relevant optional harmonization/chimera inputs.
- A globally empty harmonization candidate stream no longer errors for a singleton/short-read-only run. A candidate table emitted for an eligible group is still required to contain candidates.
- Chimera evidence attaches candidate/reference information by species, and the cutting branch preserves the short-read bypass. Coordinate correctness is not claimed fixed; the checkpoint disables cutting.
- A successful MitoHiFi command with no complete final mitogenome now emits a no-result status instead of an artificial nuclear-assembly failure. Nonzero tool exits remain errors. Missing GenBank annotation is no longer represented by an empty successful annotation file. Optional no-result completion participates in bait handling; this join can wait for the HiFi recovery branch to finish.
- BUSCO database and Meryl tasks are not launched when their consumers are disabled. The nonfunctional BlobTools evidence option is rejected with an explicit explanation instead of entering its broken branch.

### Reporting

Input decisions are published before assembly. A compact `assembly_run_summary.md` is generated from current-run finalized outputs even with QC and post-assembly plotting disabled. It contains collapsed sample/output sections and relative links. This is the foundation for the final consolidated report, not the completed reporting redesign. On an aborting task failure, the final summary may not be generated; use the current job log and early input report, not an old report in the shared results directory.

## What was checked, and what remains

Static checks covered moved call sites, input/output shapes, active include paths, module output declarations, configuration selectors, conditional outputs, launcher text and Git whitespace errors. The one missing path reported by the broad include-text scan was an already-commented-out legacy include, not an active dependency. No runtime verification has occurred.

The full rewrite is not finished. Outstanding work includes complete per-species required-failure handling, the broader tool/database pinning audit, organelle reference fallback and failure isolation, full chromosome-scale voting eligibility, coordinate-safe chimera operations, comprehensive QC-mode integration and the final consolidated report. The existing task retry/termination policy remains; independent samples are not yet guaranteed to finish after an exhausted required-task failure. The provisional input validator also blocks all graph construction when any assembly input is invalid, rather than implementing the final per-species policy. Pangenome construction and analyses remain outside this checkpoint.

After launch, return the input report and, if a task fails, its process name, exit status, current Nextflow log excerpt and task `.command.err`/`.command.log`. For successful multi-library processing, check that both Ex2_run1 and Ex3_run1 appear in the mapping read-set manifest and both library tags are represented. Those checks establish the next debugging step; they do not establish phasing/scaffolding accuracy by themselves.
