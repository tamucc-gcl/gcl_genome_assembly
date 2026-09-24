# Shared databases and assembly restart — 24 September 2026

## What the log establishes

The run completed 94 tasks. All five HIFIASM tasks finished with exit 0. SPAdes for Sde-CMat_061 exited with status 143 twice, consistent with the reported cancellations. The first cancellation triggered the configured retry; the second exhausted maxRetries=1 and aborted the pipeline. Nextflow then cancelled ten remaining tasks. The log does not establish a SPAdes software defect or out-of-memory failure.

The other startup warnings were two selectors for inactive evidence modules and the intentional Sde-CTlk_101 Hi-C-only skip.

## Changes made

- DOWNLOAD_TAXDUMP, FCS_DB_GET, DOWNLOAD_GETORGANELLE_DB and DOWNLOAD_MITOS_DB now use storeDir at their existing configured shared database locations. Downloads produce task-local payloads that Nextflow moves into the store. Downstream directory outputs remain compatible.
- BUSCO already used storeDir. Its download now fails promptly on command failure or a missing/empty dataset.cfg.
- Taxonomy requires names.dmp and nodes.dmp. FCS requires the nine documented database files. GetOrganelle stores its SeedDatabase and LabelDatabase directories; MITOS requires every requested RefSeq directory. Missing MITOS extraction output now fails instead of merely warning and creating a ready sentinel.
- New GetOrganelle downloads check all seven seed/label FASTAs, both VERSION files, and invoke the tool's own database check.
- Database download stubs cannot create dummy payloads in shared stores.
- Force-download flags now produce a clear configuration error. Refresh a database by selecting a new shared snapshot directory. Existing defaults and directory locations are unchanged.
- FCS downloading now requests 2 CPUs and 8 GB on normal; cleaning requests 2 CPUs and 16 GB on normal (32 GB on retry). Neither needs to load the large screening database. These are initial resource allocations to verify from the next trace. Screening retains its existing 500 GB allocation.
- Removed DIAMOND_BLASTX and FCS_BLOB_EVIDENCE_REPORT configuration selectors. Invalid sample rows still produce warnings; expected Hi-C-only skips remain visible as INFO messages and in the input report.
- The normal Nextflow progress dashboard is retained. The earlier -ansi-log false workaround was removed. Input validation now produces a plain capability map before workflow registration; impossible assembly/read-type branches are not called.
- Startup contact-map and pairwise-synteny settings now respect run_post_assembly=false.
- Explicit LF checkout rules were added for Nextflow/config/shell/Slurm sources.

## Shared-store behavior and limits

storeDir reuses declared outputs independently of -resume and the work-directory cache. It does not compare the download URL or script against the stored content and does not automatically refresh a database.

Existing complete stores should therefore be reused without another download. Old sentinel files are harmless and are no longer consulted. An incomplete store causes its download task to run again; a fresh large FCS download needs temporary space in the work filesystem as well as the destination store.

This is existence-based reuse, not a checksum audit on every run. Existing BUSCO, GetOrganelle and MITOS directories are assumed to contain complete, usable datasets. New-download checks execute only when a download actually runs. Do not manually alter a snapshot while it is in use. Coordinate the first initialization of a shared directory: storeDir is not a cross-run locking service. Use one FCS profile/database prefix per configured directory.

The disabled bespoke DIAMOND/BlobTools branch is not rebuilt here. Species-specific MitoHiFi reference selection remains an ordinary resumable task rather than a shared database service.

## Restart

Sync the repository changes to the cluster, preserving both the existing work directory and .nf/assembly/.nextflow cache. Use the same launcher command:

~~~bash
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch data/assembly_samplesheet.csv data/hic_readsets.csv assembly
~~~

The launcher already uses -resume. No assembler/read-preparation process command or qualified workflow name was changed. Input parsing and workflow guards changed; existing sample metadata and task inputs are preserved. Completed hifiasm tasks should be reusable provided their inputs, environments, work files and cache remain available. Confirm the cache messages rather than assuming every task will match. Database/task-path changes may rerun some preparation steps. Cancelled SPAdes work is not a completed Nextflow cache entry and will rerun.

The FCS resource split removes unnecessary pressure from downloads and cleaning, but screening and SPAdes still compete for high-memory capacity. SPAdes resource tiers and retry policy are unchanged; cancelling it remains a failure and its retry can request a larger tier. To defer a pending SPAdes task, Slurm hold/release can be used instead of cancelling it:

~~~bash
scontrol hold JOB_ID
scontrol release JOB_ID
~~~

Hold applies to pending jobs; holding a running job does not free its resources. Release held work in time for the parent Nextflow job to finish within its walltime. No scheduler changes or cancellations were performed locally.

## Verification and next remote checks

Only source text, configuration, Git diffs and the supplied logs were inspected. git diff --check passed. No Nextflow run, syntax execution, assembler, biological analysis or test pipeline was executed.

On the next cluster run, check:

1. Complete shared stores are reported as stored/skipped, without download submissions.
2. Hifiasm tasks are cached where expected.
3. FCS cleaning is submitted to normal rather than the high-memory queue.
4. Slurm stdout retains its dashboard but omits input-incompatible assembly branches and disabled QC stages; obsolete-selector warnings are absent.
5. Inspect the resource trace for download/cleaning memory and retain any new failure logs.

This patch addresses the current restart blockers and logging/database request; it does not implement the later pangenome rework or general independent-sample failure recovery.

## Supporting tool documentation

- [Nextflow process directives](https://docs.seqera.io/nextflow/reference/process): storeDir reuse is based on declared outputs.
- [Nextflow 23.10.1 CSV parser tests](https://github.com/nextflow-io/nextflow/blob/v23.10.1/modules/nextflow/src/test/groovy/nextflow/splitter/CsvSplitterTest.groovy): synchronous parsing with the same built-in parser, without implementing a second CSV parser.
- [NCBI FCS-GX quickstart](https://github.com/ncbi/fcs/wiki/FCS-GX-quickstart): database payload list; separate database screening and action-report cleaning.
- [GetOrganelle database configuration source](https://github.com/Kinggerm/GetOrganelle/blob/master/Utilities/get_organelle_config.py): database initialization/check behavior.
- [Slurm scontrol](https://slurm.schedmd.com/scontrol.html): hold and release semantics.

## Input-aware workflow registration update

Input validation now finishes synchronously before any process calls are registered. It still uses Nextflow's own CSV parser, preserves duplicate/path/mate validation, and builds the same sample metadata and read tuples. Accepted samples alone determine the branch flags. The validation-only entry point uses the same parser and report workflow.

Conditional calls now cover the assembler branches, Redundans/Pilon, purge_dups, contig Inspector, read conversion/trimming, the two organelle input branches, Hi-C scaffolding/gap filling, Teloclip, harmonization availability, chimera work, read QC and assembly QC stage aliases. BUSCO downloading and meryl building are also guarded by their settings. No workflow hierarchy or process alias was renamed.

Static expected cases (not executed):

| Accepted inputs/settings | Branches omitted |
| --- | --- |
| HiFi assemblies with no shotgun reads | SPAdes, short-read subsampling/counting/trimming, Redundans/Pilon, GetOrganelle/MITOS branch, short-read QC |
| Short-read assemblies only | HiFi conversion, hifiasm, purge_dups, Inspector, MitoHiFi, HiFi QC, HiFi scaffolding/gap filling/Teloclip, harmonization/chimera work |
| HiFi assemblies without Hi-C | Hi-C trimming/QC/mapping/scaffolding/gap filling; contigs continue to finishing |
| Mixed inputs | Both applicable assembler branches remain |
| Hi-C-only skipped row | Does not enable branches |
| qc_mode=final_only | Intermediate assembly QC aliases |
| qc_mode=none and pangenome disabled | BUSCO download, read/assembly QC and meryl DB building |

In the current mixed dataset, Sde-CMat_061 enables SPAdes and Redundans, so those should still appear. For an actual HiFi-only sheet they should be absent.

An enabled branch whose eligibility depends on computed results can still show no tasks. For example, organelle filtering needs recovered organelle sequence, and MITOS depends on resolved organelle type/recovery. This change omits branches known to be impossible from validated inputs/settings; it does not predict downstream biological results.

The changes were inspected as text and with git diff --check. Nextflow syntax/runtime behavior remains to be verified on the cluster; use the validation checkpoint before resuming assembly.
