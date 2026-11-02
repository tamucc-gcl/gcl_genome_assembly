# Core integration 01: not accepted

Review date: 2026-11-02. Source: `core-integration-01.review.tar.gz`.

## Execution and fixes

837 tasks completed, 85 were cached. SUMMARY_REPORT failed twice (initial task
and retry), with `object 'has_pangenome' not found`. The acceptance audit did not
run because the production pipeline exited with status 1.

The report renderer now initializes the graph-section switch before using it,
handles missing BUSCO denominators and unavailable input-tier metadata, and sets
a default timezone to avoid an incidental timedatectl/OpenSSL conflict. Graph
report fragments now close into one deterministically ordered report input;
previously the fragment queue could feed separate summary invocations.

Real R execution passed graph-enabled and graph-disabled fixtures and a replay
using this run's actual QC measurements. That replay supplies a minimal manifest
and omits optional metadata; it does not certify a fully integrated native report.
Nextflow 23.10.1 passed empty and multiple-fragment tests using the actual main.nf
collection expression. No alignment or biological rules were changed by these fixes.

## Unexpected upstream recomputation

The earlier promise that accepted upstream tasks would all remain cached was
incorrect. FIND_MITO_REFERENCE changed in the Oct 28 taxonomy safeguards work
after the accepted operational baseline. Its new task hash propagated downstream.
Hifiasm remained cached, but organelle and subsequent stages reran; Cactus took
14h18m50s. These reruns must not be attributed solely to report/provenance changes.
The selected mitochondrial accession remains NC_009588.1.

Nine of ten final sequence dictionaries (names and lengths, not base hashes)
match the accepted comparison packet. CLim hap1 changed from 408 to 1,050 final
records, with total length 990,919,049 to 990,869,851 bp. Its first differing QUAST
stage is scaf.base: 420 to 1,076 records and N50 72,480,792 to 67,444,140 bp.
Earlier contig-stage QUAST measurements match. Equal QUAST summaries do not prove
identical input sequence or mappings. Old/new scaffolding commands, logs, versions,
AGPs and input links are needed to establish the cause.

`collect_scaffolding_context.py` is a removable, read-only diagnostic for that
comparison. It collects small retained task files and records missing context;
it does not reconstruct cache state or rerun expensive work.

## Preserved decisions and graph checks

CTlk hap1 retains the three accepted cuts; hap2 retains one. Both cut verification
records confirm exact parent reconstruction. Their final dictionaries match.
The chosen nuclear reference remains CMat hap1 and its recorded FASTA SHA256 is
unchanged. Cached source paths changed, which is distinct from reference sequence
identity changing.

Graph clipping, sequence accounting and coordinate reconciliation pass. However,
the changed CLim input increased graph source records from 6,182 to 6,824.
Source representation remains REVIEW, with seven chromosome-scale paths below
the 80% retained-path screen, previously five. Added flags are CBau hap2 chr8
(79.685%) and CLim hap1 chr8 (79.862%). Existing flags remain CBau hap1 chr8,
CLim hap2 chr8, CMat hap2 chr8, CTlk hap2 chr8 and chr9_1. These screen flags do
not by themselves establish biological sequence loss. They must not be silently
accepted by rebasing the integration audit.

## Provenance and next gate

The production provenance collector completed: 323 normalized version records,
240 source files, 28 conda environment snapshots and three local image snapshots.
Thirty explicit metadata issues include unlabeled tool-version lines; these are
metadata limitations, not thirty failed tools. The Windows archive extraction
hardlink issue for software_versions.tsv was local to inspection.

Keep category 7 integration acceptance open. First compare the retained old/new
CLim scaffolding records. Then resume the normal pipeline with the report fixes
and review any deliberate change to the accepted baseline explicitly. Do not
launch another unrestricted graph rebuild to diagnose a missing task log.

## Follow-up: retained scaffolding context identifies a skipped resolution

The `core-integration-01-scaffolding-context.review.tar.gz` packet establishes
that all inspected CLim refinement/scaffolding commands are identical, as are
their recorded tool versions. YaHS is 1.2.2 in both runs. Round-one stdout is
identical through input link loading and initial contig correction: 72,220,048
read pairs, 34,904,189 intra links, 37,315,859 inter links and ten initial breaks.
These matching summaries are strong evidence but not byte identity checks.

The first differing scaffolding measurement is YaHS's detected RAM limit:
18.715 GiB previously, 3.038 GiB in integration. At 20 kb, the new run estimates
5.633 GiB required against 3.020 GiB remaining and reports:
`No enough memory. Try lower resolutions... End of scaffolding round.`
The old run executes this resolution. The new run skips it and moves to 50 kb.
Both return success. Final round-one scaffold counts are 420 versus 1,076, with
YaHS N50 72,461,392 versus 67,431,040 bp. The downstream changes are therefore
consistent with a documented memory-dependent branch, not unexplained changes
to the scaffolding command or evidence policy.

Upstream sources for the tested version confirm that this is intentional YaHS
behavior, not a proposed patch to YaHS itself:

- [YaHS 1.2.2 memory option documentation](https://github.com/c-zhou/yahs/blob/v1.2.2/README.md)
- [Memory estimates and fallback](https://github.com/c-zhou/yahs/blob/v1.2.2/yahs.c)
- [Detected available memory](https://github.com/c-zhou/yahs/blob/v1.2.2/asset.c)

The production wrapper now uses `--no-mem-check` in both rounds, so transient
available memory cannot silently choose different requested resolutions. This
does not disable OS/scheduler limits or supply additional physical memory. Real
allocation failures must fail the task. A defensive log check rejects reported
memory failures even if the tool exits zero. The resource request rises from
8 to 16 GB for headroom above the old run's 6.012 GiB observed peak; larger
genomes may require a user resource override. No join thresholds, restriction
enzymes, contig correction, MAPQ or resolution list were changed.

Controlled tests execute the production Nextflow module: normal outputs pass,
the required flag is checked, and a zero-exit simulated fallback is rejected.
The exact native replay entry also passes using controlled tools. These tests
validate orchestration and error handling, not restoration of the real assembly.

Before a full resume, the removable native replay hashes both old/new retained
FASTA and BIN inputs and requires them to match. It then runs the production
YaHS module once using the existing BIN, compares the resulting AGP exactly to
the accepted round-one AGP, and packages the result. No mapping or graph jobs
are launched. The old YaHS execution took 265 seconds, so this should be a small
test plus hashing, scheduling and environment startup; native runtime is not
guaranteed. Its actual result is still pending. A full run after the module fix
will legitimately invalidate scaffolding and downstream tasks; do not promise
that the graph will remain cached.

### First native replay: inputs match, worker bootstrap failed

`core-scaffolding-replay-01.review.tar.gz` confirms byte-identical old/new
decontaminated FASTAs (SHA256
`0c69489b9bb9ba015c33d9105dd1ad4017a8b0feb5fdcafe2fe61852c1becefc`)
and YaHS BIN files (SHA256
`5a62319ee517a8ad1e06f765ba15e4c701890f3138067ba35856df2b960f6f14`).
Both submitted attempts failed before the task script started, with
`conda: command not found` in the worker wrapper. YaHS did not execute, and
there is no AGP comparison result. This is a diagnostic launcher defect, not
evidence that the memory-policy fix failed.

The launcher now loads Crest's miniconda3 module before launching Nextflow and
checks that conda resolves. It also stages the retained BIN using a distinct
basename: real YaHS reuses BIN input, whereas the earlier controlled executable
unrealistically created a new prefix BIN. The original input basename collided
with the production module's final BIN output and copy destination.

The controlled executable now models BIN reuse without creating another BIN.
The exact replay entry passes that case. A launcher-level test starts without
Conda on PATH, supplies it through the module bootstrap, verifies the distinct
input name and small review archive, and passes. Shell syntax and whitespace
checks pass. A successful real YaHS replay remains pending.
