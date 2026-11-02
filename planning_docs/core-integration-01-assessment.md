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
