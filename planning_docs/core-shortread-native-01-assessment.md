# Native short-read pilot review — 2026-10-27

Packet: `core-shortread-native-01.assembly.review.tar.gz`, run 1553512.

## Outcome

Native assembly stages completed; overall pilot acceptance remains pending a
completion-wiring repair. Nextflow reports SUCCESS, nine completed and six cached
tasks, no failed/ignored/retried tasks. The launcher exits 1 because its required
run summary is absent. This is not a successful end-to-end pipeline acceptance.

| Stage | Status | Task duration | Peak RSS |
|---|---|---:|---:|
| Nuclear subsampling | completed | 37m19s | 1.6 GB |
| Read counting | completed | 1m55s | 56.2 MB |
| SPAdes | completed | 3h27m15s | 156.6 GB |
| Organelle filtering | completed | 24.5s | 3.6 MB |
| Redundans | completed | 5m40s | 10.3 GB |
| FCS-GX screen | completed | 12m45s | 456 GB |
| FCS cleaning | completed | 4.9s | 3.2 MB |
| Finalization | completed | 24.8s | 3.6 MB |

Trimming, genome-size estimation, taxonomy, GetOrganelle and MITOS are cached
successful native tasks. Resource measurements include database mappings and are
not statements about how much additional physical memory a rerun would require.

The revised SPAdes command is K55, automatic coverage cutoff, read correction
enabled, 64 threads and 680 GB requested. Subsampling uses the 15x target and keeps
64,717,862 read pairs from 184,022,746 pairs. Genome size implied by the target is
1,116,535,535 bp. Rasusa's final warning rounds actual/target to the same 15.00x
and says output is unchanged; its earlier selected-pair count and the subsequent
SPAdes ~64.66 million corrected pairs demonstrate reduction really occurred.

## Native organelle result

The wrapper log explicitly confirms verified empty-target pruning; the primary
reads result remains available. Six linear records survive in GenBank, lengths
392, 462, 670, 1,618, 1,205 and 598 bp, totaling 4,945 bp. Annotation contains nine
CDS, eleven tRNA and three rRNA features across those records. The stats table's
23 `gene_count` is the count of annotated gene features, not evidence of 23 unique
complete genes. This is fragmented partial recovery, not a complete mitogenome.
The filter removed zero nuclear assembly contigs according to its status table;
this alone does not validate organelle-removal sensitivity.

## Nuclear output limitation

The final name map records 1,973 sequences, totaling 6,358,369 bp; longest 21,551 bp,
all labelled unplaced. That is far below the genome-size estimate. This bounded
15x/single-K pilot demonstrates tool execution and routing, not a usable complete
species assembly. The packet does not establish which portion of the deficit is
assembly settings, read properties or subsequent reduction. Biological purging
optimization stays deferred. No chromosome-scale or graph-quality acceptance is
claimed from this run.

## General completion defects and fixes

1. `main.nf` gathered harmonization references with `collect(flat:false)`. Empty
   source channels emit no item under that operator. Shotgun-only and other
   unharmonized runs therefore never scheduled `ASSEMBLY_ELIGIBILITY`, blocking
   shared references and the run summary. Use `toList()` to supply the legitimate
   empty reference list. This does not manufacture a reference; selection still
   uses actual finalized assemblies.
2. Testing that repaired path exposed a disabled-chimera input collision in
   `ASSEMBLY_RUN_SUMMARY`: coordinates and report both use `NO_PAF`. Stage them in
   separate input directories so the summary executes with chimera detection off.

Two new integration cases use the actual final-inventory wiring extracted from
`main.nf` and real eligibility, shared-reference and summary processes. Both empty
and present harmonization-reference cases pass, including the collapsed short-read
representation and a graph-disabled singleton reference. Assembly input files in
these local cases are synthetic; no native assembly quality claim follows.

The pilot archive now also includes compact final FAI dictionaries and FCS action
and taxonomy reports, which this packet omitted. The launcher's missing-summary
check now emits an explicit diagnostic. Bash syntax and diff whitespace checks
are required before handing back the resume command.

After syncing, rerun the existing `sbatch ... core-shortread-native-01 assembly`
submission. This should cache the successful native assembly tasks and add only
the missing eligibility/reference/summary work. Categories 2/3 remain open pending
that end-to-end check and their separately recorded remaining contracts.
