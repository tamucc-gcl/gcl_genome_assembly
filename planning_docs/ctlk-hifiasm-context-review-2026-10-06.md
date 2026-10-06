# CTlk hifiasm context review — 6 October 2026

Source: user-supplied ctlk-hifiasm-context-02.tar.gz, extracted under comparisons/local-hifiasm-context-02. Static inspection only; no assembly, mapping or tests executed.

## Verified findings

- Original CTlk task is work/b4/4c0468bc10b9ab362064193cb14324, native job 1499288. Version is 0.25.0-r726. The saved command and log agree on Ex2+Ex3 paired trimmed FASTQ inputs and the HiFi input Sde-CTlk_104.fastq.gz. This verifies the command, not input-content identity with the historical run.
- Post-joining was enabled with `-u 1`. Coverage-triggered contig breaking was disabled with `--b-cov 0 --h-cov -1`; `--m-rate 0.75` does not independently enable it. Hi-C unitig-misjoin detection was enabled with `--l-msjoin 500000`. Dual scaffolding was not requested. Telomere motif CCCTAA was supplied.
- Homozygous read coverage threshold was 29 and purge-duplication coverage threshold 37. The initial k-mer-count maximum was at count 10. These measurements use different definitions/stages; their difference does not establish an incorrect coverage estimate. Compare with the retained sample k-mer spectrum and actual local depth before proposing --hom-cov.
- The log reports 18 misjoined unitigs (N50 1,018,765 bp) and 36 corrected unitigs (N50 616,968 bp). It does not list identities or breakpoint coordinates, so these counts cannot be assigned to CTlk's suspect transitions.
- After Hi-C partitioning it reports 2,194,193,346 heterozygous bases and 29,732,448 homozygous bases. These are internal graph classifications, not an estimate that almost every base differs between homologues. They deserve comparison with coverage/k-mer evidence but do not justify changing purging on their own.
- Later graph processing reports 3,020 inserted edges, 41 fixed bubbles, and arc recall for both haplotypes. These are operations on assembly graphs, not evidence that any particular chromosome pair was physically joined correctly.
- The original run wrote corrected-read/overlap binaries and raw/processed unitig graphs. The returned inventory found none in either the published directory or the exact resolved task directory. Recovery by inspecting that task is exhausted; retaining those products was not part of that original output contract. The current module has since added graph/cache outputs for future tasks.
- Runtime was 8,670.793 seconds (about 2 h 25 min) and peak RSS 72.080 GB. These are historical baseline measurements, not a promise of future runtime/resource needs.

The log includes a solver score diagnostic prefixed `wrong::` whose two scores differ by about 0.003 on a scale of 6.5e12. Do not interpret this line as a fatal failure or evidence identifying a misassembly. The task completed and emitted both haplotype graphs.

## What this resolves

We now know the actual settings and that Hi-C-dependent unitig correction/partitioning occurred. We still do not know whether the disputed connections were already in processed HiFi unitigs, selected during phasing, introduced during post-joining, or are chromosome-assignment artifacts. Local historical continuity remains distinct from current chromosome-scale scaffold composition. No new cut is justified by the log alone.

The version-specific command source confirms the 0/1 post-join syntax, coverage-breaking switches and minimum-unitig-size meaning of l-msjoin: https://github.com/chhylp123/hifiasm/blob/0.25.0/CommandLines.cpp .

## Next bounded experiment

First prepare a provenance-bound CTlk-only regeneration from the same current trimmed HiFi and Ex2/Ex3 inputs and hifiasm 0.25.0-r726. Confirm their retained paths/digests; do not attempt to reuse missing scratch symlinks or invent equivalent preprocessing. Keep all recorded graph-generation settings fixed and preserve native graphs, inconsistency BEDs, corrected-read/overlap caches, Hi-C cache files, command, version and logs in an isolated comparisons output. Do not run downstream assembly finishing, YaHS, harmonization or pangenome for this regeneration.

Use the unchanged `-u 1` treatment first to recover interpretable intermediate structure and verify that the existing suspect local sequences/path outcomes are reproducible. This baseline is necessary because the original intermediate graphs are gone; regenerating only a modified treatment would confound setting effects with regeneration differences. Compare contigs by sequence/path correspondence rather than contig IDs or hap1/hap2 numbering. If a local path differs from the original, diagnose input/version/phase or reproducibility differences before assigning the difference to a parameter.

Then inspect the five raw-contig transition regions and their alternatives in regenerated raw/processed unitigs and phased graphs. If post-joining is implicated, produce one `-u 0` variant from validated compatible regenerated intermediates. Keep coverage/purging settings unchanged. If Hi-C partitioning is implicated, use the library A/B/A+B treatments described in the main plan, with treatment-specific Hi-C caches; do not run that full matrix automatically before path inspection. Because the original overlap caches are absent, the first regeneration will require the expensive HiFi correction/overlap stage; later compatible variants may reuse newly retained intermediates after verification.

In parallel, resolve chromosome-block identity from current quality-selected peers and sister sequence. This remains required even if a parameter makes a composite disappear: a shorter contig is not by itself proof that the previous adjacency was wrong.

Stop conditions: identify the responsible graph/path operation with boundary-specific evidence; establish an assignment-only explanation; or document that available read/phase evidence cannot distinguish alternatives. Each leads to a defined assembly decision or limitation rather than an unrestricted parameter sweep.
