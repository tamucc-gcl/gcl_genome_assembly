# Production graph 01: failed fragment handoff

Packet: `core-graph-production-01.review.tar.gz`, reviewed 2026-10-24.

## Execution

Exit 1. The trace contains 840 cached tasks, two completed tasks and two failed
attempts of the same Cactus graph task. The corrected regional executable built
successfully (40.9 seconds elapsed; 15.4 seconds execution). Assembly summary
generation also succeeded; all ten expected assemblies remain finalized. The
retained reference is still Sde-CMat_203_hap1 with 15 chromosome contigs.

Cactus construction failed after 5h 4m, retried, and failed after another 5h 20m.
The retained error is from `cactus_consolidated(chr13_1)`:

```
Could not match contig name in alignment to cactus cap:
id=IND_5364652d434261755f313034.1|chr13_1:38585-42241101
```

There are no completed production graph, sequence-policy summary or FULL/CLIP
audits in this packet. Chromosome-arm recovery and final representation therefore
cannot be evaluated. This is not a timeout or an assembly/chimera failure.

## Cause

The regional assignment tool emits faidx fragment identifiers, for example
`sequence:11-50`. Pinned Cactus extracts that region and sanitizes its FASTA
identifier to `sequence_sub_10_50` (0-based half-open). Its Python helper is meant
to apply the matching conversion to the PAF. The helper uses an unescaped `+`
inside a basic sed expression, so ordinary digit ranges do not match. The PAF
therefore retains the colon identifier while the extracted FASTA is renamed.

The unchanged helper reproduces the no-op locally under GNU sed. The intended
FASTA conversion is explicit in the pinned native sanitizer source:

- [PAF/FASTA helper](https://github.com/ComparativeGenomicsToolkit/cactus/blob/d7e5fc1b51327a5a4cdb494c85b0607e3d81a713/src/cactus/shared/common.py)
- [Native sanitizer](https://github.com/ComparativeGenomicsToolkit/cactus/blob/d7e5fc1b51327a5a4cdb494c85b0607e3d81a713/preprocessor/cactus_sanitizeFastaHeaders.c)

The prior local workflow test executed regional splitting but substituted the
downstream Cactus alignment/export. It did not test real extraction/sanitizing
against fragment PAF names, so its success did not cover this integration boundary.

## Correction and acceptance

The production regional adapter now normalizes fragment **PAF query identifiers**
to the sanitizer's `_sub_start_end` convention after the original-coordinate split
audit. Extraction lists retain faidx syntax. Mapping lengths, coordinates, targets,
tags and chromosome assignment decisions remain unchanged. Alias collisions with
another source fail explicitly. No installed Cactus files or vendor sources change.

Before graph construction, a tiny native check executes the installed samtools
extraction, upstream rename helper and Cactus FASTA sanitizer. It passes the same
regions through the production PAF normalizer and checks identifier/length equality.
A disagreement fails before the expensive graph work. The native check provenance
and normalization counts are recorded in the production policy evidence.

Next test: rerun the normal pipeline into `comparisons/core-graph-production-02`.
Keep diagnostic replay flags off and retain the same cuts and biological settings.
The failed graph is rerun; ordinary upstream cache reuse remains in effect.
Category 5 is still pending successful real construction and FULL/CLIP review.

Local validation after the correction: the 141-test pangenome suite passes across
Windows and Linux (the platform-specific skips are covered in the other environment).
All 12 Nextflow integration scenarios pass, including deliberately incompatible
fragment sanitizing, which fails before the Cactus graph command. A production
resume check also passed with the header-check integration. External Cactus tools
are doubles locally; the new installed-tool check must still pass on Crest, and
actual chromosome-arm recovery is not established by these local checks.
