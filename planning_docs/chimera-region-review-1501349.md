# Targeted chimera review: job 1501349

Reviewed 2026-09-30 from `chimera-region-20260930-080406-1501349.tar.gz` as static output. No pipeline or analysis was run locally.

## Conclusion

The chr7/chr12 sequence transition in CTlk hap2 scaffold_3 is supported by substantial alignments. A technical misassembly remains a plausible hypothesis, but the available results do not yet establish an exact, justified cut. Retain REVIEW. This is not evidence that the assembly is correct.

## Verified observations

- The pre-finishing assessment FASTA SHA256 is `2cbed25444ba9e50e49aea7751d6068a22da49c89602413f4c5976ecf21acf23`; scaffold length is 76,699,408 bp. AGP/FASTA validation passed.
- Hi-C pairs were rebuilt from retained source pairs through the last-round AGP. The audit reports zero ambiguous ends, zero unmapped ends and zero distance mismatches. This validates that projection, not biological correctness or every earlier coordinate transformation.
- Large chr7 alignments extend to 33,665,002; a large chr12 alignment begins at 33,826,222. A shorter chr7 alignment reaches 33,802,889, so the large-anchor interval must not be treated as a precise breakpoint.
- The final regional N-gap is 33,004,681–33,004,781 (zero-based, half-open). There is no N-gap at the chromosome transition. The round-two join at approximately 30.202 Mb is also upstream. Selecting either gap would not resolve the observed transition at its actual location.
- At 33.7 Mb, cross/combined-within Hi-C ratios are approximately 0.54, 0.65 and 0.70 for the three distance bands. However, medium- and long-range cross-contact means exceed the left-flank within-contact means; the right flank has much stronger contacts. These ratios are exploratory, unbalanced band summaries, not a calibrated breakpoint test.
- The strongest short-range dip near 33.0 Mb coincides with only 1,989 read ends in its 100-kb bin, compared with 5,697–7,803 in several following bins. The transition bin at 33.7 Mb has 4,037. Coverage differences complicate interpreting contact dips as misjoins.
- Companion alignments place substantial chr7 and chr12 sequence on separate scaffolds in CTlk hap1 and CMat hap2. This supports an unusual arrangement in the candidate assembly, but does not distinguish a technical error from a real heterozygous rearrangement.
- Inspector's retained scaffold BAM was aligned to the original round-one `Sde-CTlk_104_hap2_scaffolds.fa`, before Inspector correction. It is not a BAM against the current assessment FASTA. Its input was staged from `work/8c/9e136b913f87e4137296d88f80e42b/Sde-CTlk_104_hap2_scaffolds.fa`. Do not query it with current scaffold_3 coordinates or pair it with the corrected FASTA merely because both were published together.

## Next targeted check, before another full pipeline run

1. Establish the local sequence correspondence between the assessment interval and the exact Inspector input FASTA. Reuse the old BAM only after demonstrating coordinate-compatible sequence across the locus; otherwise map the individual's HiFi reads against the complete assessment assembly, retaining supplementary alignments. Mapping exclusively to the small suspect region would hide competing placements.
2. Examine coverage, confidently anchored reads crossing candidate boundaries, clustered clipping/split alignments, and alternative placements. Use matched nearby controls and distinguish repeat ambiguity from absence of support. A 160-kb alignment interval is not something one HiFi read must span; first localize the putative discontinuity.
3. If evidence localizes an unsupported join, propose a cut with explicit uncertainty and supporting records. If reliable reads support continuity, keep the sequence and report the discordant chromosome arrangement. Low coverage or ambiguous placement remains REVIEW.
4. Only after this experiment should the pipeline gain a general contig-internal adjudication route alongside its gap-join route. Do not weaken the gap checks or add sample-specific exceptions to obtain a cut here.

No full assembly rerun or second Hi-C library is needed for this next diagnostic. Any new comparison scripts and outputs remain under the existing comparison folders.
