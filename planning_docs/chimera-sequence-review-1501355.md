# Sequence-context review: 1501355

Static review of the uploaded diagnostic archive and log, 2026-09-30. Job status SUCCESS; assessment checksum matches the preceding diagnostics. tidk version 0.2.65. No pipeline or biological analysis program executed locally.

## Decision

Do not cut the investigated CTlk hap2 transition on current evidence. Telomere-mediated false overlap is not supported by the tested canonical motif, and sequence comparison now shows corresponding continuity in the sister haplotype. This is enough to reject the proposed automatic cut, but not to declare the whole scaffold correct or prove the evolutionary origin of its chromosome assignments.

## Telomere hypothesis

Across the 500-kb extracted region, the largest 100-bp tidk window contains six motif matches, at regional window end 159,900 (approximately scaffold position 33.6599 Mb). Most positive windows contain one match. The three previously cited clean midpoint-spanning reads contain scattered motif matches: eight for read 98045650, eight for 164236248 and seven for 102176180. These are approximately 14–16-kb alignments, not alignments dominated by long canonical telomere-repeat tracts. These counts do not exclude other repeat classes or all degenerate telomeric sequence.

## Important correction to previous interpretation

The sister comparison places the whole 500-kb query in ordered pieces on scaffold_9, around 13.118–13.636 Mb. One MAPQ-60 record spans query [158981,500000), corresponding to hap2 scaffold_3 [33,658,981,34,000,000), against sister scaffold_9 [13,286,285,13,635,811). It therefore crosses the entire previously defined transition interval [33,660,739,33,830,451). This is a gapped alignment, not proof of perfect basewise identity, but it is much stronger local context than chromosome labels.

The earlier statement that the sister lacked this connection was inferred from scaffold-level labels and broad reference assignments. It is not supported by the direct local alignment and should not be used in adjudication. Sister support is also not an independent individual's replication, and pooled HiFi reads do not assign every molecule to a haplotype.

## Remaining uncertainty

Reference comparison still splits the local sequence among chr7/chr12-associated scaffolds and smaller matches. The independent CMat hap2 comparison is fragmented (124 records; no query span exceeds 15 kb). For example, query [259976,267463) has similar placements on numerous scaffolds, supporting ambiguity for that segment. The independent comparison is not evidence against continuity merely because it fails to produce a long alignment under this preset.

The chr7 reference segment is around 40.7–41.0 Mb within an 81.4-Mb scaffold, whereas the chr12-associated segment approaches the beginning of a 45.9-Mb scaffold. Thus the available reference context is not a simple pair of chromosome ends. Naming, structural variation, repetitive sequence and possible reference assembly problems remain separable hypotheses; none currently authorizes a technical-error cut here.

## General pipeline changes supported by these diagnostics

1. Use direct interval-level sister/peer alignment evidence rather than absence of a composite scaffold name to assess concordance. Never count the sister as an independent individual.
2. Retain raw tidk windows as declared/published evidence outputs, with window resolution and offsets; a motif peak far from a candidate must not be called junction-local evidence without its distance.
3. Offer optional HiFi support diagnostics against a checksum-verified assessment reference, reusing compatible BAMs and retaining alternative alignments. Midpoint-spanning reads alone must not stand in for evaluation of a whole uncertain interval.
4. Separate alignment-discordance detection from permission to cut. Report supported continuity, unresolved repeat context and localized unsupported connections distinctly; retain a conservative REVIEW outcome where appropriate.
5. Use this case as a no-cut regression example. Add a separate positive technical-join example before expanding automatic contig-internal cutting.

These reporting changes can proceed without another CTlk assembly rerun. The biological explanation can remain recorded as unresolved while work returns to the general pipeline revamp.
