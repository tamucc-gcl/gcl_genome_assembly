# CTlk next adjudication batch: parallel evidence, then one correction validation

The tested library combinations, hifiasm post-joining switch, YaHS scaffold EC, resolution cap and minimum length did not repair the principal internal transitions. This does not eliminate every possible assembler parameter. End generic parameter sweeps and adjudicate exact disputed joins.

## Scope

Five internal cores: J01 h1tg000004l 109308-433057; J02 h1tg000133l 4022347-4192059; J05 h1tg000153l 719979-857651; J06 h2tg000028l 660967-839165; J07 h2tg000298l 638671-777533 (raw coordinates, 0-based half-open). These are uncertainty intervals, not prescribed cuts. Include the leading hap1 gap at scaffold_1 63051225-63051325, the 94-95 Mb gap cluster, the two gaps around h1tg002235l_1, and hap2 scaffold gap candidates from the existing packet.

Use current exact assessment FASTAs and provenance lifts for scaffold evidence; raw contig coordinates for graphs. Never attach existing scaffold BAMs to raw-contig FASTAs. Preserve existing controls and previous evidence; rerun only newly required measurements.

## Four independent lanes submitted together

### 1. Chromosome correspondence and candidate localization

One task per CTlk haplotype. Use existing full-reference PAFs and anchor results first; add tiled 10/25 kb anchors at increasing distances through +/-1 Mb around the five internal cores and implicated scaffold gaps. Compare against all six better haplotypes, with chromosome identity from verified name maps. Count independent individuals, not six biological replicates. Test both orientations, retain competing placements, mask/reject ambiguous anchors rather than forcing a chromosome call.

Require multiple non-overlapping uniquely placed anchors for each chromosome block. Resolve chr4_2 in both and chr4_3 only hap2, and chr7 pieces, by reference interval and reciprocal sequence correspondence. Test whether sister junctions have matching chromosome endpoints and orientation. Output block/anchor TSVs, narrowest justified transition intervals, ambiguous regions, and sister/peer alternatives. The lane can disqualify a false label transition before any cutting.

### 2. Exact graph-path mechanism

One task per haplotype, processing both and HiFi-only retained native graphs in comparisons/ctlk-fusion-batch-1511216. Use GFA path/segment records and sequence alignments to trace each suspect contig core through raw/processed unitigs and phased contigs. Produce a small subgraph containing the traversed segments plus immediate competing branches, orientation/overlap lengths, coverage tags if present, and edge/path origins.

Distinguish a junction within one unitig from a graph edge or a later path connection. Export native graph tag values without interpreting coverage as a count of uniquely spanning molecules. If GFA lacks required path provenance, report that specific limitation and export local segment sequences and candidate edges; do not resume historical cache hunting or whole-assembly reconstruction. Graph topology localizes a hypothesis; it does not alone prove correct or incorrect adjacency.

### 3. Competitive local HiFi evidence

One task per haplotype, using reads retrieved from existing matching BAMs over the candidate cores plus adequate flanks. Include locally unmapped/supplementary/low-MAPQ records where retrievable; report ascertainment bias from BAM-based retrieval. Map the same read set competitively to both sister regions, six peer homologous regions, and competing chromosome-piece sequences. Retain supplementary/secondary alignments and identify unique anchors on both sides of each hypothesized breakpoint. Inspect gap/clip positions, phase-consistent spanning molecules, depth changes, and alternative path support.

Evaluate all candidates and matched intact controls together. Do not repeat the failed broad zero-bridge rule. If no uniquely anchored assay distinguishes intact controls, mark that assay non-informative. Where reads cannot span the entire uncertain/repetitive interval, use tiling and graph overlap evidence to narrow the test; no arbitrary midpoint cut. Output per-read evidence, support/alternative counts, exact localized positions and informative/non-informative flags.

### 4. Library-separated Hi-C topology

One task per haplotype. Reuse retained pairs and previous focused contact results; one additional pairs scan per haplotype at most. Examine Ex2 and Ex3 independently and combined across the four composites, both chr4_2 pieces, hap2 chr4_3, and chr7 pieces. Use uniquely anchored blocks from existing correspondence for the initial scan, so this lane can begin without waiting for lane 1. Retain endpoint/block counts to relabel with refined correspondence afterward without rescanning.

Compare distance-normalized contact decay across candidate junctions with intact matched controls; test chromosome-block competing partners rather than raw cross-gap counts alone. Check broad abrupt insulation/off-diagonal block structure and library agreement. Record ambiguous mapping and molecule/library identities. Hi-C supports large-scale placement; it cannot by itself specify an internal nucleotide cut.

## Joint decision, no additional generic grid

Combine all four lanes once. A biologically surprising asymmetric fusion requires affirmative physical/phase support; sibling disagreement increases scrutiny but never automatically causes a cut. Agreement of two sister pair labels is insufficient without homologous endpoints. Reproducibility across hifiasm treatments shares the same HiFi input and is not independent molecular replication.

RETAIN: reliable correspondence and informative molecular/topology evidence support the actual connection, or the apparent transition is explained by ambiguous correspondence.
BREAK_PROBABLE_MISJOIN: independent chromosome correspondence plus informative opposing/support-loss evidence identifies a defensible physical boundary. Preserve every base. Internal sequence cuts require localization beyond a broad reference-transition interval.
UNJOIN_UNSUPPORTED: remove an inadequately supported existing scaffold gap adjacency while retaining all component sequence; this need not disprove a biological fusion.
UNRESOLVED: explicit evidence/localization limitation; provide the minimum necessary manual view and exact interval. Never turn this label into silent approval for exceptional chromosome-fusion claims.

After decisions, run a single correction validation batch on isolated copies: AGP/base accounting, component retention, before/after chromosome correspondence, and Hi-C placement support. Re-scaffold only if needed, with adjudicated adjacency constraints; explicitly test that removed suspect joins are not recreated. The correction validation depends on the adjudication and cannot honestly run simultaneously before its cut/constraint coordinates exist.

## Implementation status

This document defines the next batch. New graph and competitive-read runners still need implementation; the existing junction_focus runner already provides parts of the anchor and library-separated contact lanes. No new jobs were submitted by this plan and no production sequence/default was changed.
