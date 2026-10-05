# Junction assessment review: job 1506302

Source: `junction-assessment-20261005-084623-1506302.tar.gz`, reviewed 2026-10-05.
The job exited 0 and its status is SUCCESS. Both CTlk haplotypes have measured
HiFi and Ex2/Ex3 Hi-C evidence; all twelve sequence comparisons completed.
The assessment hashes agree with the investigation registry.

## Findings

The eleven tracking intervals resolve to twelve suspect physical gaps and five
internal-contig transitions. Eighteen additional physical gaps are nearby
comparison joins, not known-correct positive controls. All thirty recovered
gaps were already present in round 1, using literal full-scaffold correspondence.
Round 2 therefore did not introduce these particular physical joins.

Coordinates below are zero-based, half-open. Counts are crossing UU pairs in
100 kb flanks immediately outside each gap, reported as Ex2 / Ex3. They are
descriptive measurements, not calibrated error probabilities.

| Tracker | Current location | Physical gap(s) | Ex2 / Ex3 crossing pairs | Interpretation and next action |
|---|---|---|---|---|
| J03 | hap1 scaffold_1 | 116057447–116057547 | 0 / 0 | Strong priority for misjoin review. Nearby comparison gaps have 34–41 / 68–95 crossing pairs. Assess together with J04 and the intervening component. |
| J04 | hap1 scaffold_1 | 116171488–116171588 | 4 / 19 | Weaker contact support than those same comparison gaps. Flanks of the roughly 114 kb inserted component are not yet validated as correct adjacencies. |
| H01 | hap1 scaffold_1 | 63051225–63051325 | 0 / 0 | Exact round-1 join at the historical chr7_2 / chr6_3 switch. Left-flank contact opportunity is very low: 319 / 51 external contacts versus 8196 / 4622 on the right. Zero crossing alone is consequently not decisive. |
| H02 | hap1 scaffold_1 | Seven gaps, listed below | Variable | This is a complex interval, not one cut site. Prioritize its weak joins and changes of historical correspondence individually. |
| H03 | hap2 scaffold_5 | 54123576–54123676 | 18 / 75 | Ex3 is within the nearby comparison range of 61–118; Ex2 is below the comparison range of 34–52. Historical separation alone does not justify breaking. |
| H04 | hap2 scaffold_7 | 46160478–46160578 | 33 / 86 | Counts are within nearby comparison ranges of 30–90 / 75–178. Not currently supported as a break by these measurements. |
| J01 | hap1 scaffold_1 | No gap in the tracked interval | Not a gap measurement | Every sampled eligible internal position has screened HiFi support; no local zero-support trough. |
| J02 | hap1 scaffold_1 | No gap in the tracked interval | Not a gap measurement | Same result as J01. |
| J05 | hap1 scaffold_5 | No gap in the tracked interval | Not a gap measurement | Local mapping-ambiguity interval at about 35.473–35.489 Mb; inspect alternatives, not an automatic break. |
| J06 | hap2 scaffold_5 | No gap in the tracked interval | Not a gap measurement | Every sampled eligible internal position has screened HiFi support. |
| J07 | hap2 scaffold_7 | No gap in the tracked interval | Not a gap measurement | Every sampled eligible internal position has screened HiFi support. |

## H02: seven distinct gaps

| Start (100 bp gap) | Ex2 | Ex3 |
|---:|---:|---:|
| 94151503 | 42 | 83 |
| 94272911 | 25 | 55 |
| 94361581 | 35 | 49 |
| 94416419 | 32 | 81 |
| 94575157 | 2 | 7 |
| 94667558 | 33 | 71 |
| 95221756 | 3 | 5 |

The last gap is at the historical old-hap2 chr6_3 / chr12_1 transition.
The gap at 94.575 Mb is another weak internal join. Neither should be selected
simply because it is the minimum. The other five joins have appreciably more
crossing contacts. Current and historical alignments also show local orientation
and placement differences within this interval, requiring component-level review.

## What the HiFi results mean

All twelve suspect physical gaps have zero primary molecules bracketing the
two 1 kb anchors. Of eighteen comparison gaps, thirteen also have zero; five
have one each. No physical gap passes the stricter chain/MAPQ/NM screen.
Thus this experiment does not make zero gap bridges a discriminator of error.
The 100 Ns are a scaffold placeholder, not a measured physical separation.
Nearby comparison gaps are not validated true joins, and are not positive
controls for a read-spanning assay.

For broad raw-contig transitions, asking one read to span the entire interval
is inappropriate: the intervals exceed the observed read lengths. The local
1 kb-grid profiles are the useful measurement. At each eligible position the
screen asks for a chain extending at least 1 kb on each side, MAPQ >=20,
NM/alignment-columns <=2%, and no internal indel >50 bp or skipped region.

| Raw-contig tracker | Eligible sampled positions | Zero screened support | Minimum / median / maximum screened support |
|---|---:|---:|---|
| J01 | 324 | 0 | 5 / 16 / 31 |
| J02 | 170 | 0 | 5 / 16 / 37 |
| J05 | 138 | 17 | 0 / 9 / 25 |
| J06 | 179 | 0 | 7 / 17 / 29 |
| J07 | 139 | 0 | 5 / 12 / 23 |

J05 has seventeen consecutive zero-screened-support grid points from 35472804
through 35488804. However, 13–24 primary alignments bracket each of those
positions before the MAPQ filter; none reaches MAPQ20 there. There are no Ns
in those local sequence windows. This is a mapping-confidence trough, not an
absence of mapped reads. Competition between homologous haplotypes and other
repeat placements must be distinguished before interpreting it biologically.
The existing historical continuity findings also remain relevant.

Local support is evidence against a simple unsupported assembly discontinuity;
it does not establish chromosome identity, phase continuity or a biological fusion.

## Peer and graph evidence limitations

An exploratory CIGAR-based inspection used the immediate 5 kb outside each gap,
requiring >=4 kb of aligned query bases per flank. This is a descriptive screen,
not a pipeline threshold or exhaustive mapping test. Under this screen, none of
the eight other-individual targets provided qualifying placements for both
immediate flanks of a suspect gap. Wider, more uniquely alignable anchors are
needed; absent placements do not demonstrate absent sequence.

The same screen confirms the historical endpoints for H01, H03 and H04 on
different old chromosomes, and the final H02 gap between old-hap2 chr6_3 and
chr12_1. Old assembly separation can represent either an old fragmentation or
a new misjoin. Self-alignments are technical checks, not independent support.

Native GFA extraction returned 3 segment records and 2705 A records for hap1,
and 2 segment records and 2676 A records for hap2, with no L records. These cover
the five implicated raw contigs. Purged component aliases were deliberately not
guessed, so this is not complete graph coverage of the scaffolding components.
A records are not in themselves verified read overlaps across a breakpoint.

## Bounded next assessment

1. Prioritize J03/J04, H01, and the H02 gaps at 94575157 and 95221756.
   Use retained alignments and exact component coordinates to inspect end
   mappability, wider unique anchors, orientation/order, and alternative neighbors.
2. Assess Hi-C over multiple flank distances, including end coverage and local
   distance-matched comparisons. Explicitly compare the current arrangement
   with any sequence-supported alternative; account for H01's depleted flank.
   Do not treat Ex2 and Ex3 as independent of the scaffolding that used them.
3. Investigate J05's already-retained alternative read placements to distinguish
   homolog competition from repeat ambiguity. Validate any graph-based bridge
   interpretation against actual overlaps/unique flanks.
4. Keep H03/H04 and the four continuously supported internal transitions under
   review, without cutting them solely because chromosome assignments change.
5. Record per-junction reviewed keep/break/defer decisions and exact coordinates.
   No cut is authorized by this report. Do not rerun the complete assembly or
   full HiFi mapping to perform these follow-ups; reuse this evidence packet.

Primary evidence: status.json, seed_gap_links.tsv, both assembly folders'
selected_gaps.tsv, hifi_support.tsv, hic_by_library.tsv, graph_status.json,
CTLK-J*.support_track.tsv, and sequence_comparisons/peer_*_windows.paf.
