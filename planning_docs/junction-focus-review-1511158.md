# Focused junction review — job 1511158

Source: `junction-focus-20261005-142725-1511158.tar.gz`.
The job exited 0 and status.json reports SUCCESS. All twelve anchor comparisons
completed. An independent comparison of the exported 100 kb contact rows against
job 1506302 found zero discrepancies in crossing pairs, projected UU totals or
left/right external-contact counts.

## Contact evidence

Counts below are Ex2 / Ex3 UU pairs crossing each gap, using symmetric flanks.
These nested windows are not independent tests and do not estimate the physical
gap length. Coordinates are zero-based gap starts; each gap is 100 placeholder Ns.

| Junction | Gap start | 25 kb | 100 kb | 250 kb | 500 kb |
|---|---:|---:|---:|---:|---:|
| J03 | 116057447 | 0 / 0 | 0 / 0 | 13 / 38 | 125 / 246 |
| J04 | 116171488 | 1 / 7 | 4 / 19 | 14 / 64 | 147 / 274 |
| H01 | 63051225 | 0 / 0 | 0 / 0 | 0 / 2 | 38 / 46 |
| H02, internal | 94575157 | 0 / 0 | 2 / 7 | 18 / 36 | 62 / 96 |
| H02, final | 95221756 | 0 / 2 | 3 / 5 | 46 / 76 | 153 / 187 |

The three nearby H02 comparison gaps have 500 kb counts of 165–180 / 256–330.
The 94.575 Mb gap remains depressed at that scale. The 95.222 Mb gap recovers
considerably, particularly in Ex2. H01 also remains depressed at 500 kb, even
though left-flank contact opportunity increases substantially at that scale.
There is no close matched control around H01, so this is a descriptive comparison.

J03/J04's 250/500 kb windows include the other boundary of the roughly 114 kb
inserted component. Their wider-window contact recovery can involve the outer
sequence; it does not separately validate the inserted component's two joins.

## Sequence correspondence changes the interpretation

Only paired anchors meeting the diagnostic >=80% aligned-query screen are
summarized here. MAPQ and structural correspondence are evidence, not proof of
unique homology. Target separations use alignment endpoints and can differ from
the expected reference separation because of biological or assembly differences.

### H01: strongest remaining scaffold-break candidate

- Both flanks map to different chromosomes in old hap2 (chr7_2 and chr6_3) at
  every tested offset: 0, 25, 100 and 250 kb.
- At 250 kb offset, CMat hap2 maps the flanks to different current scaffolds,
  scaffold_1 and scaffold_6, both MAPQ60. An alternative left-flank placement on
  scaffold_6 has MAPQ0 and should not be treated as equally compelling support
  for the current join.
- Neither 250 kb-offset anchor crosses another verified scaffold gap on the way
  to H01. The cross-assembly split therefore is not attributable simply to having
  moved the anchors past another scaffold join.
- Contact support remains weak across scales. The immediate left flank is
  depleted, which limits the interpretation of its zero counts.

Working disposition: high-priority provisional misjoin candidate. The packet
does not yet validate a final manual cut or rule out a real individual-specific
connection. Final adjudication should inspect the wider contact structure and
sequence/assembly-graph context at these exact component ends, not repeat the
whole evidence batch.

### H02 at 95.222 Mb: historical separation is not a sufficient break argument

- Both displaced flanks map in the expected order to current CTlk hap2 scaffold_5
  at offsets 25, 100 and 250 kb, with MAPQ60 on both sides.
- They likewise map in the expected order to old hap1 chr12_1.
- Expected current-assembly separations are 50100, 200100 and 500100 bp;
  current hap2 has 43296, 161762 and 473604 bp. These are not literal identical
  representations, but show consistent regional order.
- Neither side's displaced anchors crosses another verified gap before this
  join. Old hap2 instead separates the anchors between chr6_3 and chr12_1.

Working disposition: defer cutting; regional correspondence supports keeping
this connection for now. This does not certify the exact sequence or phase.
It does show that disagreement with one historical haplotype cannot be used
as an automatic breaking rule.

### H02 at 94.575 Mb: local arrangement remains unresolved

- Hi-C stays weak at wider scales.
- In old hap2, the 25 kb-offset anchors are ordered on chr6_3 but are separated
  by 642232 bp, versus 50100 bp in the current representation.
- At 100 kb offset the historical separation is 855414 bp versus 200100 bp;
  the right-hand anchor passes another current gap.
- At 250 kb offset the anchors have discordant orientation in old hap2 and
  old hap1, but they also pass other current joins.

Working disposition: component-order/copy-content review, not an automatic
chromosome-fusion cut. The 25 kb-offset discrepancy is useful because those
anchors do not pass another gap; its cause still needs resolution.

### J03/J04: outer continuity does not validate the inserted component

- At J04's 250 kb offset, both anchors map in order to current hap2 scaffold_5
  and old hap1 chr12_1, separated by 167541 bp versus 500100 bp in current hap1.
- Old hap2 likewise has ordered chr12_1 anchors separated by 171155 bp.
- The left anchor passes J03 and lies outside the inserted component. This
  supports regional correspondence of the outer material, not the two joins.
- J04's immediate anchors map to old hap2 unplaced_122 in opposite orientations.
  Together with the distance discrepancy, this warrants an explicit examination
  of component orientation and repeated/copy-specific material.

Working disposition: review J03/J04 as one local structural problem. Do not
delete the inserted component or assume that removing only it reconstructs the
peer arrangement; the separation discrepancy exceeds its length.

## J05: correction to the working explanation

Thirty-seven primary reads align at least 1 kb inside the focal interval. The
retained file contains only primary placements for these reads: there are no
reported alternate placements, including none covering 80% of the focal read
segment. Thus the current result does not demonstrate competition with the
other haplotype or another repeat copy. Low MAPQ alone did not establish that.

One primary alignment brackets the complete focal interval plus 1 kb on both
sides: read `m84066_260114_212201_s2/217186588/ccs`, position 35467535–35491335,
MAPQ18, NM30. Bracketing is not itself a verified unique bridge, but this is not
an interval devoid of spanning alignments.

Working disposition: no cut justified here. If needed, inspect retained mapper
tags and the exact alignment/sequence of the spanning read before requesting
additional targeted mappings. Lack of emitted alternatives is not proof of
unique placement; the saved output is not an exhaustive search.

## Next decision-oriented work

1. Prioritize H01's exact component ends for final keep/break/defer review.
2. Resolve J03/J04 and H02 at 94.575 Mb by explicit component order, orientation,
   and copy correspondence against the already identified peer locations.
   Preserve all sequence until the local structure is understood.
3. Defer cutting H02 at 95.222 Mb and J05 on the current evidence.
4. Reuse retained alignments and coordinate manifests; do not repeat full HiFi
   mapping or assembly. A final ruling must distinguish a defensible assembly
   correction from a claim that a biological fusion is impossible.

No executable break instructions are produced by this review.

Sources: status.json, exit_status.txt, assembly_0/contact_scales.tsv,
assembly_0/CTLK-J05-ambiguity.read_summary.tsv and read_alternatives.tsv,
anchors/windows.tsv, anchors/flank_correspondence.tsv, and peer_sources.tsv.
