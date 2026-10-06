# CTlk adjudication batch 1511274

## Main findings and action

All eight jobs succeeded. Recommend UNJOIN_UNSUPPORTED for hap1 scaffold_1 gap [63051225,63051325), separating the chr9 and chr7 scaffold blocks. This is a scaffold-adjacency correction, not an internal cut of h1tg000004l and not a declaration that a biological fusion is impossible. No assembly was modified by this review.

Do not cut any of the five broad internal transition intervals at their midpoints. The raw graph provides real overlapping read/path structure across most of them, while the high-specificity anchor screen is non-informative immediately around several cores. Thus native path persistence cannot be treated as a proof of correct chromosome fusion, but neither is a label switch enough to sever the path.

The remaining composites are confirmed as chromosome-scale mixtures; their nucleotide cut positions remain unlocalized by this packet. Gap decisions and internal sequence decisions must remain separate.

## Chromosome evidence

The 1 Mb spaced 25 kb tile screen independently identifies the following blocks across CBau, CLim and CMat haplotypes. Ranges are the first/last successful tile starts, not exact block boundaries. Hits represent overlapping/haplotype observations; counts must not be presented as independent biological replicates.

| Assessment scaffold | Strong tile chromosome labels and ranges |
|---|---|
| hap1 scaffold_1 | chr9 at 29-62 Mb; chr7 at 65-77 Mb; chr12 at 102-125 Mb |
| hap1 scaffold_5 | chr4 at 0-12 Mb; chr14 at 39-72 Mb |
| hap2 scaffold_5 | chr12 at 24-37 Mb; chr7 at 58-83 Mb |
| hap2 scaffold_7 | chr14 at 14-37 Mb; chr4 at 50-82 Mb |

J01 flanks: 48 successful flank placements to chr9 on the left and 17 to chr7 on the right. The unique chr7 placements occur farther right, some beyond the 63.05 Mb scaffold gap. This establishes a chromosome-scale switch but does not localize that switch inside the previously designated J01 contig interval.

J02/J06 have no unique_screen anchors on either side in this assay. J05 has no left-side unique_screen anchors; J07 has no left-side anchors and only two right-side placements to chr4. Failure here is assay non-informativeness, not absence of corresponding sequence. The strict full-span/identity screen and repeat competition can reject fragmented alignments.

The pair around the hap1 113,941 bp inserted component primarily maps to chr12 on both sides. It is not established as a chromosome-fusion boundary. Isolated chr11/chr3 placements are insufficient to override the chromosome-wide context.

## Chr4 piece correspondence

Suffixes do not describe homologous parts. Direct sister placements show:

- Hap1 scaffold_13 (chr4_2) tile at 13 Mb maps to hap2 scaffold_15 (chr4_3), approximately 15.316-15.341 Mb.
- Hap1 scaffold_13 tiles at 34 and 35 Mb map to hap2 scaffold_18 (chr4_2), approximately 0.634-0.659 Mb and 1.563-1.588 Mb.
- Reciprocal tiles from hap2 scaffold_15 at 24 and 27 Mb map to hap1 scaffold_13 near 21.891 and 24.873 Mb.

This is evidence that hap1 chr4_2 corresponds in part to BOTH hap2 chr4_3 and chr4_2. It is not yet a complete one-to-one partition, because unique tiles are sparse and hap2 scaffold_18 lacks accepted reciprocal tiles in this screen.

The separate chr7_2 scaffolds share many sister placements. Several terminal matches reverse orientation; treat that as a local correspondence flag, not an automatically proven inversion or a new cut instruction.

## Graph mechanism

The sequence-bearing primary-unitig and raw-unitig outputs are present. The .noseq counterparts correctly report no segment sequences; these are companion files, not failed biological tests.

In both-library graphs:

| Core | Processed unitig | Raw unitig evidence |
|---|---|---|
| J01 | utg000012l | Full core is internal to raw utg000012l |
| J02 | utg000483l | Two exact raw-segment placements border a local seam; competing divergent raw path also present |
| J05 | utg000340l | Full core is internal to raw utg000353l |
| J06 | utg001917l | Full core maps inside raw utg002163l |
| J07 | utg000671l | Core starts near an overlapping raw-segment boundary; do not infer a cut from alignment endpoints alone |

For J02, exact raw mappings end at padded query 162751 and begin at 169280, corresponding to original h1tg000133l coordinates 4085098-4091627. This 6,529 bp mapping seam is the narrow graph-localization target for further review. It is not itself a validated unsupported edge or cut interval.

The native A records give nonzero overlapping placement coverage throughout sampled core positions: minimum 2 for J01, 1 for J05, 2 for J06 in the raw unitigs (sampled every 1 kb). Processed-unitig minimums are 3 for J02 and 2 for J07. These are graph-derived placement counts based on native offsets/read lengths, not independently remapped unique spanning molecules. They explain why the assembler can construct a local path; they do not prove which chromosome copies those reads connect.

## HiFi assay

31,730 reads were retrieved for hap1 and 11,998 for hap2. Many names have hard-clipped representations in retrieved records, although the runner chooses the longest available representation. Competitive placements and raw CIGARs are retained.

No transition core has a MAPQ>=20 alignment spanning both full 1 kb flank anchors. This is NOT evidence to cut: the cores are 138-324 kb long, beyond typical read spans, and almost all candidate and comparison scaffold gaps also have zero such bridges. Only one hap1 comparison gap and four hap2 comparison gaps have one passing molecule each. The proposed broad bridging assay remains poorly discriminating. Use the retained read placements at narrower graph-localized seams; no new genome-wide HiFi scan is needed to inspect those.

## Hi-C and decisions

At 100 kb flanks, hap1 gap 63051225 has 0 crossings in Ex2 and 0 in Ex3. At 250 kb it has 0 and 2. The flanks map primarily to chr9 and chr7 across independent peer individuals. Wider 500 kb and 1 Mb windows do contain contacts, so this is not complete absence of long-range Hi-C evidence. Low contact participation on its immediate left flank is a mappability/local-structure limitation; the zero counts alone do not justify unjoining. The recommendation combines discordant chromosome blocks, lack of localized support, and the availability of a genuine scaffold gap preserving all component sequence.

Hap1 gap 95221756 has 3/5 crossings at 100 kb, but 153/187 at 500 kb. Its immediate chromosome anchors are weak. Keep it as manual/UNRESOLVED rather than exporting an internal cut or claiming it decisively separates chr7 and chr12.

Hap2 gaps 54123576 and 46160478 have 18/75 and 33/86 crossings respectively at 100 kb, with increasing support at wider scales. The first has chr7 correspondence on both available flanks. Neither is supported as the cause of the composite chromosome transition. Retain these existing gap adjacencies provisionally; this does not authorize retaining the entire chr7/12 or chr4/14 composite as a biological fusion.

The two gaps around h1tg002235l_1 remain a local chr12 arrangement question. Broad contacts span the surrounding region and cannot adjudicate each small insertion adjacency independently. No chromosome-fusion cut is justified there by these results.

## Concrete decision worksheet

See ctlk-adjudication-decisions-1511274.tsv. Recommended correction is reviewable and sequence-preserving; no production edits were executed. UNRESOLVED candidates are not silently accepted as chromosome fusions.

Next work should use the retained data, not repeat all eight jobs: validate the hap1 gap unjoin on an isolated copy; inspect J02's narrow graph seam with the retained molecule placements; and localize chr4/14 and J06 endpoints using fragmented/competing flank alignments rather than imposing full 25 kb anchor matches. A gap unjoin may leave ambiguous homologous sequence at a scaffold end; that is preferable to pretending a broad internal midpoint is a known breakpoint.

## Reproducibility and limitations

Extracted packet: comparisons/local-ctlk-adjudication/comparisons/ctlk-adjudication-1511274.
Local audits: summarize.py, blocks.py, sisters.py, graph_contacts.py, graph_reads.py.
The sparse bin matrix and geometry-normalized profiles are retained for later chromosome-block annotation. No mappability-adjusted likelihood or automatic statistical misjoin classifier has been calibrated in this review. Chromosome labels and graph placement counts are reported at their actual evidential strength.
