# Assignment review: run 1505984

Source: submitted `chimera-assignment-1505984.tar.gz`, read alongside the retained
HiFi support review. This is a pre-sync baseline, not validation of the newly
edited diagnostics. No biological fusion is established.

## Main findings

- CMat hap1 is the reference; its scaffold_7, scaffold_9 and scaffold_12 map to
  consensus chr7, chr9 and chr12 respectively. The composite names represent
  membership, not a claim that every named chromosome is present in full.
- CTlk hap1 scaffold_1 is 138,446,979 bp. Its reported reference footprints are
  62,475,990 bp for chr9, 33,775,183 bp for chr7, and 42,755,475 bp for chr12.
  Another 38,363,695 bp scaffold is named chr7_2. Thus the large composite is not
  three complete chromosomes: a substantial chr7 portion remains separate.
- Long MAPQ-60 PAF records support chr9 followed by chr7 followed by chr12 along
  hap1 scaffold_1. These chromosome-scale assignments cannot be dismissed as a
  handful of small repeat matches. Individual records are chaining spans, not
  basewise unique sequence evidence.
- The chr9 record spans query [61356185,62941917), reference chr9
  [71498243,73142957), reverse orientation. The chr7 record spans query
  [62618168,64034658), reference chr7 [2797697,4264246), also reverse. Their
  query spans overlap by 323,749 bp: the first diagnostic interval is ambiguous
  between chromosome labels, not a precisely localized join.
- Hap1's chr7 record near the next transition reaches reference chr7
  [39559879,40968664); chr12 starts at [233430,2481961). Hap2 scaffold_5 reaches
  the SAME chr12 interval and then chr7 [38309262,40964151), in reverse
  orientation. This is a shared connection at the scale resolved by the PAF,
  not proof of identical nucleotide breakpoints or independent confirmation.
- Both haplotypes also link chr4 sequence around 44 Mb to chr14 sequence around
  0.47 Mb in the reference, in opposite scaffold orientations. The common
  reference frame and repeated assembly methods remain possible shared causes.
- The small hap1 chr4 excursion at 116 Mb includes a 113,920 bp reverse PAF span
  [116057549,116171469), to reference chr4 [82301655,82442781). It is distinct
  from the tens-of-megabases arms and should not carry equal weight.

## Implications

Do not use local read tiling as evidence of a real chromosome fusion. Distinctive
sequence on both chromosome arms must be linked through the ambiguous region.
Where that region exceeds available read lengths, absence of such reads is not
disproof either. Repeats, a contig misassembly, and a reference/assignment problem
remain alternatives to a biological rearrangement.

The transferred PAFs do not carry cg CIGAR tags, so their outer spans cannot give
base-resolved correspondence through the overlaps. A focused comparison needs
base-resolved alignments and competing placements, using both haplotypes and
additional peer assemblies rather than treating CMat as structural ground truth.

The last-round AGP for hap1 scaffold_1 is a single unchanged component, so the
138 Mb structure was not created by a round-2 join. First-round AGP entries at
the five major diagnostic midpoints fall inside W components, not gaps. This is
suggestive of an earlier contig-level origin; full source-to-assessment sequence
verification is still required before declaring the origin of each transition.

Next bounded investigation: verify the relevant source components; obtain
base-resolved alignments through each major ambiguous transition; assess whether
left/right anchors have competing genomic placements; then evaluate retained
HiFi reads against those anchors. Resolve structural suspicion separately from
cut localization. Keep automatic cutting disabled during that assessment.
