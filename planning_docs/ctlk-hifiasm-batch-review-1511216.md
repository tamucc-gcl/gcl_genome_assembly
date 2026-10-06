# CTlk hifiasm batch 1511216 review

## Decision

The five investigated local transition cores persist in all five treatments: A-only (Ex2), B-only (Ex3), both libraries, both with post-joining disabled (-u0), and HiFi-only. Neither removing the added library nor disabling post-joining resolves these local transitions. Do not change these production settings as a claimed fix. The tested YaHS changes also failed to resolve the principal composites.

This rules out a requirement for the second Hi-C library for the existence of these particular local sequence paths. It does not establish that the full chromosome-scale scaffold layouts are invariant: these experiments generated contigs, not complete newly scaffolded/harmonized assemblies. Library choice changes phase labels and some contig extensions. Local continuity is not biological proof of a fusion.

## Validation and targeted comparison

All jobs report SUCCESS, hifiasm 0.25.0-r726. The both-libraries treatment exactly reproduces the complete sequences of all five original targeted contigs in single MAPQ60 alignments (matches equal alignment block length). This validates the positive control at the investigated loci.

Both original haplotypes were compared to both new haplotypes; phase labels were not treated as stable homology identifiers. For example, A-only reproduces original hap1 J01 in new hap2, and original hap2 J06 in new hap1.

The origin windows include 100 kb padding around each transition uncertainty interval. The core audit uses the original contig coordinates obtained by removing that padding:

| Junction | Original contig | Core interval, 0-based half-open | Result in all five treatments |
|---|---|---|---|
| J01 | h1tg000004l | 109308-433057 | Continuous MAPQ60 alignment |
| J02 | h1tg000133l | 4022347-4192059 | Continuous MAPQ60 alignment |
| J05 | h1tg000153l | 719979-857651 | Continuous MAPQ60 alignment |
| J06 | h2tg000028l | 660967-839165 | Continuous MAPQ60 alignment |
| J07 | h2tg000298l | 638671-777533 | Continuous MAPQ60 alignment |

CIGAR inspection of the best corresponding alignment confirms every core query base is aligned with M/= /X operations, with no insertion/deletion of 100 bp or larger within the core. Smaller substitutions/indels are permitted. Whole-alignment identities are not local read-support metrics and do not prove a join correct.

In -u0, J02 no longer spans its entire padded window: the original-path alignment ends at 4,283,758, 8,301 bp before the padded-window end. It still covers the entire transition core ending at 4,192,059, plus 91,699 bp beyond it. Thus the full-window screen alone would misleadingly call this a repair. Post-joining changes an extension outside the core; it does not remove the transition.

## Reference checks and limits

The supplied mappings use CMat hap1 as the reference. High-MAPQ blocks on the corresponding HiFi-only contigs still include both chr4 and chr14 for J05 and J07. Thus at least those pair assignments remain apparent without any Hi-C input.

J01's largest reference blocks overwhelmingly map to chr9, while the J02/J06 contigs have highly fragmented mappings dominated by chr12 in the largest blocks. The presence of the old transition sequence is established, but the exact chromosome assignment immediately around every transition requires independent flank correspondence and unique anchors across the six better haplotypes. A rank/label transition by itself must not prescribe an internal cut. The reference PAFs are included in the extracted results for this review; no native read-level graph support is included in the archive.

These observations reconcile the earlier historical confusion: a local path can exist in a one-library assembly while its large chromosome scaffold combination is absent. This experiment establishes the local-path persistence using current controlled inputs, independently of historical assemblies.

## Native graph artifacts

Full native graphs and caches were retained on Crest under comparisons/ctlk-fusion-batch-1511216/{both,A_only,B_only,both_no_postjoin,hifi_only}/, but the review archive deliberately excludes .gfa/.bin/.fa. Logs report corrected-unitig counts and recalled edges; they do not attribute the disputed paths to specific supporting reads or graph edges. They cannot settle breaking versus retaining a join.

## Next finite adjudication batch

Stop the tested library/post-join and YaHS parameter sweeps. Run the following evidence lanes together for the five existing transition cores, rather than another whole-assembly parameter exploration:

1. Unique flank correspondence against CBau/CLim/CMat hap1 and hap2, preserving reference-block coordinates and orientation. Resolve whether each chromosome transition is genuine correspondence or repeat/shared-sequence assignment, and whether sister junctions are actually homologous. Track chr4_2 in both and chr4_3 only in hap2 by sequence, not suffix.
2. Trace each exact original transition through the retained both and HiFi-only native graphs, localizing traversed unitig/edge boundaries and alternate paths. Record whether the disputed interval is internal to a unitig, an overlap edge, or a later path connection. Only boundary-localized read support can justify an internal cut; do not cut at a broad interval midpoint.
3. Reassess localized HiFi support and per-library Hi-C support for the competing paths, using unique anchors and matched intact controls. The previous zero-bridge test failed to distinguish controls, so zero counts alone cannot authorize a break. Sequence-alignment reproducibility across treatments is not an independent replicate of molecular support.

Assign each candidate RETAIN, BREAK_PROBABLE_MISJOIN, UNJOIN_UNSUPPORTED, or UNRESOLVED with explicit coordinate/action and sequence preservation. Existing unsupported scaffold gaps can be unjoined without forcing an unlocalized contig cut. The leading hap1 63,051,225 scaffold gap remains a separate gap-level candidate; it is not one of the five raw-contig transition cores.

## Reproducible artifacts

Extracted files: comparisons/local-hifiasm-batch-1511216/comparisons/ctlk-fusion-batch-1511216/.
Full-window audit: audit_junctions.py / junction_audit.json.
Core/CIGAR audit: audit_cigars.py / core_cigar_audit.json.
Reference inspection: audit_reference.py.
No production defaults or sequence cuts were changed by this review.
