# HiFi diagnostic review: 1501351

Static review of the submitted reports and Slurm log, 2026-09-30. No analysis program or pipeline was executed locally. Tables were read and filtered for inspection.

## Mapping validation

The diagnostic reports SUCCESS and the expected assessment SHA256 `2cbed25444ba9e50e49aea7751d6068a22da49c89602413f4c5976ecf21acf23`. Minimap2 2.28 used map-hifi against the complete assessment assembly. Its mapping stage took 3,365 seconds and reported 11.812 GB peak RSS (not total job resource use). Of 1,912,668 primary records, 1,908,541 mapped (99.78%). Secondary and supplementary records are additional alignments, not additional molecules.

## Main finding

The diagnostic midpoint at 33,745,595 has direct HiFi continuity support. It should not be used as a cut merely because chromosome-level alignments change nearby.

Examples of distinct MAPQ-60 primary reads crossing the midpoint, with 1-based alignment starts:

| Read suffix | Start | CIGAR | NM |
|---|---:|---|---:|
| 98045650/ccs | 33,736,075 | 2380M1D12759M | 1 |
| 164236248/ccs | 33,736,841 | 6205M1I4347M1I1538M1I3781M | 3 |
| 102176180/ccs | 33,738,793 | 81M1D1622M1I356M1D5757M1I3411M1D1992M1I584M | 6 |

All have prefix `m84066_260114_212201_s2/`. These are inspected examples, not an exhaustive supporting-read count or haplotype assignment. The first read spans approximately 33.736–33.751 Mb, with substantial aligned sequence on both sides of the midpoint.

MAPQ>=20 primary mean depth in the ten 1-kb windows spanning 33.740–33.750 Mb is approximately 43–48x, with zero uncovered bases. Within the inspected 33.660–33.840 Mb window range, the lowest mean-depth bin is 33.697–33.698 Mb at 10.544x, also with zero uncovered bases. Coverage does not independently prove correct assembly.

Some reads also show recurrent large insertions or supplementary alignments (including placements on scaffold_21). Because reads come from both homologues and the mapping reference contains only hap2, these observations could reflect allelic variation, repeats, or assembly problems. They are not an automatic technical-error diagnosis, especially alongside clean spanning reads.

## Decision and remaining question

Retain REVIEW and do not cut at the midpoint or move the cut to the upstream Hi-C trough/N-gap. The result weakens the unsupported-join hypothesis specifically at the midpoint; it does not prove that the entire broad chromosome transition is biologically correct.

The unresolved question is whether the chr7/chr12 label transition reflects a real rearrangement, ambiguous/repetitive assembly-to-assembly alignment, a haplotype switch, or a different localized error. Before adding any automatic contig-internal cutting rule, examine read support across the complete anchor interval and inspect the corresponding sequences in the sister/reference assemblies. Use the retained BAM for further inspection; do not remap or rerun assembly just to obtain another local view.

Full BAM: `comparisons/chimera-hifi-20260930-082608-1501351/hifi.bam` on the cluster. Report source: `chimera-hifi-20260930-082608-1501351-reports.tar.gz`.
