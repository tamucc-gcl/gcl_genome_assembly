# CTlk round-1 YaHS prevention batch 1511225

## Decision

Neither tested change resolves the main chromosome composites. Keep the current minimum contig length for now; do not adopt 100 kb as a fusion-prevention fix. Stop this YaHS cap/minimum-length sweep and use the already-running hifiasm experiments to decide assembly-stage changes and junction adjudication.

## Controlled comparison

All six jobs completed with exit status 0 using YaHS 1.2.2. Each BAM/FASTA reference-name and length check passed. Per-haplotype input paths and metadata agree across treatments. Inputs are the retained decontaminated, unscaffolded round-1 contigs and the matching contig-stage filtered Hi-C BAM. All treatments enabled both contig and scaffold EC and used the same enzyme list.

Both rerun baseline final AGPs are byte-identical to the retained foundation round-1 AGPs. Both 10 Mb capped AGPs are byte-identical to their corresponding rerun baseline. This is a particularly useful control: the current experiment reproduces the original scaffold layouts exactly.

The runs ended after an extra scaffolding round at 5 Mb (round 9). None reached 10 Mb, much less 20-200 Mb. The 10 Mb cap therefore never constrained these runs. The suspect joins cannot be attributed to the unexecuted coarse rounds.

## Principal composites

| Original round-1 scaffold | Baseline/cap length | 100 kb treatment corresponding scaffold | New length | Internal transition windows |
|---|---:|---|---:|---|
| hap1 scaffold_1, chr7/9/12 | 138,446,979 | scaffold_1 | 133,206,902 | Both survive continuously in h1tg000004l_1 and h1tg000133l_1 |
| hap1 scaffold_5, chr4/14 | 74,765,583 | scaffold_5 | 71,311,943 | Survives continuously in h1tg000153l_1 |
| hap2 scaffold_5, chr7/12 | 83,430,767 | scaffold_7 | 69,720,518 | Survives continuously in h2tg000028l_1 |
| hap2 scaffold_7, chr4/14 | 82,879,935 | scaffold_5 | 79,922,528 | Survives continuously in h2tg000298l_1 |

Correspondence uses exact input component identities, not newly ranked scaffold numbers. All five previously mapped full transition windows remain within uninterrupted AGP component intervals in every treatment. This is placement continuity, not independent biological support for these transitions. No new sequence/reference alignments were supplied in the review archives.

The 100 kb treatment removes smaller components from scaffolding participation and changes attachments, while preserving large composite scaffolds. In hap2, original scaffold_5 material also separates into 8,557,524 bp and 1,665,454 bp scaffolds plus smaller pieces. These changes do not remove its chr7/12 transition.

## Gap-level findings

The leading hap1 suspect gap at original coordinate 63,051,225 (h1tg000004l_1 / h1tg000036l_1) persists even with the 100 kb minimum, now at 71,494,816 on new scaffold_1. In the baseline this pair first appears in round 4 at 200 kb resolution, well below the cap.

Both gaps around the hap1 113,941 bp inserted component h1tg002235l_1 persist under the 100 kb minimum. Several small-component adjacencies in the original 94-95 Mb cluster disappear, but the flanks are reconnected within the same 133.2 Mb composite; the final h1tg000376l_1 / h1tg000133l_1 pair persists. The disappearance of an exact pair is not enough to declare a repaired chromosome separation.

Hap2 original gap 54,123,576 (h2tg000248l_1 / h2tg001483l_1) persists under the 100 kb treatment, now on scaffold_7. Original gap 46,160,478 (h2tg002190l_1 / h2tg001271l_1) disappears because the former component becomes separate, but both large sides remain in the same new scaffold_5; the internal chr4/14 transition also persists. This illustrates why endpoint separation and chromosome correspondence must be checked alongside exact adjacency absence.

The audit JSON records each selected candidate pair, first observed intermediate AGP, final adjacency, and endpoint scaffold placements. First-observed rounds describe the sampled exported AGPs, not private internal edge decisions. Pair matching is by component identity; duplicated/split components require interval-level review before a cut is prescribed.

## Fragmentation and sequence accounting

| Haplotype | Baseline/cap scaffold count | 100 kb count | Baseline final AGP N50 | 100 kb final AGP N50 |
|---|---:|---:|---:|---:|
| hap1 | 429 | 848 | 71,708,076 | 70,204,243 |
| hap2 | 412 | 895 | 82,879,935 | 69,720,518 |

N50 here is calculated from all final AGP objects, including small unplaced outputs. YaHS log statistics differ slightly; these values are explicitly the final AGP lengths.

Summed AGP component bases remain identical across treatments: hap1 953,537,721; hap2 1,065,284,217. Total output-length reductions of 46,700 bp and 53,400 bp respectively are explained by fewer 100 bp scaffold gaps. There is no loss in summed component-base accounting. FASTA was excluded from archives, so sequence identity was not checked by checksum.

The increased fragmentation is demonstrated. Whether particular newly separated pieces represent formerly correct chromosome attachments requires reference correspondence; this batch does not establish those breaks as improvements.

## Reproducibility

Extracted archives: comparisons/local-yahs-round1-1511225/comparisons/ctlk-yahs-round1-1511225/.
Audit script and structured results: comparisons/local-yahs-round1-1511225/audit_layouts.py and audit.json.
No production defaults or sequence cuts were changed by this review.
