# CTlk YaHS batch 1511217 review

All four attached jobs completed with exit status 0, using YaHS 1.2.2. The two treatments used the same corrected round-1 scaffold FASTA and filtered Hi-C BAM per haplotype. Input component names below belong to that FASTA; they are not interchangeable with newly ranked output scaffold names.

## Result and decision

Neither enabling scaffold error correction nor combining it with a 10 Mb maximum resolution repaired the chromosome-scale composite scaffolds. Do not adopt either treatment as the solution to these fusions. Continue the already-running hifiasm comparison; no further YaHS parameter sweep is justified by this batch alone. Survival through YaHS is not sufficient evidence to retain a disputed junction.

Full-resolution and capped final AGPs are byte-identical within each haplotype. Both full-resolution runs ended after round 8 at 10 Mb, before attempting 20/50/100 Mb. Therefore the cap was not actually exercised relative to the full run: these results do not establish what the larger resolutions would do.

| Suspect input scaffold | Input length | Full-resolution EC | 10 Mb capped EC |
|---|---:|---|---|
| hap1 scaffold_1 (chr7/9/12 composite) | 138,446,979 | Entire input retained as one component | Same |
| hap1 scaffold_5 (chr4/14 composite) | 74,765,583 | Entire input retained as one component | Same |
| hap2 scaffold_5 (chr7/12 composite) | 83,430,767 | Entire input retained, reverse orientation | Same |
| hap2 scaffold_7 (chr4/14 composite) | 82,879,935 | Entire input retained, reverse orientation | Same |

Because each suspect input survives as an uninterrupted component covering its entire length, all internal reference transitions and pre-existing gaps survive. Reversing an entire scaffold does not repair its internal joins.

## Changes relative to retained foundation baseline

Baseline: comparisons/local-review-1506189/comparisons/chimera-foundation-01/assembly/scaffold/yahs_round2/*_round2_scaffolds_final.agp.

Layouts were compared by input component names, intervals, order and orientation, normalizing whole-scaffold reversal and ignoring output renumbering.

Hap1: the only changed layout is baseline output scaffold_15. Its 34,003,065 bp input scaffold_15 remains separate, while the attached 103,114 bp segment of input scaffold_56 and 200,136 bp input scaffold_32 become separate scaffolds. Scaffold count rises from 362 to 364; total output length falls by 200 bp, exactly two 100 bp scaffold gaps. AGP component bases remain 953,569,870.

Hap2: an 8,615,447 bp input scaffold_19 is detached from input scaffold_4 (83,503,982 bp). A 297,584 bp input scaffold_35 is detached from input scaffold_15. The joined pair input scaffold_29 + scaffold_30 (430,410 + 410,834 bp) is attached to the 90,340,288 bp input scaffold_1. Scaffold count rises from 345 to 346; total output length falls by 100 bp. AGP component bases remain 1,065,349,276. These are additional assembly changes requiring evidence before promotion, not demonstrated improvements.

The initial contig-EC logs report 3 breaks in hap1 and 5 in hap2. Their initial AGPs show split input components hap1 scaffold_12/54/56 and hap2 scaffold_16/28/32/34/58, none of the four principal composite scaffolds. Some components have multiple cuts; the log count should not be interpreted as the total number of final separated junctions. Initial contig EC also ran in the baseline and is not uniquely attributable to enabling scaffold EC in this experiment.

## Implication for the core pipeline

Round-2 YaHS accepts the pre-existing composite scaffolds as whole input components even with EC enabled. The next decision must use the hifiasm treatments and retained graphs to distinguish library-dependent assembly paths, then combine exact local junction support with chromosome correspondence across assemblies. Comparative discordance must remain an independent trigger for adjudication; native EC alone missed all four headline cases.

No production defaults or sequence cuts were changed by this review. Archives omit FASTA/BAM, so component-base accounting was verified from AGP rather than by sequence checksums.

Reproducible local layout audit: comparisons/local-yahs-batch-1511217/review_layouts.py. Extracted experiment files: comparisons/local-yahs-batch-1511217/comparisons/ctlk-yahs-batch-1511217/.
