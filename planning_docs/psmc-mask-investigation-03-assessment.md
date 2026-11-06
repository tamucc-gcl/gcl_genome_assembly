# Dipcall mask investigation 03: alignment-length bottleneck identified

The supplied archive passes: exit 0, ten complete tasks, no missing tasks and
zero disagreement between independently predicted and actual paftools masks.
The regenerated paired masks also exactly match the original pilot-02 diploid
BEDs. All PAFs required regeneration because the original raw intermediates had
been removed with scratch. The diagnostic commands used the recorded environment
and standard PAF/mask Make targets; production data and filters were unchanged.
The actual cohort figure was inspected.

Hap1 self-alignments retain 100% of selected non-gap scope in all five individuals.
Every paired-mask exclusion therefore originates in hap2-to-hap1 alignment.

| Individual | Unique >=50 kb, retained % | Unique 10-50 kb, excluded % | Only <10 kb, excluded % | No observed alignment span % |
| --- | ---: | ---: | ---: | ---: |
| CBau | 23.59 | 32.21 | 24.06 | 19.66 |
| CLim | 22.11 | 33.88 | 23.32 | 20.18 |
| CMat | 23.13 | 32.56 | 23.70 | 20.19 |
| CPla | 27.51 | 32.39 | 18.00 | 20.65 |
| CTlk | 21.50 | 34.70 | 24.93 | 18.11 |

Other categories are small. Overlapping qualifying alignments account for only
0.25-0.42% of scope, and low-MAPQ alignments at least 10 kb for 0.03-0.05%.
This rules out overlap filtering as the dominant direct source of mask loss.
Repeats could still contribute to fragmented or absent mappings; the diagnostic
does not identify those regions as biological absence or prove repeat identity.

In CMat, 18,457 accepted primary alignment records are 10-50 kb and 105,586 are
shorter than 10 kb, versus only 2,165 accepted >=50 kb records. Fragmented hap2
alignments and the fixed 50 kb threshold explain most exclusions.

## Next bounded method test

Dipcall 0.3 hardcodes minimap2 asm5. Minimap2 2.28 documents asm5 for average
divergence not much above 0.1%, asm10 around 1%, asm20 several percent. Current
callable-subset heterozygosity is around 0.4%; asm5 being over-stringent is a
testable hypothesis, not an established cause. Test asm10 and asm20 against the
retained asm5 baseline while keeping the 50 kb/10 kb/MAPQ 5 mask rules fixed,
using one representative individual first. Measure coverage, block lengths,
uniqueness and mask retention. Retain the ongoing complete-HiFi comparison.

If testing a lower alignment-length threshold, update variant-call filtering
and masks consistently. Both paftools call and dipcall-aux samflt default to
50 kb. Enlarging only the BED while retaining the old VCF could incorrectly
impute uncalled differences as homozygous reference. No production default has
been changed based on this review.

## Ongoing HiFi comparator

User-reported scratch evidence: task starts around 12:10, depth table and scope
BED completed around 12:17, compressed all-site VCF is 371 MB with 15:14 mtime.
This is consistent with calling progress, not a completed fit or proof against a
later stall. The current module runs one unsharded mpileup/call stream across all
long scaffolds, including covered homozygous-reference positions, and writes a
compressed VCF. Four requested CPUs do not create four regional callers.
Chromosome-level calling parallelism is a future performance improvement;
leave the ongoing baseline intact while checking CPU/file growth.

Sources: [minimap2 2.28 presets](https://github.com/lh3/minimap2/blob/v2.28/minimap2.1),
[Dipcall SAM filtering](https://github.com/lh3/dipcall/blob/master/dipcall-aux.js).
