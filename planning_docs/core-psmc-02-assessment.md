# Own-reference PSMC pilot review

The supplied `core-psmc-02.review.tar.gz` passes execution and production wiring
checks. Exit status is zero; 935 tasks were cached and 91 completed, with no
failures/retries or Nextflow warnings. The core integration check passes,
software-version issues are empty, and report validation passes 644 local links.
Cactus construction remained cached. Four changed-reference Dipcall calls reran;
CMat's call stayed cached and its observed fit is byte-identical to the first pilot.

All five individuals use their own verified final hap1 reference. Reference
checksums match input QC. Each has one complete observed fit and ten complete
bootstraps reaching iteration 25; iteration 20 is exported. All 50 resampling
records are present, with ten distinct derived seeds per individual. All manifest
fit hashes agree with the retained models. Scaling remains unset.

## Data retention

| Individual | Selected non-gap scope, Mb | Dipcall confident, Mb | Complex-mask loss, Mb | Callable before conversion, Mb | Callable on retained scaffolds, Mb | Retained scaffolds |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| CBau | 1066.60 | 251.63 | 22.54 | 229.09 | 215.58 | 15 |
| CLim | 969.86 | 214.40 | 19.04 | 195.36 | 108.16 | 7 |
| CMat | 1055.68 | 244.22 | 19.89 | 224.33 | 165.05 | 10 |
| CPla | 749.34 | 206.18 | 16.11 | 190.07 | 160.66 | 103 |
| CTlk | 933.04 | 200.56 | 17.94 | 182.62 | 92.37 | 6 |

Simple-genotype and duplicate-position losses are zero for all five. The
sequential loss counts balance exactly. The main initial loss occurs in the
upstream Dipcall confidence mask: 72–79% of selected sequence is excluded before
our complex masking. This identifies the responsible stage, not the reason an
alignment failed confidence criteria. The archive contains final diploid BEDs,
but not each haplotype's separate BED/PAF needed to separate missing coverage,
short/low-MAPQ alignments and competing placements. Do not label these losses as
proven repeats, assembly defects or biological absence.

The fixed upstream 20% whole-sequence converter rule causes additional losses,
especially CLim and CTlk, which retain seven and six scaffolds. It is intentional
upstream behavior, not a converter bug. About 53–78% of bins on retained inputs
remain unknown. There is substantially more input than in pilot 01, but scope
selection remains relevant when comparing individual histories.

## Figures and model information

The combined R figure displays all five observed histories in one species/scaling
panel despite distinct own-reference IDs. Companion panels retain bootstraps.
Rendered PNGs were inspected. All 3,190 combined curve rows exactly preserve
their corresponding individual exported values and identities. Bau's extreme
recent-time spike from pilot 01 is absent with the corrected reference/scope.
Recent-time bootstrap uncertainty remains substantial, especially CPla.

Observed parameter groups below the upstream advisory of approximately ten
expected recombinations: CBau 0, CLim 1, CMat 0, CPla 0, CTlk 1. These are
information diagnostics, not proof that demographic histories are unbiased.

## Disposition and bounded follow-up

Reference policy, scheduling, masks/audits, seed publication, exports, figures,
links and core cache preservation worked as intended. Scientific input-method
acceptance remains open. Next inspect each cached Dipcall haplotype BED/PAF to
attribute confidence-mask losses, alongside a single own-reference complete-HiFi
route comparison. Reuse existing Dipcall mappings; do not change aligner settings
or relax masks before attribution. Existing histories are exploratory outputs,
not validated differences in population history.
