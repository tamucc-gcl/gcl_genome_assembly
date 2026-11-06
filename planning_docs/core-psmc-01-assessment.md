# First PSMC pilot: execution passed, input scope unresolved

The supplied `core-psmc-01.review.tar.gz` completed with 934 cached tasks,
74 completed tasks, no failed/retried tasks, and no core graph rebuild.
Four observed fits and ten bootstraps per fitted individual completed. Tlk was
explicitly not fitted. Mutation rate and generation time remained unset.

| Individual | Callable before conversion (bp) | Callable on retained scaffolds (bp) | Retained scaffolds |
| --- | ---: | ---: | --- |
| CBau | 146629023 | 10153730 | chr11_1 |
| CLim | 125609486 | 16840010 | chr10_1 |
| CMat | 224327982 | 165051153 | 10 scaffolds |
| CPla | 133094123 | 18232212 | chr10_1 |
| CTlk | 95537557 | 0 | None |

The selected non-gap reference scope was 1,055,684,996 bp for every individual.
Thus initial callability was only 9.05–21.25%. The later whole-scaffold losses
match upstream fq2psmcfa's fixed >=20% callable-base rule exactly. This rule is
intentional upstream behavior, not an upstream bug, and is separate from `-g`.
The original dense native fixtures did not exercise this sparse-input boundary.

All five used assembly mode. HiFi mapping depth, GQ and allele-balance thresholds
were not applied to these assembly calls. The archive did not retain the raw
Dipcall confident BED, so it cannot establish whether the main initial loss was
the confident mask or subsequent complex/genotype masking. Do not attribute the
loss to biological absence, or relax quality thresholds, without those counts.

## Implemented follow-up

- Per-scaffold mask audit with additive confidence-mask, simple-genotype,
  complex-variant and duplicate-position losses. Original consensus calculations
  and unknown-region lengths remain unchanged.
- Conversion audit records retained names, fractions and rejection reasons.
- Individual/cohort/main reports distinguish pre-conversion callable sequence
  from callable sequence on retained scaffolds. Execution completion is not
  acceptance of representative genomic scope.
- R/ggplot observed-curve overlay plus individual/bootstrap panels, retaining
  unchanged curve values in a combined table. Species/scaling groups are
  separated. Unfitted individuals are listed without invented curves.
- Explicit bootstrap resampling JSON output declarations fix missing publication.
- Compact archives include the raw caller BED when smaller than 20 MB, in
  addition to mask/conversion audit tables.

These diagnostics require a resumed Crest run. Reference selection has also been
corrected: each individual now uses its own hap1 (or sole assembly), instead of
the shared CMat hap1 reference used in this pilot. Changed-reference Dipcall calls
must rerun; accepted assembly/graph work remains eligible for normal cache reuse.
The combined figures have already been rendered from the supplied
pilot curves for review; they do not make those curves biologically accepted.

## Verification

Thirty-three PSMC contract/masking tests pass with actual pinned pysam/numpy and
fq2psmcfa, including 19.9%, 20.0% and 20.1% callable scaffolds. Native Nextflow
23.10.1 fixtures complete both assembly and HiFi routes, and publication checks
verify combined figures/tables, audits, independent seed records and agreement
with upstream numerical scaling. The own-reference fixture uses two different
individual FASTA identities/checksums and verifies their distinct provenance.
All-skipped and not-fitted cases are also
tested. Actual pilot PNGs are visually inspected. Full-genome loss attribution
remains pending Crest evidence.

## Sources

- [Dipcall confidence-mask criteria](https://github.com/lh3/dipcall#output)
- [PSMC converter sequence-retention rule](https://github.com/lh3/psmc/blob/0.6.5/utils/fq2psmcfa.c)
