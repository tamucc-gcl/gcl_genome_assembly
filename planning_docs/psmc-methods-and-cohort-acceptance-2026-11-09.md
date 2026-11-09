# PSMC method state, 2026-11-09

The Mat assembly-length pilot supports asm20 / 10 kb as the working default:
767.23 Mb callable (72.68%), all 15 selected scaffolds retained, and sampled
new heterozygotes supported at 452/478 observable sites (94.6%). Both shorter
thresholds preserve all 50 kb baseline genotypes exactly. Defaults are now
ordinary configurable params; there are no sample-specific biological rules.

The complete HiFi pilot succeeded: 853.42 Mb callable (80.84%), 10,078,285
heterozygous SNPs (1.181%), all 15 scaffolds retained, observed fit and ten
bootstraps. The serial calling step took 52 h 43 min. The production implementation
has since been parallelized by whole scaffold, with native serial-equivalence
checks. Its history differs substantially from assembly-based PSMC despite
similar aggregate heterozygosity. Do not claim demographic method equivalence.
The reported 102.64 Mb duplicate-position mask loss in the HiFi route includes
conservative masking around same-position records; this is not proof of duplicated
input files and should be interpreted alongside the indel/multirecord semantics.

The full main.nf cohort run passed: five individuals, observed fits and ten
bootstraps each (55 verified fit hashes/complete iterations), ten unchanged final
assemblies, 940 cached and 86 completed tasks, and the core graph reused from
cache. All-individual sampled heterozygote read support was 95.2–96.8% among
observable sampled sites. The focused 10–20 kb increment retained 61.04 Mb and
all 20 kb genotypes; 437/478 observable sampled new heterozygotes were supported.
Pla lost 1.21 Mb at the existing per-scaffold conversion guard; other individuals
lost none. Functional cohort acceptance is complete.

## Accepted decision

On 2026-11-09 the user accepted closing this implementation pass with both
published methods retained as defensible supported options. Auto routing remains:

- Two finalized haplotype assemblies from one diploid individual: Dipcall,
  with the accepted configurable asm20 / 10 kb defaults.
- One finalized assembly representation from a diploid individual: complete-HiFi
  mapping/calling against that individual's own reference.
- Biologically haploid samples remain unsupported for standard diploid PSMC.

The completed pilot's method-history discrepancy remains an open biological/
methodological question. The shared 655.09 Mb had 95.985% heterozygote-focused
agreement; the same mapped reads used for the read route favored that route at
487/500 sampled disagreements versus 2/500 for Dipcall. This is read consistency,
not independent truth. Do not claim method equivalence or identical accuracy.
A future matched-scope fit comparison is optional follow-up, not required work
for the current pipeline implementation acceptance. Do not schedule more PSMC
investigations automatically merely because this question remains open.

Core assembly/QC/graph acceptance is unchanged. Category 8 is accepted. Removable
acceptance helpers remain documented in docs/psmc-mask-investigation.md; normal
fresh runs require only standard pipeline inputs and PSMC params, not the
project-specific comparison directories.
