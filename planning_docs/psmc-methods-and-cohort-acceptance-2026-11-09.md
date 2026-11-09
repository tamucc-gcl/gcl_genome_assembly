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

The next acceptance is the normal main.nf run on all eligible individuals, with
own references, ordinary cache reuse and ten pilot bootstraps. A single-individual
manifest inherited from prior tests is explicitly cleared by the cohort launcher.
All expected individuals must appear in planning/results or explicit skips.
All planned individuals get sampled HiFi support from their verified published
final mapping. Concurrent supplementary audits isolate the pilot's 10–20 kb
increment and compare completed HiFi versus assembly genotypes. Method overlays
use original callable scopes. If histories still disagree, matched-scope fits
are a separate unresolved sensitivity check.

Core assembly/QC/graph acceptance remains unchanged. Category 8 is not marked
complete until full-cohort results have been reviewed. Removable acceptance
helpers are documented in docs/psmc-mask-investigation.md; normal fresh runs
need only the standard PSMC params, not these comparison directories.
