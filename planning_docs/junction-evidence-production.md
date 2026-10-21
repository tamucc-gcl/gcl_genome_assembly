# Production chromosome-junction evidence

The normal CHIMERA workflow now runs MISASSEMBLY_JUNCTION after MISASSEMBLY_ASSESS
and before MISASSEMBLY_REVIEW. It uses staged, checksum-checked pipeline inputs;
there is no project focus sheet, cache-directory search or known cut coordinate.
Previous assessment tasks can remain cached. The mapping/discovery/assessment helper is unchanged, so adding or revising this
supplement does not itself change its task hashes. The new stage follows the existing
misassembly_assessment resource/environment configuration and runs per assembly.

## Evidence and decisions

- Reproject retained assembly alignments into chromosome and unique-locus tracks.
  Build masks from at least two eligible independent individuals; both haplotypes
  of one individual contribute one vote and the focal individual is excluded.
- Competitively remap locally selected HiFi reads against complete assemblies of
  the focal individual. Comparable alternatives in the same haplotype invalidate
  confident placement; a homologous placement in another haplotype alone does not.
- A credible bridge requires MAPQ >=30, identity >=95%, >=80% aligned read sequence,
  and no emitted comparable same-haplotype alternative. Record every placement.
- Count credible molecules reaching >=1 kb assigned uniquely to each chromosome
  within a single uninterrupted chain. Report anchor distance and read-length reach.
  Zero arm bridges with unreachable anchors is not an absence test. N gaps retain
  the original N-aware spanning assay; uninterrupted counts are not used against them.
- Report opposite-side islands with chromosome-assigned versus locus-unique bases,
  individual eligibility and largest unique block. A small island is not a cut veto.
- Positive unique-arm bridges oppose a site-level cut recommendation and lower that
  alternative's review priority. Missing new measurements do not become zero or a
  new cut recommendation. All cutting continues through explicit manual review.

## Contact comparison

Rescan retained pairs once per assembly, preserving library identity. Produce exact
5 kb counts, 10/25 kb local maps, fixed-flank scans, and unique-arm contacts. Existing
100 kb / 300 kb-flank counts must replay exactly wherever geometries overlap.

Select nearby unique 5 kb anchor bins with at least five within-bin contacts. An
unobservable nearest bin can be replaced by an observable bin farther along the arm;
retain exact selected intervals. Raw zeros remain in the scan and count matrices.

For each arm-crossing bin pair, compare same-chromosome bin pairs with separation
within 10% or one bin and each end's diagonal contact depth within a factor of two.
Scale their mean count by relative diagonal depths. Require at least five controls
for every candidate pair and positive total expectation; otherwise report no baseline.
This is a descriptive distance/depth comparison, not a significance or phasing test,
and is not used to authorize automatic cutting. Overlapping controls are correlated.

## Published outputs

Each assembly report includes junction placement checks and per-library arm counts,
plus 25 kb contact maps linking 10 kb views and exact 5 kb matrices. The `junction/`
subdirectory contains read placements, reach, chromosome tracks/islands, contact
counts, provenance and matrices needed to reproduce the measurements. Diagnostic
SAMs and duplicate read FASTAs are removed; 50/100 kb and pooled products are not
created in the production supplement. Plot generation obeys chimera_plots; evidence
generation does not depend on plotting. No events produces a short not-needed status.

The standalone investigation imports these same assays rather than maintaining a
second implementation. Its own cache-input preparation remains a diagnostic utility.

## Validation

General unit tests cover control matching, observed zeros, missing controls,
read quality/ambiguity, reach, single-individual exclusion, N-gap handling and reports.
`tests/junction_production.nf` exercises the production module with real alignment
programs and synthetic positive/negative junctions. The cohort workflow integration
fixture exercises staged singleton inputs, unavailable evidence, manual cutting,
finishing, refreshed names, final mapping and the generated report.

Local acceptance: 133 focused tests passed; real minimap2/samtools positive and negative
junctions passed with one Hi-C library; both new tasks resumed from cache. Full CHIMERA
channel/report/finishing/final-mapping integration passed with a manual cut and a single
library (external-tool doubles for workflow wiring, separate from the real-tool assay test).
Crest remains the real-data acceptance run.
