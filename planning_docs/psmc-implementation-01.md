# PSMC implementation and local acceptance — 2026-11-05

Status: integrated implementation; **Crest biological acceptance pending**.
The [public interface and first pilot command](../docs/psmc.md) are current.
The accepted core freeze remains complete.

## Implemented scope

- Default-off `PSMC_ANALYSIS` in the real main pipeline, independent of graph,
  QC and visualization toggles; verified final shared linear references.
- Paired same-individual assembly route with Dipcall 0.3 and its confident BED;
  complete-HiFi route with minimap2 2.28 and bcftools 1.21 all-site diploid calls.
- Long scaffolds only; sex/autosome identification is not required. Optional BED
  inclusion/exclusion remains available. Unknown sequence is never filled as
  homozygous reference. Explicit GQ versus PL-separation confidence bookkeeping.
- Scope/callability, per-contig heterozygosity, post-fq2psmcfa retention guards,
  reasoned skips/not-fitted outcomes and reference/genotype provenance.
- Native PSMC 0.6.5 observed fits and independently scheduled, seeded block
  bootstraps. Splitfa blocks and the upstream nearest-length stopping rule are
  retained; `-I0` avoids unseeded initialization perturbations. No tool source patch.
- Explicit complete iteration parsing, finite interval/parameter checks,
  per-group expected recombination diagnostics and exact requested bootstrap sets.
- Nullable mutation-rate/generation-time parameters; unscaled first run,
  fractional-generation scaling, R/ggplot PNG/PDF and Markdown/main-report integration.
- Compact review-archive collection includes fits, seeds, QC and figures, omitting
  large BAM/VCF/consensus inputs. Those scientific inputs remain under work/results.

## Verification completed

- General suite: **348 tests, no failures**, five optional native-reader checks
  explicitly skipped in the general environment.
- All **27 PSMC contract/masking tests pass** in the pinned native-reader
  environment, including those five known-answer tests.
- Real Nextflow 23.10.1 native fixture completes both assembly and HiFi routes,
  observed fits, two bootstrap replicates each, exports, R figures and reports.
  A true haploid is skipped. Synthetic fixtures validate mechanics, not biology.
  The assembly route recovers all 1,444 planted heterozygotes over 250,000 callable
  bases. The read route recovers 1,431 (99.1%) over 247,608 callable bases from
  approximately 30x error-free simulated HiFi input; these fixture values do not
  establish sensitivity on real reads.
- Native all-skipped case schedules only plan/index; disabled case schedules no
  PSMC tasks. Poor-callability case retains evidence/reports without fits.
  Zero-bootstrap case produces only the observed fits, not phantom replicates.
- Rescaling to mutation rate 1e-8 and generation time 0.5 preserves cached assembly
  alignments, HiFi mappings/calls, masked inputs, observed fits and bootstraps.
  Only plan/exports/rendering/index work changes. These are test values only.
- Increasing the pilot from two to three replicates reuses the previous observed
  fits and four previous bootstrap fits; only the two added fits are scheduled.
- Upstream `psmc_plot.pl` numerical table exports match both adapter routes for
  iteration selection, terminal-parameter exclusion, nondefault 50-bp scaling and
  a fractional generation time. Upstream gnuplot rendering is not claimed; the
  actual R PNG/PDF outputs are rendered and visually inspected. Negative population
  display values are excluded without changing retained estimates.
- Existing multiple-species main-report checks, PSMC-enabled main-report rendering
  and production entry-point compile/input checks pass. Shell syntax, source/link
  inventory and default inventory checks pass. The additive extension introduces
  29 declared parameters; defaults now total 329.
- Version observation and report-index helpers are separate from scientific
  calculations. Native records pass the production version collector with zero
  issues. Report wording changes do not become inputs to mappings or fits.

Native test artifacts are under the ignored `comparisons/psmc-native-review/`
and local temporary fixture directories. No project-specific identifier, taxid,
mutation rate or generation time is embedded in production code.

## Crest acceptance still needed

Run the documented cached full-pipeline pilot with the five paired individuals,
long reference scaffolds, null scaling parameters and ten replicates. Inspect
actual callability, heterozygosity, interval information, curves, resources,
published links and retained core cache identity. Then compare a selected complete-
HiFi route against its paired-assembly result before considering the two modes
biologically interchangeable. Do not mark category 8 accepted from local fixtures.
