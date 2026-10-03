# Junction adjudication checkpoint — 2026-10-03

Status: diagnostic implementation prepared; remote tests and current HiFi evidence
review pending. No automatic cut policy has been expanded. Keep `--chimera_break false`.

## Implemented in this pass

- Native transition tables include contiguous chromosome-arm support measured as
  union aligned bases, fractions of assigned sequence, and an out-and-back label.
  Unaligned spans and duplicate alignments do not inflate support. No new absolute
  chromosome-size threshold or sample-specific rule is introduced. Fractions are
  descriptive, not a calibrated major-arm classifier.
- Older-join transition tables distinguish coordinate eligibility, structural
  assessment, read assessment, gap origin and proposed action. A recovered gap may
  be precisely located while the biological decision remains REVIEW. Read status
  explicitly says it has not been assessed at this stage.
- Hi-C evidence includes a local search over the transition plus 2 Mb either side,
  and a baseline from the next 2 Mb on either side. Only finite full-window profile
  values participate. Missing or zero baselines are reported explicitly. These
  overlapping windows are descriptive, not independent statistical replicates.
  The search radius is a diagnostic default, not a cut tolerance or a chromosome
  size threshold. Existing scaffold-wide metrics remain available for comparison.
- Global lowest-window selection now excludes non-finite bins even when fewer
  than five valid windows exist.

## Validation and remaining work

Run on the cluster after the active run finishes and code is synchronized:

```bash
export PYTHONDONTWRITEBYTECODE=1
python3 -m unittest discover -s gcl_genome_assembly/tests -p 'test_*.py'
```

Added tests cover union support, short out-and-back excursions, local versus distant
Hi-C minima, unavailable local evidence, and eligible recovered locations that do
not yet authorize cuts. Existing tests cover wrong nearby gaps, multiple gaps,
changed source gaps, reversed source flanks and complex-scaffold vote restrictions.
Tests have not been run locally, following the repository's remote-testing workflow.

Collect the currently running HiFi context results BEFORE any diagnostic rerun:

```bash
bash gcl_genome_assembly/scripts/comparisons/collect_chimera_review.sh \
  JOB_ID genome_assembly data/assembly_samplesheet.csv data/hic_readsets.csv
```

Replace JOB_ID with that run's numeric ID. Do not rerun mapping just to obtain these
additional descriptive fields: first inspect the existing read-context output.
Changed join-table contents can invalidate downstream context tasks despite unchanged
FASTA sequence; separate reusable mapping from interpretation before requesting
another context run if this would repeat an expensive mapping.

Next, join read-context evidence to these coordinate-bound junctions, assess unique
anchors and whether the available reads could span each junction, and compare arm
patterns with historical positive and negative examples. Then implement and test
per-junction authorization for complex scaffolds and validated older gaps. Neither
absence of a spanning read nor presence of a telomere motif is sufficient alone.
Historical coordinates are test context only and must not be applied to new FASTAs.

A manual cut through contiguous sequence still needs its own explicit, checksummed
route and rationale; the existing gap-only route must not be bypassed or given a
fabricated AGP gap. Full automatic adjudication remains unfinished until the read
evidence and those policies are validated.
