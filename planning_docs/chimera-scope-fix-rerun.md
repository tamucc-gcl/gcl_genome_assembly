# Chromosome inference scope fix: next checkpoint

The 1501024 run succeeded, but fallback naming sets were incorrectly exported as
confirmed chromosome membership. That admitted CPla fragments to primary evidence
and produced 10 callable BREAK_CANDIDATE rows. Cutting was disabled.

## Correction

The same scope function now supplies both name maps and candidate tables:

- Successful dropoff inference and scaffold in the selected set: member=yes.
- Successful inference but scaffold outside the set: member=no.
- Threshold fallback, threshold-only selection, or missing inference metadata:
  member=unresolved, with an explicit chromosome_scope_reason.

The naming set itself is preserved. Voting role is not a scope criterion, and no
absolute scaffold-size threshold or sample-specific exception is added. Explicit
threshold-only chromosome selection now supports naming but does not authorize
primary chimera evidence or cuts; it is not evidence of chromosome-scale inference.

Existing detection and cutting consumers accept only member=yes. Synthetic tests
cover failed and successful inference, non-voters, renamed/rescaled scaffolds, and
rejection of unresolved scope by the detector and supplied-file cutter.

## Preserve this completed run

From the cluster project root, with the pipeline stopped:

```bash
mkdir -p comparisons/baselines
mv -T -n genome_assembly comparisons/baselines/single_hic_repair_1501024
```

The command does not overwrite an existing destination. Confirm genome_assembly
has moved before starting a fresh publication directory. Keep work/, .nf/assembly/,
databases and all earlier baselines.

## Sync, test and validate

Push all modified/new files, including tests/test_chromosome_scope.py, then pull on
the cluster. No tests, compilation or pipeline commands were executed locally.
The suite now contains 53 tests; the remote result is required.

```bash
export PYTHONDONTWRITEBYTECODE=1
python3 -m unittest discover -s gcl_genome_assembly/tests -p 'test_*.py'
```

After tests pass:

```bash
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch \
    data/assembly_samplesheet_original.csv data/hic_readsets.none.csv validate
```

After validation succeeds and confirms one Hi-C readset per accepted sample:

```bash
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch \
    data/assembly_samplesheet_original.csv data/hic_readsets.none.csv assembly
```

Keep cutting off and retain the single-library inputs. The launcher uses -resume.
Assembly/scaffolding scripts are unchanged; harmonization and affected downstream
tasks may rerun. This checkpoint still does not validate automatic cutting or
implement distinct-individual per-junction concordance.

## Collect after SUCCESS

```bash
bash gcl_genome_assembly/scripts/comparisons/collect_chimera_comparison.sh \
    comparisons/baselines/single_hic_repair_1501024 genome_assembly \
    data/assembly_samplesheet_original.csv data/hic_readsets.none.csv
```

Upload the resulting comparisons/chimera-comparison-*-reports.tar.gz; retain the
alignment archive. The collector now includes the supplied input filenames.

Expected checks: fallback assemblies show unresolved scope and no primary profiles
or callable cuts; successfully inferred chromosome candidates retain their scope;
the CTlk hap2 interval remains REVIEW with no compatible gap; cached upstream work
and the input libraries remain unchanged. Compact older-gap audits may still cover
all scaffolds. These are acceptance criteria, not results established locally.
