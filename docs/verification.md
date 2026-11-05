# Core verification and native acceptance

## Local checks

Use an environment with Python 3.11, NumPy, Matplotlib and Biopython 1.86 for the
general tests. R checks use the pinned reporting/comparison environments. No pipeline
task installs packages dynamically.

```bash
python3 scripts/check_repository.py
python3 -m unittest discover -s tests -p 'test_*.py'
bash tests/run_chimera_channel_reproduction.sh main_compile_checks.nf enabled
bash tests/run_chimera_channel_reproduction.sh chimera_cohort_wiring_check.nf
bash tests/run_chimera_channel_reproduction.sh main_report_checks.nf multispecies
bash tests/run_chimera_channel_reproduction.sh version_record_checks.nf
```

The historical runner filename is retained as a stable test entry point; its default
is the maintained full-cohort wiring check. Integration fixtures use controlled tool
doubles for scheduling/identity checks and do not establish biological accuracy.
The runner's `/tmp/chimera-nxf-runtime` locations are test-only. It can fetch a
Nextflow/Java test runtime; this is not task-time production dependency installation.

Actual retained-data rendering checks are `tests/run_main_report_checks.py`,
`tests/run_junction_placement_plots.py` and `tests/run_misassembly_figure_checks.py`.
The last verifies every retained interval/profile row, close-up coordinates, missing
evidence behavior, independent primary/all-placement scales and detailed-report links.

## Native cached check

From the cluster project root, use the removable
`scripts/comparisons/run_core_integration.sbatch` with sample CSV, Hi-C CSV, accepted
baseline results, reviewed cut TSV and a new output directory. Its baseline settings
are read from the accepted run's provenance and filtered to current configuration
keys; old unused options are not restored. Invocation-specific paths come from the
caller. CLI overrides remain explicit.

It launches the real `main.nf`, uses normal Nextflow resume and symlink publication,
then checks report validation, final assembly/reference identities, graph admission
and Cactus cache reuse. No custom cache reconstruction or deleted-work recovery is
implemented. The post-run audit reports unexpected rebuilds; it is not a scheduler veto.

After this freeze the version-format fixes can rerun tidk and final contact-map
producers plus reports/provenance; those are presentation/metadata tasks, not assembly
or read-remapping changes. Existing cached version files are not rewritten or inferred
from bare numbers. Accept clean new records after those producers run.

Expected output: launcher exit 0, report validation PASS, core integration PASS,
no unintended scientific-output changes, and no new version-format issues.
Use the review archive to inspect actual figures and retained evidence. Review archives
omit full assemblies and graph/index products, so full identity acceptance happens
on the cluster rather than by pretending omitted local files are present.
