# Chimera integration test

From the Crest project root, after syncing source:

```bash
mkdir -p logs
sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_integration.sbatch ctlk_analysis
```

Three jobs run concurrently with no array concurrency cap. Each uses 2 CPUs, 8 GB and a 30-minute limit. They invoke the actual CHIMERA_ADJUDICATE and BREAK_CHIMERAS Nextflow modules and verify the resulting files. Outputs are isolated under `comparisons/chimera-integration-JOBID/`; production assemblies are not inputs or outputs.

Expected results:

| Test | Expected |
|---|---|
| unresolved | Two original records retained; zero selected actions |
| eligible | Synthetic complete measurements produce one automatic gap-edge cut; three records |
| manual | Explicit reviewed action produces the same two neutral pieces; three records |

All tests require 14 total bases, unchanged unrelated sequence, matching FASTA/name-map records, and the expected gap Ns. Split pieces must have neutral names, unplaced class and forward orientation. The applicator additionally writes source-coordinate lifts, parent reconstruction verification and selected-action provenance. Successful jobs publish `results/verification/PASS.json` and have `exit_status.txt=0`.

Return the three `*.review.tar.gz` files in `comparisons/chimera-integration-JOBID`. These include Nextflow logs, trace, execution report, decisions and verification. Work directories are excluded.

This batch checks process staging, imports, channels, publishing and action semantics. It does not exercise the whole CHIMERA workflow, chromosome-transition discovery, expensive read mapping, Hi-C controls, native graphs, chromosome reassignment or final finishing. Synthetic positive measurements test application of the policy; they do not establish calibration or automatic eligibility of H01. Current real-data collectors still cannot automatically authorize a gap cut.

After this software batch passes, use isolated CTlk assessment inputs for the full evidence/workflow integration: both haplotypes plus current peer assemblies, retain all candidate transitions, compare the proposed H01 manual correction to the already validated reconstruction, and ensure the five supported internal paths remain uncut. That run needs its own complete input manifest; do not rerun the production pipeline merely to exercise these small module tests.
