# Maintained acceptance helpers

These removable helpers validate the real pipeline or prepare bounded fixtures. They
are separate from production modules and are not required for a fresh pipeline run.
Run cluster helpers from the project root containing `gcl_genome_assembly`, `data`
and `work`. Scientific outputs belong in that project's `comparisons/` directory.

| Helper | Purpose |
| --- | --- |
| `run_core_integration.sbatch`, `audit_core_integration.py` | Cached full-pipeline acceptance against an explicit baseline |
| `prepare_integration_params.py` | Explicit supported settings from baseline provenance; no hardcoded species/tuning |
| `run_chimera_pipeline.sbatch` | Full architecture run and compact review archive; normal resume and symlink publication |
| `prepare_shortread_smoke.py`, `run_shortread_smoke.sbatch`, `shortread_smoke.config` | Bounded or explicitly full shotgun input routing/assembly acceptance |
| `run_identity_reference_checks.sbatch` | Native identity/reference contract checks |
| `run_final_comparison_checks.sbatch` | Native final-assembly comparison checks |
| `run_scaffolding_replay.sbatch`, `scaffolding_replay.nf` | Controlled replay of accepted YaHS inputs/AGP, without new assembly biology |
| `make_sv_benchmark.py` | Generic known-answer fixture generator for deferred SV method validation |
| `submit_psmc_mask_investigation.sh`, `run_psmc_mask_investigation.sbatch`, `psmc_mask_investigation.py` | Removable Bash/Slurm audit of Dipcall alignment/mask exclusions; recorded environments, exact trace paths, cached PAF reuse or PAF-only regeneration; [instructions](../../docs/psmc-mask-investigation.md) |
| `submit_psmc_alignment_presets.sh` | Removable parallel asm10/asm20 comparison against completed asm5 mask controls; fixed confidence rules, paired-mask gains/losses and R/ggplot evidence; [instructions](../../docs/psmc-mask-investigation.md#assembly-preset-comparison-after-the-mask-audit) |
| `submit_psmc_full_presets.sh`, `prepare_psmc_full_presets.py`, `run_psmc_full_preset.sbatch`, root `psmc_from_final.nf` | Removable full asm10/asm20 comparisons through the production PSMC workflow on accepted final inputs; separate Nextflow sessions and unchanged mask rules; [instructions](../../docs/psmc-mask-investigation.md#full-variant-calling-and-psmc-comparison) |
| `run_psmc_reliability.sbatch`, `psmc_reliability.py` | Removable actual-input genotype comparison and deterministic sampled HiFi support using the retained completed mapping, even while calling continues; [instructions](../../docs/psmc-mask-investigation.md#genotype-reliability-on-shared-and-preset-specific-scope) |
| `submit_psmc_length_tests.sh`, `prepare_psmc_length_tests.py`, `run_psmc_length_summary.sbatch`, `psmc_length_summary.py` | Removable parallel asm20 20 kb/10 kb trials, retained 50 kb control, automatic reliability audits and combined compact review; [instructions](../../docs/psmc-mask-investigation.md#full-20-kb-and-10-kb-length-trials-with-asm20) |

Old specimen-specific fusion batches, cache-recovery experiments and obsolete report
launchers were removed. Git history retains them; do not restore them as production
workflows. See [verification](../../docs/verification.md) and
[retirement policy](../../docs/repository-layout.md).

The same integration launcher can enable `--run_psmc true` against an accepted
baseline. It runs the production PSMC workflow and packages compact PSMC evidence;
see [PSMC controls and pilot](../../docs/psmc.md). No separate production PSMC
launcher or project-specific assembly rebuild is required.
