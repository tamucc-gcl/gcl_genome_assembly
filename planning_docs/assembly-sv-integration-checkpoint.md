# Optional assembly SV integration checkpoint

## 2026-10-02 pilot failure: delta scratch paths

Run 1503186 completed alignment but ASSEMBLY_SV_CALL failed during SNP identification: show-snps tried to open reference.fa in the previous alignment task's scratch directory. MUMmer delta files retain reference/query paths in their first line; staging the FASTAs into a new scratch directory does not rewrite that header. The pandas FutureWarnings in the excerpt are not the reported fatal error.

The calling process now streams a separate calling.delta copy with its header rebased to the currently staged FASTAs, validates the NUCMER signature and preserves the remaining alignment bytes. It does not edit the cached delta or alignment process. Static review only; the fix still requires cluster validation. A future retry should reuse alignment but rerun SyRI calling; intermediate results from the failed scratch task are not automatically resumed.

Per the latest plan, do not rerun this optional pilot now merely to validate the fix. Proceed with the new-Hi-C production assembly with run_assembly_sv=false once the previous Nextflow controller has exited. Keep its work/cache and failure logs for the post-core extension phase described in [the optional analyses plan](post-core-optional-analyses-plan.md).

## Latest checkpoint: bounded execution and close-out scope

The returned plan selected CMat hap1 against CLim hap1, CLim hap2 and CMat hap2 (15 chromosomes each). CPla is not chromosome-scale. CBau and CTlk remain chromosome-scale but have split chromosome representations; CTlk hap2 additionally has a composite assignment. Skip reasons now distinguish these implementation limits from unrecognized naming. These restrictions do not alter graph eligibility.

`--assembly_sv_queries Sde-CLim_110_hap1` selects one explicit query. Empty selects all eligible queries. Unknown/ineligible requests fail, while other eligible rows become `not_selected` in the plan. Selected task metadata stays unchanged when broadening selection, allowing the pilot pair to resume.

After tests and input validation pass, execute the first full chromosome-set pair:

```bash
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch \
  data/assembly_samplesheet_original.csv data/hic_readsets.none.csv assembly \
  --run_final_hifi_mapping true --run_assembly_sv true \
  --assembly_sv_plan_only false --assembly_sv_queries Sde-CLim_110_hap1
```

In parallel, audit existing final HiFi BAM provenance using the already-tested environment's samtools (not a new mapping):

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_final_hifi_audit.sbatch \
  "$PWD/comparisons/environments/syri-1.7.1" genome_assembly
```

Return its generated `comparisons/final-hifi-audit-*.tar.gz` and job log. It reads final FASTAs for checksums; no BAMs or FASTAs are archived. Successful provenance checks do not validate biological events or raw-read identities.

After the SV run, collect compact reports:

```bash
stamp=$(date +%Y%m%d-%H%M%S)
find genome_assembly/post_assembly/sv -type f \
  \( -name 'sv_plan.*' -o -name 'sv_qc*' -o -name '*provenance.json' \
     -o -name '*sha256' -o -name 'versions.tsv' -o -name 'syri.summary' \) \
  > "comparisons/assembly-sv-${stamp}.files"
tar -czf "comparisons/assembly-sv-${stamp}.tar.gz" \
  -T "comparisons/assembly-sv-${stamp}.files"
```

Return that archive plus matching Nextflow/Slurm logs and pipeline trace for resource/cache review. Retain native calls, delta and BAM files on cluster. Do not mix older published pair outputs with proof of current execution; reconcile reports against this run's plan and trace.

First-release acceptance: successful multi-chromosome alignment/calling, unchanged call preservation, correct coordinates/provenance, explicit skipped coverage and bounded resource use; then expand to all currently eligible pairs. Mark as assembly-derived calls with advisory QC. Later focused work: split/composite chromosome support, repeat/alternative-placement disambiguation, read-supported confidence filtering, and more informative rearrangement diagrams. Do not label this release biologically validated or imply that skipped assemblies have no SVs. Unified Markdown inclusion remains part of the broader report integration.

2026-10-02. Implemented, statically reviewed, awaiting cluster tests. No local pipeline or analysis execution. Short-read work remains parked.

## Contract

`--run_assembly_sv true` enables an independent branch from the shared assembly eligibility manifest. It does not depend on graph execution, QC, final visualization or HiFi reads. `--assembly_sv_plan_only true` is initially the default: only plan/report generation runs. Missing shared reference selection is reported rather than independently choosing a conflicting reference. The current shared manifest may withhold references for singleton cohorts; this initial branch does not add a separate singleton reference policy.

The planner compares the shared species reference against each other suitable assembly, including its partner haplotype. No all-versus-all and no sample-specific rules. It requires complete species inputs, eligible representations, chromosome-scale status, and identical selected chromosome-name sets of unambiguous harmonized `chr<number>_1` sequences. Composite/fragmented/mismatched sets are explicit skips, never silently intersected. This is a conservative implementation boundary, not a claim that biological rearrangements cannot cause mismatches. The plan must show its actual coverage before accepting production use.

All selected chromosomes in a pair are aligned together, permitting between-chromosome alignments within that set. Unplaced sequence is excluded, and missing/chimeric chromosomes are not repaired to satisfy SyRI. Chromosome names remain the final assembly names. Full multi-chromosome performance and calls have not yet been validated; the successful real pilot covered one chromosome pair.

## Execution and output

`ASSEMBLY_SV_ALIGN` extracts selected chromosomes from staged copies, aligns with MUMmer4 default reference-unique anchors, filters using the tested `-m -i 90 -l 100` settings and retains delta/coords/selected FASTAs in work. `ASSEMBLY_SV_CALL` runs SyRI and advisory QC separately. Native syri.out, summary, VCF, QC tables/report, assembly identities and source/selected FASTA hashes publish to `post_assembly/sv/<taxid>/<pair_key>/`. Pair keys are reversible encoded query IDs; `sv_plan.tsv` provides readable names. Plan files publish directly under `post_assembly/sv`.

Environment pins match the successful pilot dependencies: SyRI 1.7.1, MUMmer4 4.0.1, Python 3.11, pandas 2.2.3, NumPy 2.3.5, samtools 1.21. The new environment solve remains untested. Alignment/calling each request 4 CPUs, 64 GB and 96h; each process has maxForks 1. This can still allow an alignment and a call concurrently. No blanket retries for deterministic failures. One-chromosome timing is not a full-genome runtime estimate.

Raw calls are retained. QC status is advisory and read support remains NOT_ASSESSED. No figures are generated by default; the comparison diagnostics remain the plotting tool. Workflow emits report and versions channels for later unified-report integration, and standalone Markdown is already published. The unified pipeline report has not yet been updated to consume this branch.

The QC implementation is now `py_scripts/sv_alignment_qc.py`; the comparison entry point forwards to it, avoiding parallel implementations. QC changes rerun the call/QC task but not alignment; finer separation can follow only if useful.

## Next user-run gate: planning only

After syncing, run tests and normal input validation. Do not launch a second resume into the active main session while HiFi mapping is still running. Once that run finishes:

```bash
export PYTHONDONTWRITEBYTECODE=1
python3 -m unittest discover -s gcl_genome_assembly/tests -p 'test_*.py'
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch \
  data/assembly_samplesheet_original.csv data/hic_readsets.none.csv validate
# After validation succeeds and the active mapping run has finished:
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch \
  data/assembly_samplesheet_original.csv data/hic_readsets.none.csv assembly \
  --run_final_hifi_mapping true --run_assembly_sv true --assembly_sv_plan_only true
```

Keep work/results; final mapping and upstream assembly should resume. Only the new planner should add substantive work (small reports can refresh). Return:

```bash
mkdir -p comparisons
tar -czf comparisons/assembly-sv-plan-review.tar.gz \
  genome_assembly/post_assembly/sv/sv_plan.json \
  genome_assembly/post_assembly/sv/sv_plan.tsv \
  genome_assembly/post_assembly/sv/sv_plan.md \
  genome_assembly/assembly_eligibility.json
```

Also return the corresponding Nextflow and Slurm logs. Before setting `--assembly_sv_plan_only false`, inspect pair/reference selection and skipped coverage. No Cactus rebuild is needed. Existing experimental alignment tasks are not automatically imported into Nextflow cache; first production pair alignments are new work.

Following gates: one eligible full chromosome-set pair for resource/format validation; then expand to the planned set. Add an explicit bounded pair selector before that first execution if the plan has multiple pairs. Audit final HiFi BAM dictionaries/checksums against final assemblies before read-support evaluation; selected chromosome FASTA hashes alone cannot validate a BAM against the complete reference.
