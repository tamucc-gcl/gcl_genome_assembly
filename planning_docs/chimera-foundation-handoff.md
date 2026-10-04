# Chimera foundation: implementation and cluster checks

This is the first implementation batch from the revised assessment plan. It does
not classify the current fish events as genuine fusions or misassemblies, and does
not change upstream assembly settings.

## Changes ready for testing

- Empty, reference, populated and recovered called-join tables use one ordered
  schema. Missing descriptive fields are explicit dots, not shifted columns.
- Malformed manual instructions fail rather than silently dropping selected rows.
- Recurring composites receive transition evidence even when the scaffold-wide
  vote says NOT_A_CANDIDATE. Scope remains inferred chromosome membership, with
  no new absolute size threshold. Their join verdict is REVIEW.
- Vote-only automatic cutting is explicitly unavailable. Reviewed, checksum-bound
  manual gap cuts retain their existing coordinate checks and sequence-preserving
  behavior. `callable=yes` means a coordinate is eligible, not that an edit is justified.
- Cross-scaffold Hi-C pairs touching assessed scaffolds are lifted and retained in
  `*.alternative_partners.pairs.gz`; within-scaffold profiles are unchanged. These
  contacts are diagnostic and have not yet been scored as independent support.
- Raw hifiasm GFAs, binary caches and BED annotations, and intermediate YaHS AGPs
  are optional published outputs. Existing standardized final graphs/AGPs remain.
  Old scratch-only artifacts cannot be recovered by declaring new outputs.
- The coordinate audit explicitly distinguishes an unavailable reference alignment
  from an assessed alignment. Report prose no longer claims Hi-C/telomere validation
  authorized cuts. Sample-specific causal claims were removed from cut-policy comments.

## Run the regression suite and validation

Run from the cluster project root after syncing this batch. Local Python/R/Nextflow
execution was not performed; the changes need these cluster checks.

```bash
export PYTHONDONTWRITEBYTECODE=1
python3 -m unittest discover -s gcl_genome_assembly/tests -p 'test_*.py'
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv validate \
  --outdir "$PWD/comparisons/chimera-foundation-01"
```

New checks cover mixed empty/populated aggregation, recovered-row normalization,
malformed instructions, recurrence that must not suppress evidence, a fully eligible
vote candidate that must not trigger auto cuts, and reverse-oriented cross-scaffold
liftover. Existing tests still check valid reviewed cuts, checksum mismatches,
unresolved chromosome membership, ambiguous placements and exact older-gap recovery.
Input validation alone does not exercise chimera tasks or establish biological validity.

## Diagnostic run after both checks pass

Use the current two-library input sheets. Keep existing results and work/cache.
This launcher accepts extra Nextflow parameters after its first three arguments.
The comparison output overrides its normal production output path.

```bash
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv assembly \
  --outdir "$PWD/comparisons/chimera-foundation-01" \
  --chimera_break false \
  --run_pangenome false --qc_mode none --run_post_assembly false \
  --run_assembly_sv false \
  --run_final_hifi_mapping false --chimera_hifi_context false
```

The final two switches avoid additional HiFi mapping for this foundation check.
Existing HiFi evidence from the preceding run remains part of the next assessment;
this comparison does not replace it. Verify the launcher prints the intended output
path in the actual Nextflow run (its preliminary echo still shows its default).
Inspect the trace for assembly/scaffolding cache reuse. No upstream command changed
in this batch, but cache reuse must be confirmed from the actual run rather than assumed.
Optional retention outputs do not force recovery of missing files in old cached tasks.

## Collect files for review

### If reference scoring is killed with exit 137 on the current CR_CORE cluster

The observed scheduler reserves cores without accounting for requested memory.
Use the optional comparison override below for a controlled diagnostic: exclusive
nodes, one harmonization task at a time, unchanged eight alignment threads, and a
300 GB memory envelope within the normal nodes' reported capacity. This does not
assert that the task needs 300 GB. It removes other scheduled jobs and supplies
headroom so that actual usage can be measured. Successful execution alone will not
distinguish node contention from a former per-task limit; inspect accounting and
the cluster's kill logs before selecting production resource settings.

```bash
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv assembly \
  -c "$PWD/gcl_genome_assembly/scripts/comparisons/harmonize_isolated.config" \
  --outdir "$PWD/comparisons/chimera-foundation-01" \
  --chimera_break false \
  --run_pangenome false --qc_mode none --run_post_assembly false \
  --run_assembly_sv false \
  --run_final_hifi_mapping false --chimera_hifi_context false
```

This override is for comparison runs only; keep production resource defaults.
An exclusive task can wait until a whole node becomes available. The launcher
already uses resume; do not remove existing results or work. After the run,
replace JOB_ID below and return the accounting output with the evidence archive:

```bash
sacct -u "$USER" --starttime today \
  --format=JobID,JobName%60,State,NodeList,Elapsed,ReqMem,MaxRSS,ExitCode -P \
  > comparisons/chimera-foundation-01/assembly/chimeras/slurm-accounting-JOB_ID.txt
```

Replace JOB_ID with the diagnostic job number. The archive includes small evidence,
coordinate audits, name maps, intermediate/final AGPs, input sheets and logs. Keep
large pairs, FASTAs, graphs, binary caches and BAMs on the cluster.

```bash
python3 gcl_genome_assembly/scripts/comparisons/inventory_join_artifacts.py \
  comparisons/chimera-foundation-01 \
  comparisons/chimera-foundation-01/assembly/chimeras/join_artifact_inventory.tsv
bash gcl_genome_assembly/scripts/comparisons/collect_chimera_review.sh \
  JOB_ID comparisons/chimera-foundation-01 \
  data/assembly_samplesheet.csv data/hic_readsets.csv
```

Return the resulting `comparisons/chimera-review-*.tar.gz` and the unittest/validation
results. The large alternative-partner pairs are intentionally retained separately;
their audit counts are included. Avoid deleting or moving linked work outputs.

## What follows this batch

1. Check changed diagnostic scope, canonical headers, coordinates, pair counts,
   actual cache reuse and available graph/AGP artifacts. Reference assemblies with
   no reference PAF remain unassessed by this transition detector, not cleared.
2. Build the join-origin ledger from retained hifiasm graphs and each YaHS AGP,
   handling corrections by exact anchors rather than assuming stable coordinates.
   Add reference-independent assessment for the reference itself.
3. Identify oriented junction anchors and alternatives; compare library-specific
   Hi-C contacts and competitive, uniquely anchored HiFi support with opportunity
   to span each junction. Separate structural adjacency from phase consistency.
4. Calibrate breakpoint localization and decisions on positive, negative, repeat-rich
   and heterozygous controls. Only then enable evidence-driven automatic editing.
5. Select upstream experiments from the first stage introducing each suspect join:
   version-checked hifiasm post-join/purge alternatives or YaHS filtering/resolution/
   correction alternatives. Keep these as comparison settings until outcomes justify
   a general default change.

Manual cuts still split inside a verified N-gap and preserve its halves. Exact
artificial-gap removal and protection of newly cut ends from telomere extension are
a separate coupled change to implement before using the new automated edit policy.
