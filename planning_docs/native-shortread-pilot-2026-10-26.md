# Native short-read pilot — 2026-10-26

## Scope

Run the normal `main.nf` assembly route for the user-supplied CMat_061 shotgun
mates. Species is Spratelloides delicatulus, taxid 373251, as in this project.
This is an isolated native execution/routing pilot, not assembly-quality acceptance.
The main HiFi/pangenome run and its cache are untouched.

Use full original files without a prefix subset. Normal trimming, genome-size
estimation and the configured 60x short-read depth cap remain active. Organelle
recovery/annotation use real tools and databases. The test config reduces SPAdes
to k=21,33, disables its read-correction stage (`--only-assembler`), disables isolate
mode, sets coverage cutoff off, and uses one Redundans iteration. Reduction,
scaffolding and gap closing remain enabled. These settings reduce work but do not
guarantee a short runtime on this approximately gigabase genome. QC, graph and
post-assembly analyses are off for this bounded route check.

## Crest invocation

Sync the repository first. From the project directory containing
`gcl_genome_assembly/` and `raw_fastq/`:

```bash
mkdir -p comparisons logs
python3 gcl_genome_assembly/scripts/comparisons/prepare_shortread_smoke.py \
  raw_fastq/ssl/Sde-CMat_061_1.fq.gz \
  raw_fastq/ssl/Sde-CMat_061_2.fq.gz \
  comparisons/core-shortread-native-01 \
  --full-input --sample Sde-CMat_061 --taxid 373251

sbatch gcl_genome_assembly/scripts/comparisons/run_shortread_smoke.sbatch \
  comparisons/core-shortread-native-01 assembly
```

The preparation command rejects overwriting an existing prepared directory. For
an ordinary retry, rerun only the `sbatch` command; it uses Nextflow resume in the
pilot's own launch/work directories. The sheet records absolute original FASTQ
paths and the species taxid; no read copies are made. Its `subset.json` explicitly
records full-input mode despite the existing manifest filename.

Return `comparisons/core-shortread-native-01.assembly.review.tar.gz`, generated
on success or failure after launch. It includes exit status, trace, Nextflow log,
preparation manifest, sample sheet, test config and compact result reports/logs,
annotations and figures. Large reads, FASTAs and graphs remain on Crest.

## Acceptance checks

1. Input resolves to one collapsed primary assembly with biological ploidy 2;
   taxid and genetic code propagate correctly.
2. Actual GetOrganelle database/setup/recovery and any required annotation complete,
   or verified no-recovery explicitly releases nuclear processing. A tool failure
   must fail the run, not become a successful no-recovery result.
3. Native trimming/subsampling/SPAdes and enabled Redundans stages execute with
   expected resources; final assembly and run-summary outputs exist.
4. No sample disappears from required context joins; no HiFi/Hi-C/graph work is
   invented for this shotgun-only sample.
5. Review native warnings, task durations and any unexpected omissions before
   claiming either input or organelle category complete. Controlled taxonomy/name
   lookup and retry/dependency checks remain separate outstanding work.

Local preparation tests: three passed. Modified launcher: Bash syntax check passed.
Native execution is pending Crest evidence.
