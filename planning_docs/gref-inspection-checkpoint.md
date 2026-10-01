# Bounded GREF inspection

The follow-up now writes standard/raw.biological_slots.tsv. Identity graph suffix .1/.2 defines the expected first/second GT positions; an unsuffixed reference expects its sole first slot. Missing unused positions are counted separately from biological missing observations. Called unused positions and unphased partial-slot assignments are flagged rather than repaired or silently discarded. Short GT fields that cannot represent the expected position remain unresolved. This checks the export convention; it does not establish phase accuracy or independent population sampling. Original VCFs and AC/AN checks are unchanged. Rows without GT are reported as no_GT and cannot contribute slot counts. Interpret slot totals only alongside inspection.json flags.

The next comparison commands and acceptance gates are in [sv-benchmark-checkpoint.md](sv-benchmark-checkpoint.md).

The catalog audit passed: 32,176,170 records, six VCF columns correctly mapped to five individuals/ten input haplotypes. The next job inspects exported genotype and allele metadata before choosing decomposition or structural-variant benchmark commands. It does not rerun Nextflow or alter graph/VCF files.

The collector chooses up to twelve populated coordinate sequences in each naming stratum (reference-named, _alt, other), spread across sequence-length ranks. It takes 10 kb windows at start/middle/end and merges overlaps. These are deterministic diagnostics, not a random or representative sample and not a mechanism-based SV sensitivity benchmark. A VCF record may itself span more than the window size; exact allele strings are retained in compressed extracts. Raw GREF is inspected only when an existing raw index is available; an expensive new index is not created.

Outputs include compressed standard/raw window VCFs, full headers, genotype-slot/missingness distributions by original VCF column, example allele lengths and nesting/coordinate tags, AC/AN consistency against exported GT, AT cardinality, the actual segment-map schema prefix, ledger and reproducible command manifest. GT slots are not assumed to prove biological ploidy. AC/AN disagreements are diagnostic flags, not grounds for silently fixing genotypes or filtering records. Equal AT cardinality does not prove traversal correspondence after decomposition.

All scripts stay in scripts/comparisons and outputs in comparisons. Static code review only locally; tests and analysis are for the cluster.

## Run after syncing

From the cluster project root:

```bash
export PYTHONDONTWRITEBYTECODE=1
python3 -m unittest discover -s gcl_genome_assembly/tests -p 'test_*.py'
sbatch gcl_genome_assembly/scripts/comparisons/run_gref_inspection.sbatch \
  --pangenome-dir genome_assembly/pangenome/373251 --taxid 373251 \
  --conda-cache /work/birdlab/.conda_builds
```

The optional cache argument locates an existing bcftools 1.21 executable if bcftools is not already available. No environment or package is installed. To choose an exact executable instead, pass --bcftools /absolute/path/to/bin/bcftools. The executed path and tool version are recorded.

The job prints the comparisons/gref-inspection-*.tar.gz path to return, together with logs/gref-inspection-JOBID.out. No full pipeline rerun or baseline move is needed. Retain the current work/results.

## Follow-on method selection

Inspect these examples before building the known-answer INV/DUP/translocation cases. Keep standard catalog statistics distinct from fine allele decomposition and assembly rearrangements. Tool tests must explicitly account for reference versus graph-derived coordinate space, fragmentation, nested variants and genotype representation. The method shortlist and acceptance requirements remain in pangenome-variant-checkpoint.md and assembly-pangenome-agreed-plan.md; this job does not declare a winning classifier or automatically launch a full-cohort benchmark.
