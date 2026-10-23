# Supported regional chunk-setting comparison

The gate diagnostic established that restored chromosome arms pass density and final
coverage gates but fail the strict >1 Mb individual-interval requirement. This batch
tests existing CLI settings before any additional source-policy changes.

## Matrix and cache

Compare `-u 1000000`, `-u 100000` and `-u 25000`, with original and corrected grouping,
under native, chromosome-isolated and disabled deletion filtering: 18 cases total.
The six verified 1 Mb cases are read directly from the existing gate-run outputs.
Twelve new independent assignment jobs run the smaller settings, followed by one
summary task. Existing graph construction, filtering, compiler builds, gate runs and
other assignments retain their task bodies/inputs for normal Nextflow resume.

The only assignment option changed by each new case is `-u`. All other alignment,
quality, coverage, uniqueness and softclip settings remain the same. The executables
are the verified unmodified original and corrected binaries already built in cache;
the logging executables are not needed for these new assignments. No production
default, input assembly or reviewed cut changes.

Every case covers the complete existing diagnostic selection, including missing and
low-retention sources, composites, fragmented-assembly sources and ordinary chromosome
controls. This selection is not every scaffold in every assembly. Its membership is
published in the case folders and must match exactly across all 18 cases.

## Evidence and acceptance

- Mapped source-coordinate union coverage and changes from the matching 1 Mb baseline.
- Newly assigned or lost chromosomes and changes to ambiguous sequence.
- Explicit source-bound, non-overlap and mapping-to-fragment audits.
- Chromosome-level added/removed spans and exact source intervals, including fragment
  counts, largest assigned region and mapped fraction for changed chromosomes.
- Review flags for new chromosomes on ordinary controls, newly assigned chromosomes
  represented only by regions below 1 Mb, and removed regions. Flags identify side
  effects for inspection rather than proving biological errors.

Do not select a setting solely because the motivating arms recover. Assess ordinary
control losses, new small secondary assignments and fragmented-assembly behavior.
Both grouping versions remain in the comparison so the need for a source correction
can be assessed separately from the chunk parameter.

## Run after syncing

From the Crest project directory:

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_pipeline.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv \
  comparisons/core-graph-chunks-01 \
  --chimera_break "$PWD/data/chimera_reviews/Sde-CTlk_104.v8.reviewed.tsv" \
  --chimera_hifi_bam_manifest "$PWD/comparisons/chimera-observable-contact-v7/results/assembly/chimeras/mapping_manifest.json" \
  --run_pangenome true --pangenome_validate_only false \
  --pangenome_assignment_replay false \
  --pangenome_deletion_filter_replay true \
  --pangenome_regional_patch_replay true \
  --pangenome_regional_gate_trace true \
  --pangenome_regional_chunk_replay true
```

Read `results/pangenome/373251/deletion_filter/regional_chunk_report.md`.
Return `comparisons/core-graph-chunks-01.review.tar.gz`.
The report links complete per-source comparisons, chromosome changes, review flags
and compressed exact added/removed intervals. New `chunk_<policy>_<variant>_<size>`
folders contain assignment coordinates, mapping summaries, coordinate audits, command
provenance and native assignment logs. Intermediate PAFs are audited and removed.

## Local validation

91 focused Python tests pass. Six real native-CLI cases demonstrate the expected
behavior for a long secondary run made of 40 kb intervals: neither grouping version
assigns it at 1 Mb or 100 kb; at 25 kb the original tool assigns 60 fragments and the
corrected tool assigns one coherent region. Every output passes source-coordinate
ownership checks. Nextflow 23.10.1 integration verifies all 18 matrix cases, exactly
12 new assignment jobs, published report tables, and cached filtering, compiler,
gate and earlier assignment tasks. External deletion filtering uses small doubles
in this integration check; the regional assignments use the actual pinned binaries.
Real cohort acceptance still requires the Crest run above.
