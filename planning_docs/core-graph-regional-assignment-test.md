# Regional chromosome assignment correction

The deletion comparison reproduced the cached native PAF exactly and showed that
chromosome-isolated filtering restores mapping evidence for both arms of retained
composite scaffolds. Native regional assignment still discarded one arm.

The pinned upstream `rgfa-split` groups intervals using the aligned-base count where
it should compare reference chromosome IDs. The correction changes that comparison
only. Coverage, uniqueness, smoothing and minimum chunk settings remain unchanged.

## Next run

Sync the repository, then run from the Crest project directory:

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_pipeline.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv \
  comparisons/core-graph-regional-01 \
  --chimera_break "$PWD/data/chimera_reviews/Sde-CTlk_104.v8.reviewed.tsv" \
  --chimera_hifi_bam_manifest "$PWD/comparisons/chimera-observable-contact-v7/results/assembly/chimeras/mapping_manifest.json" \
  --run_pangenome true --pangenome_validate_only false \
  --pangenome_assignment_replay false \
  --pangenome_deletion_filter_replay true \
  --pangenome_regional_patch_replay true
```

The pipeline creates a small compiler environment and builds the pinned original and
corrected tools once. Real CLI tests must demonstrate recovery in both scaffold
orientations, with unchanged ordinary chromosome controls, before the tool is used.
Three corrected regional assignments compare native, chromosome-isolated and disabled
deletion filtering. Expensive graph construction, filtering and existing assignment
tasks retain their original inputs and process bodies for ordinary Nextflow resume.

The new flag requires deletion-filter replay and cannot be combined with validation-only
mode. This remains a diagnostic; it does not rebuild the production graph.

## Acceptance evidence

Read `results/pangenome/373251/deletion_filter/regional_assignment_report.md`.
Return `comparisons/core-graph-regional-01.review.tar.gz`.

- Compare both chromosome arms with the input filtered mappings and native assignments.
- Inspect every ordinary chromosome control decrease, rather than hiding regressions.
- Every assigned or ambiguous region must be in bounds and non-overlapping within its
  original assembly/scaffold. Adjacent regions are allowed; duplicate regions fail.
- Each output PAF record must belong to its assigned chromosome fragment and use valid
  fragment-relative coordinates. Source coverage uses coordinate unions.
- Ambiguous and unreported sequence remains explicit. It is never copied into another
  chromosome to manufacture representation.
- Reports retain regions, surviving mapping intervals, assignment logs, build checksums
  and the exact one-line patch. Intermediate PAF/GFA/binaries stay in the work cache.

After accepting these results, integrate the supported filtering/assignment policy into
production construction and audit FULL/CLIP source representation and distinct-node
accounting. Passing intermediate coordinate checks alone does not certify the final graph.

## Local validation

57 focused Python tests pass, including 17 new source-ownership/provenance/report tests.
The real original and corrected upstream CLIs reproduce the grouping defect in both
orientations: affected-arm assignment rises from 2.0 Mb to 3.4 Mb; the ordinary
chromosome control remains 4.5 Mb. Nextflow 23.10.1 integration passes with the real
corrected assignment executable and small deletion-filter doubles. It retains both
4 Mb arms under isolated/disabled filtering, retains only the surviving 4 Mb under
native filtering, passes coordinate audits, and confirms preparation/filtering/old
assignments resume from cache when the correction flag is enabled. These are local
software checks; real cohort acceptance still requires the Crest run above.
