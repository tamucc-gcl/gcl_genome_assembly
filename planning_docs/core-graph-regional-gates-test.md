# Regional assignment: observe rejection gates

`core-graph-regional-01` passed execution and ownership audits, but the grouping
correction changed none of the cohort's assigned chromosomes or surviving mapped
coverage. This diagnostic measures why the existing rules reject regions before
considering additional source or parameter changes.

## Scope

`--pangenome_regional_gate_trace true` adds a removable observation branch to the
existing regional-patch diagnostic. It requires both deletion-filter replay and
regional-patch replay. Existing compiler, graph, filtering and assignment tasks
retain their original inputs and process bodies for ordinary Nextflow resume.

One new compiler task builds an unmodified original tool and two observed tools:
original grouping with logging, and corrected grouping with logging. The original
vendor sources remain unchanged. Observation adds statistics/log writes without
changing existing decision expressions or thresholds. Applied observation diffs,
source checksums, binary checksums and compiler commands are retained.

Three independent jobs inspect native, chromosome-isolated and disabled deletion
filtering. Each runs original and corrected grouping, with and without observation.
Every `.paf` and `.fa_contigs` output file must be byte-identical between observed
and unobserved runs of the same variant. Any difference fails the batch. Output
region ownership and mapping-to-fragment checks must also pass.

The logging records:

- Merged interval counts, longest interval, covered span sum and aligned-base counts
  per chromosome before and after overlap resolution.
- Each secondary candidate's original source coordinates, interval count, longest
  interval, configured minimum, aligned-base density and the actual Boolean outcome
  of both local gates. End extensions are included as they are in the tool.
- Dominant regions removed by the small-region filter.
- Whole-query adjusted coverage, threshold, masked bases, softclipped bases and
  accepted secondary span; local acceptance can be overruled here.
- Final assigned regions and coordinate audits, including explicit ambiguous parts.

The tools must first pass actual native-CLI known-answer cases for the strict 1 Mb
boundary, just above that boundary, low density, overlap removal and whole-query
coverage rejection. These measure the upstream code rather than reproducing its
algorithm in Python. No configuration sweep or production graph rebuild is included.

## Run after syncing

From the Crest project directory:

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_pipeline.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv \
  comparisons/core-graph-gates-01 \
  --chimera_break "$PWD/data/chimera_reviews/Sde-CTlk_104.v8.reviewed.tsv" \
  --chimera_hifi_bam_manifest "$PWD/comparisons/chimera-observable-contact-v7/results/assembly/chimeras/mapping_manifest.json" \
  --run_pangenome true --pangenome_validate_only false \
  --pangenome_assignment_replay false \
  --pangenome_deletion_filter_replay true \
  --pangenome_regional_patch_replay true \
  --pangenome_regional_gate_trace true
```

Read `results/pangenome/373251/deletion_filter/regional_gate_report.md` and return
`comparisons/core-graph-gates-01.review.tar.gz`. The report links the complete
chromosome summaries and byte-equivalence checks. Each `trace_<policy>` folder
contains compressed, per-candidate events for both grouping variants and final
regions/coordinate audits. Large mappings, compiled executables, source copies and
raw instrumented logs remain in the ordinary work cache. The review archive includes
the compressed event tables, compact compiler provenance and observation diffs.

Interpret the evidence before changing policy: a substantial run rejected because
its individual intervals are short is different from low density, conflicting
overlaps, a whole-query coverage failure or missing input mappings. Do not infer
chromosome presence from scaffold names or relax thresholds to force expected results.

## Local acceptance

79 focused Python tests pass. Five native-CLI gate regressions pass with byte-identical
observed/unobserved outputs. The Nextflow 23.10.1 integration and resume check passes:
all six policy/grouping equivalence checks pass, final coordinate checks pass, compact
reports/events/compiler provenance are published, and the previous preparation,
filtering, compiler and assignment tasks remain cached when observation is enabled.
External deletion filtering uses small tool doubles in that wiring check; regional
assignment and gate tests use the actual pinned C++ source. Real-data acceptance
requires the Crest packet above.
