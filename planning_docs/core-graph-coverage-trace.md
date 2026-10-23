# Cached graph coverage trace — 2026-10-23

Core-graph-v8 completed successfully (838 cached / 23 completed tasks). The separated
CTlk chr4 and chr14 paths now survive in CLIP. Core graph acceptance remains open:
the retained chr7/chr12 composites have 75.39/83.43 Mb input lengths but only
40.64/40.81 Mb retained CLIP path lengths. Existing regional and relaxed assignment
tests do not recover both portions. The missing coverage is largely already absent
in the filtered mappings, before assignment; raw mapping versus filtering is the
next question. This does not change the reviewed assembly cuts or retained joins.

## Implementation

The existing optional `--pangenome_assignment_replay true` branch now consumes the
cached raw GAF as well as the raw PAF and SV graph. It selects all missing/FULL-only
sources, composites, chromosome-scale sources retaining less than 80% of input path
length, and two ordinary represented chromosome controls per assembly. Selection
uses source tables and identities, with no sample names, chromosome pairs or cut
coordinates embedded in the implementation.

The 80% value is a general screening threshold for investigation, not an acceptable
loss rule, automatic graph-policy change or biological deletion criterion. The
always-on source audit records this threshold and flags partially retained
chromosome-scale paths as REVIEW while preserving their name-presence status.

Query-coordinate coverage is unioned separately for raw GAF, raw PAF, pinned Cactus
filtered PAF and each of the four assignment cases. Exact losses and gains are
coordinate set differences, with each assignment compared to the common filtered
input. Merged intervals and 1 Mb bins retain original source coordinates. Assignment
region labels provide context, including ambiguous regions. CLIP traversal lengths
are shown separately; neither those lengths nor FULL W-coordinate spans are treated
as actual aligned query-coordinate coverage.

## Run after syncing

From the Crest project directory, with no other controller using the same resume
session:

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_pipeline.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv \
  comparisons/core-graph-coverage-01 \
  --chimera_break "$PWD/data/chimera_reviews/Sde-CTlk_104.v8.reviewed.tsv" \
  --chimera_hifi_bam_manifest "$PWD/comparisons/chimera-observable-contact-v7/results/assembly/chimeras/mapping_manifest.json" \
  --run_pangenome true --pangenome_validate_only false \
  --pangenome_assignment_replay true
```

Ordinary Nextflow resume should reuse Cactus construction and assembly work. The
source audit, diagnostic preparation, four parallel replay cases and summary rerun;
their report consumers may rerun too. Construction process code and its inputs were
not changed. The HiFi manifest is reuse of earlier chimera mappings for this project
resume, not a new graph-trace requirement or a normal fresh-run requirement.

## Outputs and next decision

Start with `results/pangenome/373251/assignment/assignment_report.md`. Its main tables
now include retained composites and low-retention chromosome-scale paths, followed
by a raw/filtered/assigned coverage comparison. The full underlying products are:

- `coverage_trace.tsv`: every selected source and stage, covered fraction, and lost/gained bases.
- `coverage_bins.tsv`: 1 Mb source-coordinate coverage and assignment-region context.
- `coverage_changes.tsv.gz`: exact lost/gained query intervals.
- `coverage_intervals.tsv.gz`: merged covered query intervals at every stage.
- Existing assignment comparison, decisions and region tables remain available.

Exact interval tables are compressed and included in the review archive even if
they exceed the ordinary compact-file size cap. Multi-GB mapping intermediates and
the prepared graph stay in Nextflow work directories. No new graph is constructed
by this diagnostic branch.

If coverage is absent in raw GAF, inspect graph mapping. If it survives GAF but not
raw PAF, inspect conversion. If it survives raw PAF but not filtering, inspect the
pinned filter's rejection conditions. If it survives filtering but not assignment,
inspect chromosome partitioning. Compare query intervals as well as total coverage.
Select the corresponding correction only after this batch establishes the cause.

Validation: focused Python tests cover partial-retention alarms, general selection,
union coverage, losses versus gains, coordinate offsets, ambiguous assignment and
retained-composite report inclusion. Production Nextflow 23.10.1 modules pass the
tiny diagnostic integration, including all four cases and compressed evidence
publication, plus the two-species construction-to-diagnostic input fixture. External
tools in these fixtures are doubles; the pinned biological run on Crest remains
necessary to establish the cause.
