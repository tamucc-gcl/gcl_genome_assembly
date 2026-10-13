# General site review policy (v5)

The policy uses measured evidence, not sample names, chromosome numbers, event
IDs, expected answers or genomic coordinates. Manual recommendations remain
advice: generated selections are NO and automatic cutting remains deferred.

- Discover chromosome-scale transitions using independent individuals, keeping
  sister haplotypes as context rather than independent votes.
- Assess each actual gap at MAPQ 20, MAPQ 30 and all primary placements. Unknown
  low-confidence support is not zero. Independently observe both neighboring
  regions and check native contig continuity.
- Recommend an exact manual gap cut only when broad chromosome arms are
  independently confirmed, at least one individual separates the local flanks,
  no qualifying primary molecule spans the gap, both neighbors are observable,
  and suitable Hi-C loss or an observable regional contact trough corroborates
  it. Local chromosome opposition, native continuity or calibrated library
  disagreement prevents this recommendation. Multiple qualified gap options
  require a manual choice. Stricter automatic-cut prerequisites are recorded
  separately and do not trigger application.
- A local chromosome switch in at least two independent individuals, a HiFi
  trough and corroborating contact loss can favor breaking a region without
  selecting a coordinate. Sustained strict-MAPQ losses are separately considered
  while remaining MAPQ20 bridges and competing placements stay counterevidence.
- A regional contact trough is a ratio at most 0.35 of the lower median shoulder
  ratio, using two or more bins per shoulder, 300–700 kb away. Shoulders and
  observable centers need at least 100 within-side pairs. This is a review
  heuristic, not statistical calibration. Low-count centers remain labeled
  limitations and cannot authorize an exact manual cut.
- Preserve separate MAPQ20 troughs instead of merging them into a broad MAPQ30
  loss. Five scattered weak probes do not establish a sustained loss: five
  consecutive 1 kb probes or a recorded weak run of at least 5 kb are required.
- Compatible intervening chromosome order does not establish an exact switch.
  Same-chromosome and conflicting locations remain auditable but are excluded
  from actionable options. Regional assays overlapping a gap share one editable
  gap option with evidence aliases.

## Validation and limits

Shifted-coordinate synthetic tests use unrelated chromosome, scaffold and event
names. Negative cases cover supported gaps, competing placements, unavailable
measurements, same-chromosome context, library disagreement, low-count contacts,
scattered weak probes and multiple qualifying gaps. Production Nextflow wiring
checks cover a single Hi-C library, manual application, reassignment, finishing,
final mapping and cache reuse. These tests establish software behavior; they do
not estimate biological false-positive rates.

A local replay of the ten retained cohort packets checks report presentation.
It cannot fill missing all-primary gap assays or measure newly split locations.
The next real pipeline run must collect those measurements; no expected CTlk
verdict is hardcoded into production logic or treated as a test oracle.

Local validation on 2026-10-12: 329 tests passed. Both production Nextflow
wiring checks passed, including single-library evidence without plots and cached
manual cutting through reassignment, finishing and final mapping. Six candidate
plots were regenerated from the ten-packet replay; the revised layout was visually
checked. Biological decisions still require the refreshed Crest measurements.

## Crest rerun

Run from the project directory after syncing this repository:

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_pipeline.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv \
  comparisons/chimera-site-review \
  --chimera_hifi_bam_manifest "$PWD/comparisons/chimera-regional-review/results/assembly/chimeras/mapping_manifest.json"
```

This runs the normal pipeline with resume and symbolic-link publishing. The BAM
manifest is a reuse shortcut for this rerun, not a requirement for a fresh run.
Outputs include `results/assembly/chimeras/README.md`, per-assembly reports and
`review.tsv`; the launcher creates `comparisons/chimera-site-review.review.tar.gz`.
No cuts are applied by this command.
