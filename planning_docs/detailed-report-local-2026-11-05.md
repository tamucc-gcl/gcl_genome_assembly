# Detailed reporting batch — local acceptance and Crest test

## Completed presentation changes

- Final combined chimera review figures and their close-ups now render with R,
  ggplot2 and patchwork. Both the combined and placement-confidence figures remain
  in the detailed report. Peer intervals, independent-consensus shading, HiFi
  confidence/depth tracks, separate primary-placement scales, per-library contact
  profiles, telomeres and assay/gap/submitted-cut markers are preserved.
- The close-up window matches the accepted Python selection function, including
  user-selected positions. Assigned intervals and every profile row reconcile to
  the retained evidence. Missing tracks are labelled unavailable.
- Junction contact maps at 10 and 25 kb now render in R from a lossless sparse
  export of the retained matrices. Dimensions, raw totals, zero/nonzero cell counts
  and matrix coordinates are checked. The two resolutions are unchanged; 5 kb
  matrices remain retained evidence rather than extra figures.
- Contact maps have per-library titles, exact site/coordinate keys, distinct assay
  and submitted-cut markers, readable logarithmic count ticks and a common colour
  scale across the event's libraries at each resolution. This is raw-count
  comparison, not normalization for sequencing depth or a new cut rule.
- PNG and PDF outputs, R/package versions, render row counts, window/marker/axis
  specifications and source tables are retained. The plotting dependencies are
  pinned; there are no task-time package installations.
- Pangenome sharing bars now display their measured percentage and Mb values.
  Main detailed sharing/regional figures are visible without expanding collapsed
  sections. The grey-cell/composite explanation precedes the regional figures.
  PDF downloads and source/audit links remain available.
- The per-species pangenome report now includes the measured individual-sharing
  table before its biological and technical report links. Thresholds, composite
  selection, graph evidence and source-review flags are unchanged.
- Published-report validation follows linked Markdown detail reports recursively,
  checks their images/downloads, and terminates on cyclic navigation. Existing
  directory links are allowed; images and registered artifacts must be files.
- Review archives now include compact GenomeScope figures/estimates. The prior
  archive omitted those figures even though the detailed report linked them.

## Rendering scope and limits

The final review overview and junction maps are rendered in R. Existing whole-
scaffold pooled contact images remain intact because their matrices were not
retained in the evidence packets. Earlier assessment/enrichment tasks still make
their Python plot versions; the final review replaces those particular overview
and junction images with the accepted R versions. Their calculations and task
scripts were deliberately left unchanged to preserve evidence cache reuse.
This is not a claim that every pipeline figure is R-only.

No assembly mapping, candidate discovery, assay, ranking, biological threshold,
manual decision or graph-construction rule changed. Automatic cutting remains
deferred. The new sparse exports retain measured zeros through explicit matrix
dimensions; they do not omit low-count evidence.

## Local verification

- End-to-end actual review CLI: six CTlk events, four verified cuts and two
  retentions, with both figure types and separate contact maps.
- Six actual events: all retained interval/profile rows and accepted close-up
  windows reconcile; the R renderer does not modify evidence.json.
- 157 local links/images in the detailed chimera preview resolve.
- Generic missing-assay, ambiguous-depth-spike, user-cut window, reversed focal
  chromosome orientation and single-library matrix round-trip checks pass.
- Main report renderer: eleven general cases plus the actual ten-assembly run pass,
  including disabled, missing outputs, first pass, reviewed and manually added cuts.
- Nextflow 23.10.1 cohort/review wiring, multiple-species report wiring and main
  input/compile checks pass. The cohort fixture uses external-tool doubles and
  plots=false; actual plotting was tested separately with the real R environment.
- Nine report-validation tests, six existing manual-review tests, two placement
  adapter tests, a matrix-export test and twelve sharing/context tests pass.
- Sharing/regional figures were rendered using R 4.3.3, tidyverse 2.0.0, ggplot2
  3.5.2 and patchwork 1.3.0 and visually inspected. Revised contact/combined figures
  were also inspected. PDF font warnings were resolved using Cairo exports.

Local previews (not production inputs):

- `comparisons/detailed-report-local-02/assemblies/`
- `comparisons/detailed-pangenome-local-01/`
- `comparisons/main-report-validation/accepted_run/`

## One normal cached Crest integration run

Sync the complete repository changes, including new helpers, R scripts and the
review environment. From the project directory on Crest:

```bash
cd /work/birdlab/GCL/spratelloides_delicatulus_genome
mkdir -p logs

CORE_INTEGRATION_ALLOW_REBUILD=false \
sbatch gcl_genome_assembly/scripts/comparisons/run_core_integration.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv \
  comparisons/main-report-01/results \
  data/chimera_reviews/Sde-CTlk_104.v8.reviewed.tsv \
  comparisons/detailed-report-01 \
  --run_report true --report_chromosome_painting true \
  --chimera_hifi_bam_manifest "$PWD/comparisons/chimera-observable-contact-v7/results/assembly/chimeras/mapping_manifest.json"
```

The BAM manifest is the existing test-reuse override, not a fresh-run requirement.
This uses normal Nextflow resume and symlink publication. Assemblies, scaffold
alignment/read mapping, chromosome painting and core graph producers should stay
cached. Final review, sharing/regional plots, reports and affected downstream
provenance can execute; a new small plotting environment is expected.

Expected native acceptance: Nextflow success, report_validation PASS and
core_integration_check PASS; unchanged final names/reference/source-audit evidence
and cached Cactus construction. Do not enable allow-rebuild just to bypass a
failure. The audit detects unexpected rebuilds after execution; it is not a
pre-execution scheduler veto.

Review `results/assembly_report.md`, the linked detailed chimera reports and the
pangenome pages. Return `comparisons/detailed-report-01.review.tar.gz`. The local
review archive excludes FASTAs and graph/index products, so full publication-link
and core identity acceptance must be checked on Crest. Category 6 remains pending
that native acceptance and presentation review; this batch does not reopen biology.
