# Manual misassembly review: completed improvement pass

Policy: cohort-manual-review-v6. No sample, event ID or coordinate determines a verdict.
Automatic application remains deferred; generated cut selections are NO.

## Implemented

- R assay coordinates have purple dashed lines across all five evidence panels;
  shaded ranges retain uncertainty. Red dashed lines identify nonexcluded exact
  gap options. Submitted user-selected cuts are green solid lines.
- Each event has an overview and a close-up, centered on submitted cuts, the
  recommended gap or priority region. Figure labels are set explicitly and
  exported with tight bounds. Broad original brackets have faint shading.
- Summary tables show MAPQ20/MAPQ30/all-primary counts and per-library contact
  findings. Probe counts are actual reads at the assay; original-interval minima
  are labeled separately. Low-count troughs are limitations, never positive
  contact evidence. Hap2 C02 remains a possible misjoin under review, rather
  than receiving stronger advice from unobservable contact loss.
- Every event has a blank exact-coordinate row, even if gap options exist.
  Human-selected coordinates can lie outside the suggested review range.
  An undetected cut can be added as a separate provenance-bound row.
- Reports distinguish submitted cuts from detector advice, list exact submitted
  coordinates and include sequence-reconstruction verification counts.
- Shared figure rendering lives in misassembly_plot.py. Assessment and report
  processes stage it; the review process declares its plotting environment.
  Packets without figures remain without figures.

## Verification

336 tests passed. Production Nextflow modules passed single-library evidence
and mixed manual-cut tests with resume. Mixed input applied a literal gap cut,
an internal cut outside the detected interval, and an added undetected cut;
four resulting pieces passed reassignment, finishing, final naming and final
read mapping. External mapping/finishing tools in these integration tests are
software doubles; this validates orchestration, not biological alignment accuracy.

All ten retained real packets were reclassified and rendered locally. Six
original-event rows and all twenty review rows remained unselected; every report
link resolved. Figures were visually inspected, including the C01 and C03
close-ups. No real assembly was cut during local verification.

## Crest validation

After syncing, run from the project directory:

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_pipeline.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv \
  comparisons/chimera-manual-review-v6 \
  --chimera_hifi_bam_manifest "$PWD/comparisons/chimera-site-review/results/assembly/chimeras/mapping_manifest.json"
```

This is the full production pipeline with resume, symlink publishing and retained
mapping reuse. It produces `comparisons/chimera-manual-review-v6.review.tar.gz`.
It applies no cuts.

Copy the generated review.tsv outside generated output before editing. For an
internal cut set selected=YES, review_disposition=CUT, action=BREAK_PROBABLE_MISJOIN,
localization_status=localized, an exact cut_bp, reviewer and reason; leave gap
fields blank. Preserve the source checksum and original scaffold. To apply the
reviewed file, add `--chimera_break "$PWD/data/chimera-reviewed.tsv"` to the same
launcher with a separate output directory. Changed assemblies are reassigned,
finished and remapped downstream as required. Do not select an event row and a
gap row for the same physical cut.
