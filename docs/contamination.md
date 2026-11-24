# FCS contamination screening and optional BlobTools evidence

## Standard production route

`--run_decon_contigs true` (default) enables both FCS-adaptor and FCS-GX. They
share the FCS stage control; there is no separate adaptor enable switch.
`--run_decon_contigs false` disables this standard route, including the final
screen. The older optional intermediate scaffold screen remains available through
`--run_decon_scaffolds true` and also includes both tools.

The contig-stage order is:

1. FCS-adaptor screens the stage input with explicit `--euk` or `--prok` mode.
2. The FCS GX `clean-genome` executable applies the adaptor report. This handles
   terminal trimming, internal splitting and sequence exclusion, including internal
   adaptors not removed by FCS-adaptor's own convenience FASTA.
3. FCS-GX screens the adaptor-cleaned FASTA against the configured GX database.
4. GX cleaning applies its contamination action report as before.
5. GC/read-depth plots use every contig in the exact GX input, including sequences
   subsequently removed by GX. This keeps classification and coordinates consistent.

The adaptor image lacks a shell. Nextflow runs a **host Python task**, which invokes
Singularity/Apptainer (or Docker) with `/app/fcs/bin/av_screen_x` directly. Nextflow's
own container wrapping is disabled for this task, so it never asks the image to
run `/bin/bash`. The host runtime must be available on PATH; its native image cache
is used. A local SIF path can also be supplied through `fcs_adaptor_container`.
The launcher does not guess a cleaned output filename or use shared temporary
filenames. `decon_container_engine` selects the host runtime (`singularity` by
default; `apptainer` and `docker` are supported). `decon_fcsadaptor_mode` is
`euk` by default; use `prok` for prokaryotes. The adaptor image defaults to
`docker://ncbi/fcs-adaptor:0.5.5`; `fcs_gx_container` controls the GX image. Set an
audited local SIF or immutable image reference as needed for the deployment.

## Final screen placement

The second FCS pass is on the **exact finalized FASTA**, after gap filling,
telomere extension and final naming/reorientation. It finishes before finalized
assemblies are released to final QC and downstream consumers. It runs both
adaptor and GX screens. `--run_decon_final false` disables this additional audit
while retaining the contig FCS route.

The final pass is an audit: it does not automatically trim, split, mask or rename
finalized sequences. That would invalidate coordinate maps, alignments, indexes
and naming evidence. Final adaptor or GX cleaning recommendations are displayed
as **REVIEW / CLEANUP REQUIRED**. A successful task exit is not a clean bill of
health. Resolve findings in a tracked matching coordinate frame, regenerate the
necessary downstream products, and screen the resulting final assembly again.

The final plot reuses final QC coverage when QC is enabled. If QC is disabled,
the FCS diagnostics branch computes the coverage needed for its plot itself.
Contig FCS diagnostics map reads to the pre-GX-cleanup frame; post-FCS BlobTools
read evidence, when enabled, uses a separate matching post-cleanup frame.

## Reports and publication

FCS outputs stay together, with assembly identifiers in their filenames:

- `contig/decontam/fcs/`: adaptor reports, cleaning logs, GX reports and FCS plot evidence.
- `contig/decontam/contamination_report.md`: all contig assemblies in one report.
- `final/decontam/fcs/` and `final/decontam/contamination_report.md`: final audit results.
- `scaffold/decontam/`: corresponding outputs for the optional intermediate screen.
- `contig/decontam/blobtools/<assembly>/`: separately enabled full-search results.

The main assembly report links to the contamination reports. Each combines FCS
and optional BlobTools evidence only when those analyses ran. Disabled analyses
are explicitly identified. Adaptor-input, GX-input and post-FCS coordinates remain
separate. FCS action counts are recommendations, not measurements of bases removed.

The combined GC/coverage figure has one panel per assembly, shared axes and
colours, and one proportional pie per complete record. Pie area represents
record length and wedges represent category base-pair proportions. Small centre
marks keep tiny records visible. It uses **FCS-GX taxonomy results**, not
protein/nucleotide search hits. The seven display categories are primary division,
repeat with expected division reported, low GX coverage with expected division
reported, repeat/low GX coverage with another division reported, other inconclusive
results, taxonomy contaminant/chimeric results, and unreported span. Expected
division is the `asserted-div` in GX metadata; when unavailable, no expected-division
assignment is invented. Repeat/low-coverage results remain inconclusive. Missing
measurements and zero-depth contigs are counted explicitly. No GC, coverage,
contig-length or coding-content threshold selects which contigs receive a screen.

The `report/` directory retains the PNG, PDF, plotting scope, FCS plot evidence,
assembly summaries and separate GX/adaptor action tables. Per-assembly
`*.fcs_metrics.json` files bind plot metrics to FASTA, taxonomy and coverage
checksums. Plotting uses pinned R/ggplot2/scales dependencies; it does not change
screening or cleaning calculations.

The same palette is shared by the GC/coverage pies, mixed-record composition
bars and assembly-specific detail pages. The overview at
`overview/decontam/contamination_report.md` combines current-run stage tables:
one assembly per row, contigs on the left and finalized assemblies on the right.
Bars partition the full span of records previously labeled mixed, including
explicit unreported intervals. They do not use aligned bp as their denominator.
Every stage report links to assembly detail pages with GC/coverage pies, the 30
longest mixed records (ties by identifier), complete tables, and independent GX
alignment/masking summaries. Masking tracks may overlap and must not be summed.
Coordinates, raw results and all records remain in the tables regardless of the
display cap. Non-overlapping segments and complete base-pair partitions are
checked before plotting. No new biological threshold or removal rule is added.

## Optional full BlobTools workflow

`--run_contamination_evidence false` remains the default. Enabling it selects
**full** mode. The heuristic staged mode and its sampling parameters have been
retired from the production interface; old staged runners stop with an explanatory
error rather than silently launching a much larger full search.

The full branch retains protein BLASTX, optional assembly-wide MEGABLAST, optional
sensitive BLASTN follow-up, local read-support plots, a BlobToolKit dataset, and
explicit reviewed cleanup. Original 100 kb windows, 5 Mb chunks, search settings,
node-worker queue and durable nucleotide chunk checksums are retained. No new
large search is needed for standard FCS diagnostics.

The full branch is implemented, with substantial component and earlier pilot
validation. A complete single-assembly acceptance run against the full production
protein and nucleotide databases remains deferred. Do not describe that deployment
as fully accepted based on the retired sampled cohort trial. Its substantial
runtime and resources must be considered before enabling it.

Add these options to a normal assembly run only when requesting full searches:

```bash
--run_decon_contigs true \
--run_contamination_evidence true \
--contamination_search_mode full \
--contamination_protein_db_dir "$protein_bundle" \
--contamination_nucleotide_db "$nt_prefix" \
--contamination_megablast_execution nodes
```

Use `--stop_after_contig_review true` to stop before scaffolding. Sensitive BLASTN
remains separately optional through `--contamination_targeted_blastn true`.
Reviewed BlobTools decisions can be applied on a resumed run through
`--contamination_review "$reviewed_tsv"`; preserve its source checksums, coordinate
stage and review IDs. This applies decisions in the post-FCS frame and does not
reinterpret final-audit coordinates as post-FCS coordinates.

Managed protein and nucleotide database downloads remain available. Their
previous command examples and benchmark provenance are retained in the
[historical trials document](contamination-legacy-trials.md); its sampled-run
commands are obsolete.

## Bounded FCS-adaptor acceptance on Crest

The local tests exercise shell-less application launch with native Singularity
wrapping enabled globally, scratch-staged helper imports, workflow wiring, final-audit byte preservation, taxonomy
segment handling, exact-frame checks, missing/zero-depth handling and real R
figures. Screen binaries in the orchestration test are stubs, not a substitute for
testing the actual containers on Crest.

Submit this small real-container test from the project root:

```bash
mkdir -p comparisons/fcs-adaptor-acceptance-01
sbatch --parsable gcl_genome_assembly/scripts/comparisons/run_fcs_adaptor_acceptance.sh
```

If an image is already downloaded, the acceptance launcher accepts
`sbatch --export=ALL,FCS_ADAPTOR_IMAGE=/absolute/path/to/fcs-adaptor.sif ...`.
It passes that local image to the host runtime directly, preserving the fixture
download cache and avoiding another image pull. The per-assembly
`*.fcs_adaptor_launch.json` records the image, runtime, exact application arguments
and source FASTA checksum.

It downloads NCBI's published approximately 12 Mb test FASTA, checks the 11 known
adaptor calls, verifies five terminal trims, five internal splits and one complete
sequence exclusion against exact expected bases after GX's uppercase conversion, then repeats adaptor
screening and requires zero remaining calls. It uses the actual adaptor and GX
cleaning containers, without a large GX/nt/protein database download or production
assembly run. Return `comparisons/fcs-adaptor-acceptance-01.review.tar.gz`.

Actual-container acceptance passed on Crest on 2026-11-24 (controller job
1575945; all five tasks completed with exit code zero, no retries; elapsed
1m 16s). The returned `fcs-adaptor-acceptance-01.review.tar.gz` records:

- Source SHA256 `276d661078d16301963b2dc6eaf50149e9bfd3fa8f58ff8019876325998de938`.
- 418 source sequences and 422 cleaned sequences; five terminal trims, five
  internal splits, one excluded sequence and exactly 1,102 bases removed.
- Exact retained bases, identifiers and lengths after uppercase conversion;
  zero adaptor findings on the cleaned FASTA.
- FCS-adaptor 0.5.5 local image launched with Singularity CE 4.4.2, and the GX
  cleaning container; reports published separately under contig and final FCS
  folders.

This accepts the real adaptor screen/cleanup/rescreen path, not a full GX
database screen or the integrated assembly/reporting route. A subsequent bounded
real assembly check remains necessary to verify the integrated contig and final
reports, exact FASTA/depth matching and FCS-colored GC/coverage plots.

## Bounded real-assembly integration check

After syncing, submit from the CREST project root:

```bash
mkdir -p comparisons/fcs-assembly-acceptance-01
sbatch --parsable \
  --export=ALL,FCS_ADAPTOR_IMAGE=/work/birdlab/singularity_cache/ncbi-fcs-adaptor-0.5.5.img \
  gcl_genome_assembly/scripts/comparisons/run_fcs_assembly_acceptance.sh
```

The launcher resolves the complete CTlk hap1 original GX input, finalized FASTA
and sample HiFi FASTQ from successful tasks in
`genome_assembly/pipeline/pipeline_trace.txt`. Ambiguous or absent retained inputs
stop preparation; no arbitrary task is selected.
Scratch-only input links are resolved through the exact source/destination
binding retained in the task's `.command.run`; the wrapper is never executed.
Missing upstream files still stop preparation rather than substituting cleaned
sequence. Explicit `--trace`,
`--source-fasta`, `--final-fasta`, `--hifi-fastq` and `--gxdb-dir` overrides are
available as arguments after the submission script. The existing full GX database
is reused; no nt or protein databases or BlobTools searches are involved.

This calls the maintained cleanup, final audit, diagnostics and reporting
workflows. It screens all sequence in one haplotype at each stage and computes
fresh full-read mapping for both exact plotted frames. It does not rerun assembly,
scaffolding or finalization: the final screen audits the actual retained finalized
FASTA, independently of the newly cleaned contigs. A full main-pipeline resume
remains separate from this bounded integration check.

Results, controller cache and work files stay under
`comparisons/fcs-assembly-acceptance-01/`. The existing stage-relative output layout
is retained inside `results/`. The launcher verifies FASTA preservation, report
scope, metrics row counts, source/taxonomy checksums and rendered PNG/PDF files,
then makes `comparisons/fcs-assembly-acceptance-01.review.tar.gz` with reports,
plots, provenance and logs. Biological findings do not fail this integration
check automatically; review the final recommendations and rendered figures.
Resubmit the same command to use this run's Nextflow cache.

The bounded real-assembly check passed on Crest on 2026-11-24 (controller
1576138; returned `fcs-assembly-acceptance-01.review.tar.gz`). Eight screen,
cleanup and mapping tasks loaded from cache; four evidence/report tasks completed
with exit code zero. The resumed run took 7m 2s, including environment setup;
the original mappings each took approximately 41 minutes.

Both stage reports and PNG figures were inspected. All 2,267 contig-stage records
and all 365 finalized records were plotted, with none omitted for missing GC or
depth. Metrics cover 953,537,721 bp and 953,579,871 bp respectively; these are
different assembly stages, not a claimed cleanup delta. Source/taxonomy checksums
and row counts passed, the finalized FASTA was unchanged, and both stages had
zero actionable GX/adaptor findings. The removed-sequence FASTA link is absent
from the review archive by design because FASTAs are excluded; it remains a
production output. Other report links resolve within the archive.

This accepts the maintained component workflows on a complete real haplotype,
including real screening, fresh exact-frame read mapping, R figures, stage reports
and ordinary cache reuse. The remaining integration gate is a normal production
entry-point resume to confirm cohort collection and main-report links. Optional
full-database BlobTools acceptance remains deferred.

## Full-cohort category/report acceptance

After syncing, submit from the CREST project root:

```bash
mkdir -p comparisons/fcs-cohort-acceptance-01
sbatch --parsable \
  --output=comparisons/fcs-cohort-acceptance-01/slurm-%j.out \
  --export=ALL,FCS_ACCEPTANCE_RUN=fcs-cohort-acceptance-01,FCS_ADAPTOR_IMAGE=/work/birdlab/singularity_cache/ncbi-fcs-adaptor-0.5.5.img \
  gcl_genome_assembly/scripts/comparisons/run_fcs_assembly_acceptance.sh --cohort
```

This discovers all finalized assemblies in the production trace and requires the
initial FCS cohort to match. Both haplotypes remain separate; read inputs are
deduplicated by individual only for channel routing. It uses complete retained
pre-FCS and final FASTAs, the existing GX database and HiFi reads. It reruns the
bounded FCS/mapping/report route for that cohort without assembly, scaffolding or
BlobTools searches. This new run has its own cache; the earlier single-assembly
cache is not automatically imported. Repeat the same command to resume it.

Return `comparisons/fcs-cohort-acceptance-01.review.tar.gz`. The combined index is
`results/overview/decontam/contamination_report.md`; individual stage reports
remain under `results/contig/decontam/` and `results/final/decontam/`. The launcher
checks exact source/taxonomy provenance, preserved final FASTAs, cohort membership,
category partitions, record counts and expected figures. Visual review of the
cohort figures and main production-resume acceptance remain separate checks.

## Method references

- [NCBI FCS-adaptor quickstart](https://github.com/ncbi/fcs/wiki/FCS-adaptor-quickstart)
- [NCBI FCS-GX quickstart](https://github.com/ncbi/fcs/wiki/FCS-GX-quickstart)
- [NCBI FCS output definitions](https://github.com/ncbi/fcs/wiki/FCS-GX-output)
- [NCBI's published adaptor fixture calls](https://github.com/ncbi/fcs/blob/main/examples/FCS_combo_test.fcs_adaptor_report.txt)
