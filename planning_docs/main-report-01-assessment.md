# Main report first-pass assessment

The retained main-report-01 run completed successfully. Its integration audit passed
with upstream assembly and graph tasks cached. The 6,182 final scaffold painting
partitions match the name-map dictionaries and cover their sequences without gaps
or overlapping intervals. Four verified cuts and two retained joins are represented
correctly. This accepts execution and evidence consistency, not the final report design.

## Accepted structure layout

The user selected Option 2: individual scaffold lengths and chromosome painting above,
with all unplaced sequence summed below on a separate scale. This is implemented in
the production R renderer, using pinned patchwork 1.3.0 with ggplot2 3.5.2.

- Select up to the configured number of longest chromosome-named scaffolds, then sort
  by chromosome and piece number. Assemblies without chromosome names show their
  longest individual pieces instead.
- Sum every unplaced piece, including those outside the upper display. State that
  these are independent pieces and that the upper and lower views are not additive.
- Keep the two haplotypes together and preserve common scaffold scales within species.
- Use discrete scaffold axes to prevent labels leaking from adjacent facets. Order
  the reference chromosome legend naturally and identify the reference dynamically.
- Preserve all inventory and painting evidence; record displayed rows and complete
  unplaced counts/lengths for reconciliation. No alignment or biological rules change.

The revised production figure was rendered and visually inspected using the retained
main-report-01 results. Nine renderer cases passed, including real-data reconciliation,
disabled sections, primary assemblies and unnamed-only assemblies. The production
Nextflow reporting fixture also passed with multiple species.

## Second main-report review

Implemented for user review with retained real results:

- Native quality table adds assembly length and chromosome-assigned percentage,
  calculated from complete final name maps, including composite names. This is
  deliberately distinct from alignment painting or chromosome-scale eligibility.
- Plot panels state the percentage represented by displayed individual bars.
- Cut locations include original scaffold names and formatted coordinates.
- Supplied/skipped/incomplete sample outcomes are visible; successful tables expand.
- The main pangenome section presents individual-level sharing from the structured
  sharing summary. A generated per-species detail index links all enabled biological
  analyses and technical validation, retaining the existing representation flags.
- Contact-map resolutions are grouped per assembly. Evidence headings and organelle
  map/annotation links are clearer. Disabled evidence categories remain absent.

Validation: nine renderer checks passed, including independent reconciliation of
assembly lengths, assigned percentages, QC metrics and exact unplaced totals with
source tables; retained artifacts are linked once. The multi-species Nextflow
fixture exercises structured sharing-summary delivery. The production main workflow
compiled and completed its input-validation checkpoint. The real figure was visually
inspected. No new assembly, alignment, graph or chimera evidence was generated.

Preview: `comparisons/main-report-01-inspection/results/assembly_report.preview.md`,
with new assets in `reports/main-preview/`. The original retained report is preserved.

## Individual-centered tables and evidence navigation

The next requested revision adds the measured combined Merqury completeness before
parenthesized haplotype values, e.g. `98.6 (78.5 / 78.4)`. Missing combined values stay
explicitly unmeasured; they are never estimated from haplotype completeness. Primary
assemblies retain one value. The exact combined source values are retained separately.

Downloads now have one row per individual, dynamic haplotype columns containing FASTA
and name-map links, and an available mitogenome GenBank column. Existing final snail
SVGs are explicitly registered from their production output channel. A native Markdown
table embeds clickable snail plots per haplotype and the existing circular mitogenome
map per individual. These figures are reused, without changing their renderers.

The evidence index is grouped by assembly quality, chromosome structure/comparisons
(dotplots and ribbons adjacent), organelles, pangenome, and additional reports. Existing
section anchors and retained links remain accessible; unavailable themes stay absent.

Ten renderer cases passed, including missing combined completeness and the retained
real run. The production multi-species reporting fixture includes final-snail delivery,
and the full main workflow passed its compilation/input-validation check.

Defer further biological pangenome highlights and postponed post-pangenome assembly
supplements to their scheduled analysis pass; revisit the main pangenome summary there.

## Focused detail reports and measurement audit

Downloads are now a normal section above assembly quality, after a generated contents
list. Species names accompany taxids in structure figures and pangenome headings.

Completed manual CUT/RETAIN/DEFER choices take precedence over original detector advice
in the chimera overview and event cards. Verified cuts require exact reconstruction
and a matching assembly, action ID, cut coordinate and assessment fingerprint.
Original recommendations remain explicitly labelled. Decision presentation changes
do not alter candidate discovery, evidence scores or cutting rules.

`assembly_stage_changes.md` gives signed differences between consecutive available
measurements, with explicit previous/current stages and source metrics retained in
`stage_changes.tsv`. It covers QUAST length/N50, BUSCO, k-mer completeness/QV, Hi-C
retention and trans:cis. Missing measurements are not imputed; measurement/filtering
differences are explained. Every reported previous/current value was checked against
the real source QC table. Existing trend plots are retained.

The mitochondrial source counters used text occurrences, counting paired gene/tRNA
qualifiers twice and gene/CDS names twice. A separate read-only ANNOTATION_AUDIT task
uses pinned Python 3.11 / Biopython 1.86 to count distinct annotated locations and
typed CDS/tRNA/rRNA features, with a retained locus audit. All five real GenBanks
contain 37 loci: 13 protein-coding, 22 tRNA and 2 rRNA. Original annotations and
assembly tasks are unchanged. Raw assembly-summary qualifier counts remain retained
as source records but are not displayed as biological gene counts. Reporting this
audit separately avoids invalidating expensive organelle assembly caches.

TIDK total lengths are inferred window spans, so the detailed report labels them
accordingly and distinguishes final-scaffold detection from preceding telomere
extension. Exact final lengths remain in the main inventory.

The old detailed report now links to main downloads/metrics, preserves contiguity
assessments, organism/read profiles, supporting maps/comparisons, QC plots and methods,
and links to focused stage, annotation and per-species pangenome reports instead of
concatenating all graph reports. Optional figure/organelle/telomere sections are omitted
when absent. The retained review replay shows four verified cuts and two retentions
before original detector classifications, without regenerating its accepted figures.

Validation includes ten general renderer cases plus the real-data replay, three
detailed-renderer cases, four annotation-locus tests, six chimera presentation/manual
review tests, the enabled multi-species workflow with the real annotation process,
the disabled workflow, and production main compilation/input validation. The annotation
workflow was tested in its pinned environment. Real structure figures were visually
inspected and local preview Markdown/figure links checked. A full native Crest rerun
of these latest reporting changes remains the next integration check.

## Final read-support and report integrity pass

The main report now includes per-individual final mapped-read depth, Hi-C retention,
trans:cis and contact-map links, with haplotype values kept separate. Read source comes
from actual accepted-input metadata; short-read depth is not labelled HiFi. Values are
copied from final QC records and independently checked against the source CSV. Captions
distinguish mean mapping depth from unique molecules and explain scaffold dependence
of the cis/trans ratio. Zero measurements remain distinct from unavailable values.

A concise Needs attention table combines unresolved/deferred decisions, incomplete
chimera assessment scope, incomplete assembly outcomes, graph representation review
flags and missing expected report inputs. Completed cuts/retentions are excluded.
These are reporting dispositions, not new biological cut recommendations or thresholds.

Expected-output checks use run controls plus each sample's available input types.
They distinguish PRESENT, DISABLED, NOT_APPLICABLE and MISSING in
`report_completeness.tsv`. They cover final QC/BUSCO/k-mer/read evidence, snail plots,
each requested Hi-C map resolution, eligible comparison pairs, recovered organelle
annotations/audits, chimera registry membership and admitted graph reports.

`py_scripts/validate_report.py --root <full-results>` additionally checks registered
published artifacts and local file/heading links in the main report, generated detail
pages and registered Markdown reports. The normal batch launchers invoke it only after
Nextflow exits, when publication has finished. Missing expected outputs or broken links
give a nonzero exit. Diagnostics are `report_validation.md/json` at the results root,
outside cached task directories. Symlink publications into work directories are valid.
Compact review archives intentionally omit some large targets and are not substitutes
for the full results tree in this publication check. Direct Nextflow users can invoke
the same general-purpose validator after their run.

Validation: eleven general renderer cases plus the real-data replay, six publication
validator cases (missing outputs, broken anchors, detailed links, disabled reports,
disabled analyses and symlinks), enabled/multi-species and disabled production reporting
fixtures, full pipeline compilation/input validation, and batch launcher syntax checks.
No native Crest run of the latest reporting contract has yet been evaluated.

## Fragmented-assembly presentation preview

The current long-scaffold view hides smaller named fragments, while its unplaced total
is a naming category rather than absence of chromosome associations. In the retained
Pla data, approximately 193.8/290.8 Mb (hap1) and 199.4/280.0 Mb (hap2) of unplaced
sequence already has qualifying chromosome labels from final-reference painting.

A temporary comparison script created a read-only design preview for
final assemblies classified as non-chromosome-scale. It summed every disjoint painted
query interval by chromosome and final naming category, with explicit Ambiguous and
Unassigned bins. Totals reconcile exactly with the complete scaffold inventories.
It does not change alignments, filters, harmonized names, scaffolds or biological
decisions. These sums include redundant fragments and are not reference coverage or
chromosome-completeness estimates. The preview is in
`reports/layout-options/fragment_assignment_totals.png/pdf/tsv` within the inspection
results. The user rejected this preview on 2026-11-04 as overlapping existing/planned
pangenome presentation. Do not implement it in the production report. The temporary
renderer was removed; generated inspection artifacts record the rejected design.
The existing chromosome-structure layout remains unchanged.

## Actionable attention entries

Per the user's subsequent review, no-chromosome-scope cases (including both Pla
haplotypes) are excluded from Needs attention. They remain in the generated
`assessment_coverage.md` page, which explains that chromosome association labels
do not establish chromosome-scale scaffolds and that the screen limitation does
not identify a repair or cut to make. Other unavailable assessments remain visible.

The generic graph warning is replaced by individual flagged assembly/path rows.
The reporter consumes the existing `retention_review` flags and retained-path
fractions from `source_sequences.tsv`, plus the actual screen cutoff from
`source_audit.json`. These are forwarded from the existing source-audit output
channel; only flagged TSV rows enter the report specification. There is no new
threshold or reconstruction of biological evidence, and no graph recomputation.

Each attention row states the path/input length percentage and review screen,
links directly to source/path and mapping evidence, and asks the reviewer to
distinguish graph representation differences from omitted source sequence. The
section explains that path-length percentages are not missing-DNA percentages.
The retained real run gives five graph-path entries and no Pla entries. Renderer
and native multi-species workflow checks verify this presentation and exclusion
of unflagged paths. Main-preview assets have been refreshed from the same evidence.

Attention rows also carry an explicit Area column. Graph-path entries are labelled
Pangenome graph, separately from assembly join decisions, assembly completion,
chimera assessment and report completeness. The introduction states that these are
checks of graph representation, not newly identified assembly misjoins or assembly
edit recommendations. Input sequence length is the denominator, not whole-genome
assembly length. No evidence calculations or audit thresholds change.

## Remaining presentation review

Keep the approved native Markdown quality table and links to all metrics and the
interactive rankings. Continue reviewing the main report before revising detailed
reports. Candidates include making displayed fractions more prominent, showing the
original scaffold alongside cut coordinates, simplifying graph links and explaining
representation flags, and grouping repeated contact-map links by assembly/resolution.
These do not reopen accepted assembly, cut or graph evidence calculations.
