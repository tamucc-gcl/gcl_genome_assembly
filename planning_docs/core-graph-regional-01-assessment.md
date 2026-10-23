# Regional assignment diagnostic: source intent and cohort result

The Crest packet `core-graph-regional-01.review.tar.gz` completed with exit 0:
871 tasks, 865 cached and 6 completed, with no failed tasks. The patched build,
three patched assignments and summary all completed. The source diff contains
exactly the proposed `.first` to `.second` comparison correction.

## Coding error versus policy

The upstream header explicitly defines a coverage interval value as
`<aligned bases, reference contig>`. In `smooth_query_intervals`, the loop comment
describes grouping runs of consecutive mappings to non-dominant reference contigs;
`ref` is assigned the preceding interval's reference-contig ID (`.second`). The
condition compares the current interval's aligned-base count (`.first`) against
that ID. This is strong evidence of a field-selection coding error. A chance
numeric equality between a count and a chromosome ID should not determine whether
two intervals belong to the same chromosome run.

Upstream history checked on 2026-10-23:

- [Initial smoothing implementation, 2021-05-04](https://github.com/ComparativeGenomicsToolkit/cactus-gfa-tools/commit/8d8e241891b24d945ae7063e19574d28567003ed)
  introduces the mismatched comparison and the separate longest-interval gate.
- [Follow-up correction, 2021-05-06](https://github.com/ComparativeGenomicsToolkit/cactus-gfa-tools/commit/dde6e887b8253fdf9194bdaeadc918fcf1e803b6)
  adds the empty-candidate guard but retains the mismatched comparison.
- Current upstream main still contains the comparison. No maintainer confirmation
  of this proposed field correction has been obtained.

The condition `max_interval_length > min_chunk` is different: it explicitly requires
at least one sufficiently long premerged mapping interval in a non-dominant run.
At `-u 1000000`, a large chromosome arm fragmented into intervals no longer than
1 Mb can fail this condition despite substantial total mapped coverage. This rule
must be treated as an existing policy, not automatically called a bug or replaced
with total run length. Density, overlap resolution and later smoothing also matter.
The packet does not record all rejected candidate statistics, so it cannot identify
the precise rejection gate for each missing arm.

## What the new run actually demonstrates

- The two synthetic native-CLI regressions improve affected-arm assigned span from
  2.0 Mb to 3.4 Mb, with the ordinary chromosome control unchanged.
- Across the real cohort's 4,806 source/filter-policy comparisons, **none** change
  surviving mapped coverage or reported chromosome membership after the correction.
- Under isolated filtering, CTlk hap1 has 72,629,202 bp of input mapped coverage,
  but corrected regional assignment still retains only about 40.485 Mb on chr12.
  Hap2 has 79,621,204 bp of input mapped coverage but retains about 40.603 Mb,
  again only on chr12. The chr7 arms remain unresolved.
- All 4,806 source-coordinate audits pass; 16 ordinary chromosome controls show
  no decreases. This certifies the audited intermediate ownership, not final-graph
  accounting or biological completeness.

The field correction is a well-supported bug candidate, but **is not the demonstrated
solution to the cohort's missing chromosome arms**. Earlier language presenting it
as that solution was too strong. Keep it diagnostic. Before further changes, log the
candidate-run statistics and each rejection gate on the cached inputs, then separate
configuration suitability from genuine implementation defects. Do not rebuild the
production graph or relax thresholds merely to recover the expected example.
