# Pangenome preprocessing changes: implementation and biological rationale

Accepted from the real production build reviewed on **2026-10-25**. This document
describes the production differences from the pinned upstream tools, the evidence
that motivated them, and their limits. It distinguishes implementation corrections
from pipeline policy choices. None of these changes makes additional assembly cuts.

## Scope and exact tool versions

| Component | Audited version | Production treatment |
|---|---|---|
| Cactus | 3.2.1, commit `d7e5fc1b51327a5a4cdb494c85b0607e3d81a713` | Installed package is unchanged; task-local preprocessing adapters |
| cactus-gfa-tools / `rgfa-split` | commit `f9b999e54c49a05ad6e9d44e1832c4a8421f28f2` | Compile a separate executable with one comparison corrected |
| hal2vg / `filter-paf-deletions` | v1.1.8, as selected by the pinned Cactus build recipe | Installed executable and its threshold arguments are unchanged; query grouping changes |
| samtools and `cactus_sanitizeFastaHeaders` | Executables supplied by the configured Cactus image | Execute a native fragment-header contract check before graph building |

The image is `docker://quay.io/comparative-genomics-toolkit/cactus:v3.2.1`.
The native executable hashes, source manifest, compilation command, effective
configuration hash and installed Cactus commit are recorded per run. The container
tag is not an immutable digest; do not treat the tag alone as complete provenance.
These findings concern the pinned versions, not every past or future upstream release.

Primary implementation references:
[Cactus build recipe](https://github.com/ComparativeGenomicsToolkit/cactus/blob/d7e5fc1b51327a5a4cdb494c85b0607e3d81a713/build-tools/downloadPangenomeTools),
[regional assignment](https://github.com/ComparativeGenomicsToolkit/cactus-gfa-tools/blob/f9b999e54c49a05ad6e9d44e1832c4a8421f28f2/rgfa-split.cpp),
[deletion filter](https://github.com/ComparativeGenomicsToolkit/hal2vg/blob/v1.1.8/filter-paf-deletions.cpp).

## 1. Chromosome-isolated deletion filtering

**Classification: pipeline policy workaround, not an upstream binary patch.**

### Issue

Native deletion filtering protects graph topology from implausible jumps between
mapping blocks. In the pinned implementation, a query's dominant reference path
governs the calculation, while the record interval selected for removal can also
contain mappings to other chromosomes. Their support need not contribute to the
dominant-path removal cost. On the investigated composites this removed a coherent
chromosome arm before chromosome assignment could examine it.

This is a concrete inappropriate outcome for these multi-chromosome queries. It
does not mean the original filter's general purpose is wrong or filtering should
be disabled for all mappings. The native-control replay reproduced the previous
exported PAF byte for byte, establishing that the effect was not a replay artifact.

### Change

`py_scripts/pangenome_cactus_guard.py` determines a reference-chromosome owner for
each usable graph node. Rank-zero reference nodes are boundaries. Only
non-reference nodes are connected into components when propagating ownership;
a graph edge between two reference chromosomes cannot merge their labels.
Used unanchored nodes or components bordering multiple chromosomes cause a failure,
rather than an inferred chromosome assignment.

Each original `(query, reference chromosome)` pair receives a temporary query
identifier. The same native deletion-filter executable runs on these lanes with
the same arguments. The original query identifiers are then restored. Original
query lengths, coordinates, target fields and alignment tags are preserved.

This also applies filtering independently inside non-dominant chromosome lanes;
it is **not** guaranteed to produce the same retained mappings as stock filtering.
It is not equivalent to a minimal repair of the native removal-cost calculation.

### Biological rationale and limitations

Reference chromosome membership is a coordinate system, not a requirement that
every input scaffold have exactly the reference's chromosome arrangement. A true
fusion, translocation or other rearrangement can put long homologous blocks from
different reference chromosomes in one assembly sequence. An erroneous assembly
join can do the same. Neither case justifies erasing one block before evaluating
its own within-chromosome alignment evidence.

Separating filter comparisons by chromosome prevents a decision about an implied
deletion on one reference chromosome from collateral removal of another's blocks.
Within-chromosome filtering remains active. This is our biological rationale for
the policy, supported by the controlled tests; it is not a claim that every rescued
mapping is orthologous or every composite is a genuine rearrangement.

Paralogs, repeats and incorrect reference anchors remain possible. Ambiguous graph
ownership is deliberately unresolved. This guard is not a general orthology solver
and can require a different validated policy for structurally complex cohorts.
Assembly join decisions remain in the separate chimera evidence/manual-cut workflow.

## 2. Correct the regional grouping comparison

**Classification: implementation defect with a one-line source correction.**

The interval value stores `(aligned-base count, reference chromosome ID)`.
The loop's `ref` variable holds the previous chromosome ID, but the original
comparison uses the aligned-base count. This splits consecutive same-chromosome
intervals into artificial runs.

```diff
- if (interval.value.first != ref || clip_candidates.empty()) {
+ if (interval.value.second != ref || clip_candidates.empty()) {
```

The surrounding code assigns `ref = interval.value.second`, and its stated purpose
is to group consecutive non-dominant-chromosome mappings. Comparing chromosome IDs
therefore restores that intended grouping. No thresholds, interval coordinates,
coverage calculations or final decision rules change in this source patch.

Biologically, a long homologous chromosome block can consist of many local alignments
with different lengths. Changes in local alignment length are not evidence that
chromosome identity switched. Grouping by chromosome ID lets the existing coverage
and interval gates judge the coherent block. Other chromosomes still separate runs.

The original vendored sources remain unchanged under
`vendor/cactus-gfa-tools-regional/`; their manifest and MIT license are retained.
`pangenome_cactus_policy.py build` verifies the source hashes, applies exactly one
replacement in a copied build directory, compiles a separate executable and runs
a real CLI regression. The installed container executable is never overwritten.

The correction alone did **not** recover the missing arms under the 1 Mb regional
setting. This negative result is important: grouping and the next gate were
independent causes, not interchangeable explanations.

## 3. Enable regional assignment at 25 kb with explicit softclip accounting

**Classification: supported parameter/policy choice, not a source defect.**

### Issue and change

The stock first-pass invocation assigns whole queries. A query with substantial
support for two chromosomes can fail whole-query uniqueness. Regional mode offers
a supported way to assign non-overlapping source subranges instead.

The regional diagnostic initially used a 1 Mb chunk setting. Its secondary-run
gate requires at least one individual interval **strictly longer** than the chunk,
as well as sufficient coverage. The motivating chr7 runs spanned 32.83/40.18 Mb at
about 95% alignment density, but their longest intervals were only 49,750/102,554 bp.
They therefore failed the 1 Mb gate despite strong aggregate evidence. The dominant
chromosome did not face the same secondary-run test.

Production now supplies `-u 25000 -s` through the named parameter
`pangenome_regional_chunk_bp=25000`. Coverage, uniqueness, MAPQ and other native
gates are preserved. This is not a new rule that any 25 kb match establishes a
chromosome arm. The setting also affects smoothing and endpoint softclipping;
those coupled effects are part of the policy, not hidden implementation details.

### Biological rationale and tradeoff

Base-level homology across diverse assemblies is often represented by many local
alignments rather than a single megabase alignment. Aggregate coherent chromosome
support should not be discarded solely because its individual records are shorter
than a megabase. Regional assignment can represent both parts of a rearranged
scaffold without pretending that the original source was two independent assemblies.
It partitions graph-alignment subproblems; it does not cut or validate the scaffold.

The 25 kb choice is empirical, not a universal biological constant. The 100 kb case
recovered only one motivating arm, whereas 25 kb recovered both. It is configurable
for future validation in other taxa or assemblies. Smaller accepted regions and
reference-dependent assignments remain possible; coverage alone does not prove
orthology or rule out paralog placement.

The chosen setting excludes some additional source ends. Unassigned ends are not
automatically sent to an "other" graph bin. They remain in the source FASTAs and
are listed explicitly with original coordinates. A graph built from selected
chromosome partitions is not claimed to contain every input assembly base.

## 4. Match fragment alignment names to sanitized FASTA names

**Classification: pinned Cactus interface defect; corrected in the pipeline adapter.**

The first production attempt reached native Cactus alignment but failed because
the PAF retained faidx names such as `sequence:11-50`, while the FASTA sanitizer
produced `sequence_sub_10_50`. The pinned Python rename helper uses an unescaped
`+` in a basic sed expression and does not convert ordinary numeric ranges. The
native sanitizer implements the conversion independently, exposing the mismatch.
See the [helper](https://github.com/ComparativeGenomicsToolkit/cactus/blob/d7e5fc1b51327a5a4cdb494c85b0607e3d81a713/src/cactus/shared/common.py)
and [sanitizer](https://github.com/ComparativeGenomicsToolkit/cactus/blob/d7e5fc1b51327a5a4cdb494c85b0607e3d81a713/preprocessor/cactus_sanitizeFastaHeaders.c).

After auditing the native split output, the production adapter normalizes fragment
PAF query identifiers to the sanitizer's `_sub_start_end` convention. Extraction
lists keep faidx names. PAF coordinates remain fragment-relative, and lengths,
targets and alignment tags stay unchanged. Alias collisions fail explicitly.

This is a lossless naming/coordinate-convention reconciliation, not a biological
threshold change. One-based inclusive `11-50` and zero-based half-open `[10,50)`
describe the same 40 bases. The original source interval remains recoverable.

Before expensive construction, a tiny probe uses the **installed** samtools,
upstream helper and native Cactus sanitizer, together with the production PAF
normalizer. FASTA and PAF identities and lengths must agree. The successful Crest
run passed this check. Local tests reproduce the original helper's failure and
confirm that other mapping fields are unchanged.

## Validation and observed biological representation

The broader motivation—representing assembly-derived structural variation while
reducing dependence on one linear reference—is consistent with the
[Minigraph-Cactus methods paper](https://www.nature.com/articles/s41587-023-01793-w).
The paper does not independently validate these particular wrappers or the 25 kb
setting. Their acceptance rests on the tests and real production evidence below.

| Evidence | Finding |
|---|---|
| Native-control filtering | Reproduced previous exported PAF byte for byte |
| Grouping/gate experiments | Corrected grouping makes coherent runs; the independent longest-interval gate still explains rejection at 1 Mb |
| 18-case filter/grouping/chunk comparison | Guarded, corrected 25 kb recovers both arms; 98.74 Mb mapped gain across the selected sources, no mapped losses or new chromosome labels on 16 ordinary controls |
| Selected coordinate audits | All 19,224 source checks pass, with no duplicate/overlapping assignment |
| Production-02 | Exit 0; 17 completed tasks, 840 cached, no failures/retries; native fragment check passes |
| CTlk hap1 composite | CLIP path length 40.64 → 73.00 Mb, 96.82% of its 75.39 Mb input |
| CTlk hap2 composite | CLIP path length 40.81 → 79.99 Mb, 95.88% of its 83.43 Mb input |
| All-source accounting | 6,182 sources; no assigned source span missing from FULL, no FULL span outside assigned partitions, no CLIP span outside FULL |

The two recovered chr7 partitions exactly match the selected diagnostic's source
bounds. Five pre-existing low-retention chromosome paths remain flagged, with
essentially unchanged CLIP lengths. They are unresolved representation limitations,
not demonstrated biological deletions. The changes do not establish the truth of
the retained chr7/chr12 joins or imply that all graph sequence is biologically correct.

The selected diagnostic examined 1,602 sources; the final production audit covers
6,182. Their excluded-base totals have different denominators. In production,
203,087,437 ambiguous bases plus 254,171,095 bases outside partitions are explicitly
excluded before chromosome alignment, 4.44% of 10,304,718,010 input bases.

FULL has 2,372,383,979 distinct graph bases; CLIP has 1,959,888,337. Distinct segment
totals count each graph node once. Haplotype traversal lengths and source-coordinate
spans are different measurements, and are not added to distinct sequence totals.
Non-overlapping source assignments prevent allocating one source interval to two
chromosome subproblems. This does not erase genuine repeated sequence or prove that
the graph cannot contain biological/repeat-induced duplication.

## Permanent implementation and review outputs

Production ownership:

- `modules/pangenome_sequence_tools.nf`: verified cached corrected-tool build.
- `py_scripts/pangenome_cactus_policy.py`: native probe, task-local wrappers,
  parameter application, fragment names and preprocessing accounting.
- `py_scripts/pangenome_cactus_guard.py`: chromosome ownership and query lanes.
- `modules/pangenome_sequence_audit.nf` and `py_scripts/pangenome_sequence_audit.py`:
  original-source coordinate integrity, stage discrepancies and distinct-node totals.
- `modules/cactus_pangenome.nf`: normal Cactus build using these adapters.

Under each `pangenome/<taxid>/` output:

| Artifact | Purpose |
|---|---|
| `sequence_policy/policy.json` and `tool_build.json` | Effective policy, native probe and exact tool provenance |
| `sequence_policy/filter_*.lanes.tsv.gz` | Temporary query IDs linked to original source/chromosome |
| `sequence_policy/split_*.regions.tsv.gz` | Original-source partition coordinates |
| `sequence_policy/split_*.excluded.tsv.gz` | Ambiguous/unassigned intervals and reasons |
| `sequence_policy/sequence_accounting.md` | Preprocessing totals and exclusions |
| `representation/sequence_graph_audit.md` | Bounds/overlap checks, stage differences, distinct-node totals |
| `representation/source_stage_accounting.tsv` | Every source through preprocessing, FULL and CLIP |
| `representation/source_report.md` and `clipping_report.md` | Mapping/path-length limitations and residual review flags |

Paths are relative to the repository or the run output as indicated. Detailed
acceptance records are in [production-02 assessment](../planning_docs/core-graph-production-02-assessment.md),
[chunk comparison](../planning_docs/core-graph-chunks-01-assessment.md) and
[gate assessment](../planning_docs/core-graph-gates-01-assessment.md).

## Diagnostic retirement and future upstream replacement

The completed assignment, deletion-policy, grouping, gate-trace and chunk-matrix
workflows, switches, adapters and experiment-only tests have been removed from the
active pipeline. General ownership/lane tests were retained against the production
helper; the native grouping regression remains in the production tool build.
The pinned vendor sources, build environment, production adapters and permanent
audits are still required. Past reports/data have not been deleted.

The pre-cleanup diagnostic code is recoverable from repository snapshot
`0b2178dc4765e28437b39bedf78fbab687e2f541` using Git history, with its original
planning documents describing the experiments. This is a source recovery point,
not a claim that every historical Crest run used that exact commit.

An upstream upgrade must repeat known-answer grouping, native fragment handoff,
source-coordinate and real-cohort arm-retention checks. Reassess whether its native
deletion policy handles multi-chromosome queries appropriately; a new version number
alone is insufficient. Retire each workaround only after its replacement passes.
Use ordinary Nextflow caching. Removing these default-off diagnostic branches does
not require changing the accepted Cactus construction command or rebuilding the
accepted graph. Keep the same work directory, source inputs and production settings.
