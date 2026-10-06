# CTlk first adjudication packet — 6 October 2026

Status: source changes and static evidence reconciliation prepared. No new mappings, assembly runs, tests, or cuts executed. Job 1511195 was cancelled and incomplete; its outputs are not available in this local checkout.

## Structures to explain

| Structure | Existing observation | Competing explanation | Next decisive inspection |
|---|---|---|---|
| hap1 scaffold_1, 138.447 Mb, chr7/9/12 | Two major assignment-transition windows already occur within raw contigs; H01 is a separate weak physical round1 gap | Correct local sequence with mistaken assignment versus wrong internal path versus wrong connection of large blocks | Map distributed distinctive blocks to current sisters/six better peers and overlay exact raw-contig/gap boundaries; inspect H01 separately |
| hap1 scaffold_5, 74.766 Mb, chr4/14 | Assignment transition lies inside h1tg000153l; J05 has low MAPQ but a retained spanning molecule | Ambiguous chromosome assignment/copy or phase versus wrong internal path | Recover unitig/read-overlap path through h1tg000153l and identify the distinctive sequence on each side against peers |
| hap2 scaffold_5, 83.431 Mb, chr7/12 | Assignment transition lies inside h2tg000028l; H03 is a separate physical gap | Same three possibilities; shared pair labels do not establish shared boundaries | Compare its actual subregions and orientation with hap1, including both separate chr7 pieces; reconstruct h2tg000028l path |
| hap2 scaffold_7, 82.880 Mb, chr4/14 | Assignment transition lies inside h2tg000298l; H04 is a separate physical gap | Shared internal ambiguity versus wrong scaffold connection versus assignment | Reconstruct local unitig/read path; establish correspondence with hap1 chr4/14 boundary and all separate chr4 pieces |

Include chr4_2 in both haplotypes, chr4_3 only in hap2, and every substantial chr7/9/12/14 piece. Number suffixes do not identify matching homologous intervals. The six better chromosome-scale comparator assemblies are CBau/CLim/CMat hap1/hap2 (three individuals). CPla may inform local sequence correspondence but is not an informative chromosome-order vote.

## What can be reused now

- Assignment/foundation name maps and reference PAFs: four composites and existing chromosome assignment evidence. Reference-alignment endpoint ranges are not unique-copy coverage.
- Origin packet 1506290: exact local raw-contig placements, verified source boundaries, stage digests, and contig graph record summaries. Local extracted graph snippets lack the internal unitig paths needed to validate these transitions.
- Assessment packet 1506302: twelve suspect physical gaps, eighteen comparison gaps, five raw-contig support tracks, measured HiFi and Hi-C with coordinate provenance.
- Focus packet 1511158: multiscale contact results and retained flank comparisons. H01 has the strongest reported contradiction; H02 final gap has comparative support; J03/J04 need joint component review. These are boundary findings, not clearance of an entire composite.

History is optional diagnostic context. It shows local continuity through the five windows and some old scaffold splits, not proof that old final assemblies retained the current chromosome-scale composites. The production method must not depend on historical assemblies.

## Prepared source changes

The all-junction helper now inventories all boundaries separately and requires an explicit selection for evidence work. `--registry` selects flagged intervals only; it no longer implicitly requests measurements for every N-run. `--selection-packet` imports the bounded physical gaps/controls from a completed assessment packet with its assembly digests; explicit `--selection` tables are also supported. Selection rejects unknown IDs, stale coordinates/digests and unspecified roles. Anchor-confounding checks still use the full inventory, including unselected physical boundaries.

Hi-C collection indexes alternative partner counts by library/junction/side instead of repeatedly scanning all accumulated partner counts. It logs the start, every ten million input pairs, and completion of each assembly scan. This addresses an observed code inefficiency; it does not establish the bottleneck of the cancelled job. Contact streaming still reads a complete relevant source pairs file once; selection does not eliminate that I/O.

No sister mapping cache/resume implementation is added in this first patch. Completed unmodified old outputs are preserved; do not claim automatic reuse of cancelled-job sister mappings. Existing successful HiFi summaries can still be imported through the existing evidence-packet validation.

## First cluster action: presence-only hifiasm context inventory

After syncing source, from `/work/birdlab/GCL/spratelloides_delicatulus_genome`, run:

```bash
python3 gcl_genome_assembly/scripts/comparisons/inventory_hifiasm_context.py \
  --results comparisons/chimera-foundation-01 \
  --sample Sde-CTlk_104 \
  --out comparisons/ctlk-hifiasm-context-01
```

This lists retained sample graphs/caches, checks whether raw/processed unitig graphs are present, and copies only small diagnostic text. It does not run hifiasm, map reads, copy large graphs/caches, or search all work directories. If a specific surviving hifiasm task directory is already known, add `--task-dir /exact/task/path` to include its command, log and intermediates. Output directories must be new; do not overwrite a prior packet. Return `context_inventory.json` and the small `root_*` diagnostic files. Missing intermediates are a result, not a reason to resubmit the assembly immediately.

The next source step depends on this inventory: reconstruct from retained processed/raw unitig graphs and actual overlap paths if available; otherwise design one explicit regeneration using compatible retained HiFi intermediates, with a fresh treatment-specific Hi-C cache. Do not extrapolate from contig S/A records.

## Cluster verification and gated evidence run

Run focused implementation checks after sync:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover \
  -s gcl_genome_assembly/tests -p 'test_junction_review.py'
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover \
  -s gcl_genome_assembly/tests -p 'test_junction_assessment.py'
```

No new sister/contact batch is requested yet: inspect the context inventory and reuse the already completed evidence first. If a bounded sister refresh later proves necessary, use the foundation assessment results, `--registry gcl_genome_assembly/planning_docs/investigations/ctlk-junction-registry.tsv`, and `--selection-packet comparisons/junction-assessment-20261005-084623-1506302`. The latter packet selects 30 physical gaps/comparison joins, and the registry adds 11 flagged intervals; it is a CTlk pilot, not complete batch discovery. Omit `--hic` unless existing contact summaries cannot answer a named question. Supplying `--evidence-packet` reuses matching HiFi gap summaries; internal tracks remain separately reviewed evidence.

## Terminal outcome of this first stage

For each composite, identify whether the current uncertainty concerns chromosome assignment, an internal unitig/contig path, a physical scaffold adjacency, or more than one. Each implicated boundary receives a named test with available/missing input status. Preserve sequence and mark biological uncertainty separately from assembly action. No executable cut is justified simply by this packet.

## Follow-up to returned context inventory

The supplied context_inventory.json confirms only the standardized hap1/hap2 contig FASTAs/GFAs, a 28,955-byte hifiasm log, and versions.tsv in the published directory. No raw/processed unitig graphs or HiFi overlap caches were found there. This establishes absence from that directory, not loss from the task directory. The log itself has not yet been supplied locally, so its coverage/phase diagnostics are not verified.

The retained foundation trace identifies CTlk hifiasm as cached task hash `b4/4c0468`, native job 1499288, originally submitted 23 September. Its cached provenance matters: the October foundation pass did not newly run hifiasm. The actual command/version of the September task is needed before attributing any output to current defaults.

The collector now resolves just that sample's hifiasm task hashes from a trace. After sync, run this small follow-up from the cluster project root:

```bash
python3 gcl_genome_assembly/scripts/comparisons/inventory_hifiasm_context.py \
  --results comparisons/chimera-foundation-01 \
  --sample Sde-CTlk_104 \
  --trace comparisons/chimera-foundation-01/pipeline/pipeline_trace.txt \
  --work-root work \
  --out comparisons/ctlk-hifiasm-context-02
tar -czf comparisons/ctlk-hifiasm-context-02.tar.gz \
  comparisons/ctlk-hifiasm-context-02
```

Return this small archive: it contains the inventory plus available small logs, versions, command and annotations, without copying large graphs or caches. If the trace used a different work root, substitute that actual root. Missing or ambiguous task paths are recorded; the collector will not choose among ambiguous matches or recursively search work. Tests are prepared in test_hifiasm_context_inventory.py but have not been executed locally.

Static inspection of the existing assessment PAF rows adds one useful limitation: at the descriptive MAPQ>=20 and >=50 kb query-span screen, J05 has an 89,115 bp alignment in CMat hap1, but its query coordinates 47,363–136,478 do not span the full transition core (100,000–237,672 in that extracted window). J07 has a 51,938 bp alignment in that same target, but query 286,912–338,850 is downstream of its core (100,000–238,862). Neither validates the disputed chromosome boundary. J01 has several much longer alignments in current peers, but span/order/CIGAR and phase still need review; a long PAF span is not a verified uninterrupted bridge. This was inspection of retained rows, not a new mapping run or calibrated uniqueness test. Source targets must be resolved from sequence_comparisons/sources.tsv and the corresponding name maps, not from peer file indices alone.
