# CTlk localized-boundary batch

Submit from /work/birdlab/GCL/spratelloides_delicatulus_genome after syncing the repository:

```bash
mkdir -p logs
sbatch gcl_genome_assembly/scripts/comparisons/run_ctlk_boundary_batch.sbatch ctlk_analysis
```

Six tasks run independently without a concurrency cap. Each requests 16 CPUs, 96 GB, and 12 hours. Dependencies are the existing Python/samtools/minimap2 analysis environment. No new assembler run, genome-wide mapping, Hi-C pairs scan or package installation is required.

Defaults: previous evidence comparisons/ctlk-adjudication-1511274; assessment comparisons/chimera-foundation-01. Supply different actual Crest paths as second and third arguments if relocated. They must be the native retained job directories, not merely an extracted review archive without its FASTAs.

| Task index | Test |
|---:|---|
| 0 | H01 isolated gap-unjoin validation |
| 1 | J01 smaller-anchor localization, raw boundaries, and adjacent H01 read test |
| 2 | J02 narrow seam endpoints/interval and local read probes |
| 3 | J05 smaller-anchor localization and local read probes |
| 4 | J06 smaller-anchor localization and local read probes |
| 5 | J07 smaller-anchor localization, raw boundaries and local read probes |

## H01 output

Creates comparisons/ctlk-boundaries-JOBID/H01/unjoined.fa as an isolated copy of the complete assessment assembly. Only scaffold_1 is split, at the right end of its [63051225,63051325) gap. The 100 Ns are preserved at the left record's terminus. All original bases and other scaffold sequences remain; indexed output record lengths and total bases are verified, and the original/reconstructed scaffold sequence SHA256 values are recorded.

Existing exact peer-anchor coordinates are relabeled as left/right/straddling. Existing 250 kb Hi-C bins are relabeled; the bin spanning the cut remains ambiguous. The source contact profiles are copied for review. Relabeling does not manufacture stronger molecular support or improve an evidence score by itself. The right record still contains chr7/12 material: this is one adjacency correction hypothesis, not a complete resolution of the composite.

This task does not rerun YaHS. Recreation of the join is tested later with explicit adjudicated constraints, after the correction is accepted. Production files are never overwritten.

## Internal-junction outputs

Each job uses an exact, full-window alignment of the original raw contig to the assessment before projecting coordinates. Non-unique/non-exact placement fails rather than exporting guessed coordinates.

Raw graph placements and native A records are reused from the previous packet for both and HiFi-only treatments. Exact segment endpoints within a core become physical-boundary TESTS, not cut instructions. J02 explicitly tests original h1tg000133l positions 4085098 and 4091627, and the 6529 bp interval between them. A raw-segment alignment endpoint does not by itself establish an unsupported assembly edge.

Each junction maps 2/5/10 kb tiles every 10 kb over its core plus 750 kb on each side against six peer haplotypes and both current sisters. Secondary placements are retained. Conservative full-span/identity/competitor screens remain in place; short anchors need consistent independent tiles and independent individuals before establishing chromosome identity. Raw PAFs permit review of rejected fragmented alignments. No midpoint is designated as a breakpoint.

Retained competitive HiFi placements are inspected at actual graph-boundary tests and at 1 kb spaced local-coverage probes across each core. The latter are assay positions, explicitly not proposed cuts. CIGAR blocks exclude insertions/deletions when counting aligned flanks; target and original query coordinates are validated including reverse orientations. Molecule names spanning >=1 kb on each side are counted separately at any MAPQ and MAPQ>=20. All emitted competing placements for selected molecules are exported for phase/repeat review. Existing nearby comparison gaps are reevaluated with the same metric.

The local-coverage probes make the read-span question feasible; no molecule is required to cover a broad 100-300 kb core. Nevertheless, BAM ascertainment, residual hard clipping, unrepresented alternative paths and low mapping uniqueness remain limitations. Zero counts require informative comparison controls before supporting a break. Existing reads are reused, not newly retrieved or remapped; no graph-boundary support is claimed independently of those limitations.

## Return and interpretation

Return the six comparisons/ctlk-boundaries-JOBID/*.review.tar.gz files and any failed-task logs. Full FASTAs remain on Crest; review packets contain status/accounting, coordinates, peer PAFs, narrow read summaries and competing placements. No automatic cut file is produced.

Interpret narrow boundaries using independent chromosome anchors, graph overlap/path structure and molecule-specific competition together. Adopt the H01 unjoin only after its before/after accounting/correspondence review. Internal cuts remain prohibited where localization is inadequate; ambiguous regions are reported explicitly rather than approved as chromosome fusions.

## Validation

Five unit tests pass locally: gap/base preservation, rejection of a non-gap cut, insertion-aware lifting, reverse PAF lifting, and rejection of inconsistent CIGAR/endpoints. Python 3.11 syntax checks pass. Coordinate accounting was also checked on 20,298 retained HiFi/graph alignments from the previous packet. The cluster launcher reruns unit tests before each task. Full external-tool integration remains a cluster check; no jobs were submitted locally.
