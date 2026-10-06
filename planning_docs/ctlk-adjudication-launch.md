# Launch CTlk parallel adjudication batch

Sync the repository to Crest, then submit from /work/birdlab/GCL/spratelloides_delicatulus_genome:

```bash
mkdir -p logs
sbatch gcl_genome_assembly/scripts/comparisons/run_ctlk_adjudication.sbatch ctlk_analysis
```

The existing ctlk_analysis environment (Python, samtools, minimap2) is sufficient; no new environment or assembler rerun is required. Eight tasks can run simultaneously without an array throttle, subject to Crest resources. Each requests 16 CPUs, 96 GB and 24 hours. The Hi-C tasks each scan their retained pairs once; progress appears every ten million pairs. Graph and peer-mapping progress also appears in logs.

Default retained directories:

- packet: comparisons/junction-assessment-20261005-084623-1506302
- assessment: comparisons/chimera-foundation-01
- native hifiasm results: comparisons/ctlk-fusion-batch-1511216

If these directories were relocated, pass their actual paths in that order after the environment name. These are Crest paths, not the locally extracted archive paths.

## Array tasks

| Index | Lane | Haplotype |
|---:|---|---:|
| 0 | Chromosome/anchor correspondence | 1 |
| 1 | Chromosome/anchor correspondence | 2 |
| 2 | Native graph neighborhoods | 1 |
| 3 | Native graph neighborhoods | 2 |
| 4 | Competitive local HiFi | 1 |
| 5 | Competitive local HiFi | 2 |
| 6 | Library-separated Hi-C | 1 |
| 7 | Library-separated Hi-C | 2 |

No lane waits for another. All evidence is read-only. Exact assessment hashes, HiFi BAM dictionaries, Hi-C source/AGP hashes, pairs dictionaries and library manifests are checked in the relevant lanes. The existing packet supplies trusted candidate/control coordinates directly, avoiding stale hand-entered registry hashes. Historical assemblies are not required.

## Returned evidence

Each task writes comparisons/ctlk-adjudication-JOBID/LANE_hapH.review.tar.gz. Return all eight archives, plus failure logs if any. Native reads/sequence files remain on Crest and are excluded from review archives; small graph topology and native A records are retained. status.json reports measured outputs and interpretation limits, never biological absence from a failed job.

Anchors: 10/25 kb flank tiles at multiple distances plus 25 kb tiles every 1 Mb across every assessment scaffold >=5 Mb. Map to the six better haplotypes and both current sisters. Retain all PAFs and chromosome name maps. High MAPQ alone does not override a comparable alternative chromosome/repeat placement. Suffixes and phase numbers do not establish homology. Name-map labels on a composite scaffold are explicitly all retained, rather than silently selecting one label. The emitted-alternative cap means the unique_screen flag is a conservative screen, not a mathematical uniqueness guarantee.

Graphs: map padded original raw-contig cores to native segment sequences in both and HiFi-only outputs; export seed placements and one-hop links/tags/read-placement A records. Missing segment sequences are reported. This is local graph evidence; cross-stage ancestry and uniquely spanning molecular support must be reviewed from placements/topology, not invented from coverage tags. Whole native graphs/caches are preserved.

HiFi: retrieve local records including supplementary/secondary records from the matching BAM, preserve the longest available molecule representation, and flag hard clipping. Map competitively against merged local corresponding regions of all eight current assemblies. Preserve PAF/CIGAR placements, scores and alternatives, exact target coordinate lifts and original flank coverage. No artificial alternative join sequences are invented. Reads absent from regional BAM retrieval are not sampled; some representations can remain clipped. Broad-core zero bridges are not a cut criterion. Reassess narrow boundaries after graph localization, using these retained reads/placements and the same controls.

Hi-C: one additional pairs scan per haplotype yields a compressed 250 kb endpoint matrix, per-library margins, candidate/control flank counts at five scales, competing partner bins, and distance-binned cross/left/right contact profiles. Profiles include library-size and coordinate-pair-area normalization to compare unequal windows and separations. This is geometry normalization, not mappability calibration. The unique-anchor lane subsequently supplies reliable block labels/mappability for joint interpretation without rescanning pairs. Upstream UU filtering is not equivalent to all endpoints being uniquely anchored to a chromosome block. Original scalars and library separation are preserved.

## Joint evaluation

Apply anchor correspondence to chromosome-piece identities and contact bins, inspect localized graph boundaries against HiFi competing placements, and compare candidate profiles with the retained comparison controls. Comparison controls are not pre-certified intact joins; ambiguous controls must be discarded from calibration. Reproducibility across treatments shares the same input HiFi molecules and does not add independent molecular replicates. Six peer haplotypes represent three individuals.

Assign RETAIN, BREAK_PROBABLE_MISJOIN, UNJOIN_UNSUPPORTED or UNRESOLVED with exact coordinates and reasons. A raw contig transition cannot be cut at an arbitrary uncertainty-interval midpoint. Gap unjoining can preserve every base without claiming a biological fusion impossible. Only after those decisions should one isolated correction-validation batch test sequence accounting, chromosome coherence, and recreation of removed joins.

## Validation status

Local Python syntax and whitespace checks passed. Meaningful repeat-ambiguity, graph-neighborhood, SAM-orientation/clipping and distance-exposure unit tests run at the start of each cluster task. Local end-to-end bioinformatics/tool integration has not been run; cluster outcomes must be checked. No jobs were submitted locally and no production defaults or assembly sequence were changed.
