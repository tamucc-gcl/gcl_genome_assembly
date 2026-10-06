# CTlk localized-boundary batch 1511285

## Decision

All six tasks completed successfully. The H01 isolated gap unjoin passes sequence preservation and is suitable for the proposed conservative scaffold correction. Do not cut the five investigated internal cores using these tests: their local sequence paths are supported, and the proposed J02 raw-graph seam has direct spanning molecules.

This establishes local continuity, not biological chromosome fusion. A supported path through homologous/repetitive sequence need not uniquely establish the chromosome-scale placement of the flanking blocks. RETAIN_LOCAL_PATH below means do not cut that tested interval; it does not accept the whole composite scaffold as an authentic fusion.

No production assembly has been modified. The isolated H01 candidate is on Crest at comparisons/ctlk-boundaries-1511285/H01/unjoined.fa; FASTA is deliberately excluded from the returned archive.

## H01 correction validation

Original hap1 scaffold_1: 138,446,979 bp. Candidate pieces:

- scaffold_1__left: 63,051,325 bp, including the original 100 gap Ns at its terminus.
- scaffold_1__right: 75,395,654 bp.

Original and reconstructed scaffold sequence SHA256 both equal 080819efa656a3deddf562ae2286d3248fbe13058e1d775068d12238bd5b20a5. Indexed assembly total remains 953,570,870 bp; record count changes from 362 to 363. Existing anchor/contact evidence is relabeled, not remeasured. The 250 kb contact bin straddling the boundary remains ambiguous rather than being assigned to either side.

Together with the earlier chromosome-block discordance and poor localized two-library contact support, this validates the sequence-preserving UNJOIN_UNSUPPORTED candidate. It is not proof a biological fusion is impossible. The right piece still contains chr7/chr12 material; H01 alone does not finish CTlk correction. Neither YaHS recreation nor final harmonization was run in this batch.

## Narrow molecular evidence

Counts below are distinct molecule names in retained competitive placements with >=1 kb aligned on both sides and MAPQ>=20. These are local mapping support, not globally unique chromosome/phase assertions.

| Core | Passing local probes / total | Minimum / median / maximum passing molecules |
|---|---:|---:|
| J01 | 304 / 324 | 0 / 9 / 37 |
| J02 | 170 / 170 | 5 / 16 / 37 |
| J05 | 138 / 138 | 2 / 11 / 32 |
| J06 | 179 / 179 | 7 / 17 / 63 |
| J07 | 139 / 139 | 5 / 13 / 23 |

Every sampled position in all five cores has at least one spanning placement at some MAPQ. J01 has 20 probes without passing-MAPQ support, but no run without any spanning placement. Do not interpret the high-MAPQ drop as a cut without examining competitors and chromosome-specific anchoring.

The same assay still has passing molecules at only 1/8 hap1 comparison gaps and 4/10 hap2 comparison gaps. Therefore absence across actual scaffold gaps is still not calibrated as a misjoin classifier. Presence at a narrow internal seam is useful affirmative evidence; zero counts remain weak negative evidence.

### J02 seam

Original raw h1tg000133l:4085098-4091627 projects, in reverse orientation, to assessment scaffold_1:[95978246,95984775).

- Lower assessment seam endpoint: 13 passing molecules.
- Upper assessment seam endpoint: 8 passing molecules.
- Entire 6,529 bp seam plus >=1 kb on each side: 3 passing molecules (10 at any MAPQ).

All three full-seam molecules have MAPQ60 and near-complete whole-read query placements:

| Molecule suffix | Assessment span | Query aligned bases | Alignment identity |
|---|---|---:|---:|
| 131598796/ccs | 95965240-95988679 | 23,436 | 99.94% |
| 261423625/ccs | 95968932-95996042 | 27,098 | 99.74% |
| 44630382/ccs | 95974074-95989355 | 15,279 | 99.92% |

Full names begin m84066_260114_212201_s2/. The raw graph mapping seam does not identify an unsupported physical join. Reject this interval as a cutting hypothesis on the current evidence.

### J07 raw-segment endpoints

The two graph-endpoint tests at assessment positions 39685322 and 39667254 have 13 and 15 passing molecules respectively. Neither is justified as a cut merely because a raw graph segment ends there.

### J05 and J06

No raw graph endpoints occur inside their tested cores. All sampled narrow positions have passing spanning molecules. Midpoint cutting of these cores is unsupported by the evidence.

## Smaller anchors

The smaller-anchor assay improves sensitivity but shows mixed/interleaved chromosome labels instead of a single clean, uniquely established boundary. Counts include overlapping tile sizes and haplotypes; they are not independent observations.

- J01 vicinity has chr9 and chr7 placements over overlapping coordinate ranges; the broad scaffold blocks are real, but the internal core is not localized as a clean chromosome-specific junction.
- J02 vicinity has chr7 labels from approximately 95.14-96.40 Mb and chr12 from 96.05-96.56 Mb, plus isolated other labels.
- J05 has chr14 correspondence supported across three peer individuals downstream; the few local chr4 placements in this dense assay are from one peer individual. Earlier broad tiles establish chr4 elsewhere on the scaffold, not a precise breakpoint here.
- J06 has overlapping chr12 and chr7 ranges and isolated placements to other chromosomes.
- J07 has chr14 and chr4 correspondence across all three peer individuals, but interleaving and outlier short-anchor labels prevent a simple nucleotide cut call.

A MAPQ-based unique_screen against one reference genome is not sufficient to prove that a 2 kb sequence identifies one chromosome across individuals. Evaluate long collinear blocks and competing homologous copies; do not aggregate isolated short hits into a forced chromosome switch. The current data do not support declaring these local paths errors solely from chromosome labels.

## Resulting action worksheet

H01: UNJOIN_UNSUPPORTED candidate validated; promotion remains an isolated assembly-correction step with downstream validation.
J01: retain tested local path; use the adjacent H01 gap for the proposed unjoin. Chromosome-scale structure beyond that remains unresolved.
J02: retain tested local path; explicitly reject the proposed 6,529 bp seam as an unsupported cut.
J05/J06/J07: retain tested local paths; no internal cut is justified at the tested cores or J07 endpoints. Chromosome-scale composites remain unresolved rather than automatically approved as authentic fusions.

The next implementation decision should preserve these distinctions in the pipeline: detect chromosome-scale discordance; identify removable scaffold adjacencies separately; and veto unsupported internal cuts where informative spanning molecules support the tested local path. Do not convert supported local sequence continuity into automatic retention of an entire chromosome-scale fusion.

No further broad bridging or library/parameter sweep is justified by this packet. The retained peer placements and molecule competitors are sufficient for focused correspondence review; any additional job must test a specific alternative chromosome path/placement rather than resampling the same cores.

## Reproducibility

Extracted results: comparisons/local-ctlk-boundaries/comparisons/ctlk-boundaries-1511285/.
Local audits: comparisons/local-ctlk-boundaries/audit.py and read_check.py.
Evidence limits: BAM-ascertained reads and emitted competitive alternatives; regional references do not enumerate every possible path; native graph placements do not prove global chromosome uniqueness. No final biological-fusion determination or production cuts were performed.
