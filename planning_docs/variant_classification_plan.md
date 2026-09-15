# Pangenome variant classification + SV/private-sequence rework

**Status:** batches 1–5 built, applied and verified running end to end on the rebuilt graph.
Every process in the pangenome arm executes and produces output. Batches 6, 6b, 7 and 8
outstanding. Batch 6 C1 has RUN and is falsified — see §4b; `--lastTrain` is kept.

**Origin:** Chris Bird, 2026-08-26 — bp-weighted SV spectra, private-haplotype size spectra,
independent mapping of private haplotypes, transposon drivers.

**Cohort:** 5 individuals × 2 haplotypes = 10 paths. Reference `Sde-CMat_203_hap2`.
Passengers `Sde-CPla_115_hap1/hap2` retained. Clip graph 1,926,892,905 bp / 151,697,494
nodes; full graph 2,631,391,946 bp / 153,003,433 nodes.

**The graph was rebuilt after batch 1.** Current clip total is **1,926,884,214 bp**
(odgi-confirmed), 8,691 bp smaller than the pre-rebuild figure, and per-haplotype private bp
moved bidirectionally by ≤0.1% — largest gains in `CBau_104#2`, `CTlk_104#2` and
`CPla_115#1`, the assemblies whose part indices the harmonization patch changed. Every
pre-rebuild number in §1 is therefore a **baseline for comparison, not a value to reuse.**

---

## 0. Findings the rework surfaces

The rework has to make candidate results visible in published outputs without hand analysis.
§0.1 is both the first candidate and the acceptance test.

### 0.1 Candidate segregating inversion, chr10 — RE-DERIVE ON THE NEW GRAPH
From `odgi untangle` on `chr10_1.full.og` against `Sde-CMat_203_hap2#0#chr10_1`, **pre-patch
graph**.

Two loci, ~1.68–2.67 Mb and ~58.2–70.2 Mb (terminal ~12 Mb of ~71 Mb), carried by
`Sde-CBau_104#2` and `Sde-CTlk_104#2` only.

**Evidence.** 2 of 10 haplotypes; the other 38 chromosome-scale queries carry nothing,
**including both carriers' sister haplotypes** — so both individuals are heterozygous, AF 0.2.
Carriers reach the same reference interval from *different* query scaffolds. Global
orientation is not the explanation: the three `rev` chr10 scaffolds project 3.8%, 15.3% and
45.9% inverted, where a failed reverse-complement reads ~100%.

**Structure is stable across a 100× resolution change**, which is the strongest argument.
Run span (gaps <200 kb merged) 14.45 / 14.48 / 14.49 Mb at `-e` 1 Mb / 100 kb / 10 kb, while
the union of inverted intervals tightens 13.95 → 13.74 → 12.85 Mb. Fixed region, sharpening
edges. Segmentation artifacts do not behave that way.

**Must be re-derived from scratch:** `Sde-CBau_104_hap2`'s chr10 pieces **swapped names**
(`chr10_1 ↔ chr10_2`) under the batch 1 patch, so piece assignment changed, not just labels.

### 0.2 Two composite scaffolds are chimeric joins, not orientation errors
| query | fwd | inv | %inv | harmonization |
|---|---|---|---|---|
| `Sde-CPla_115#1#chr4_7+chr10_10` | 10 bp | 715,908 | 100.0 | `fusion;unsupported(ref4+ref10:0f/8s)` |
| `Sde-CPla_115#1#chr10_17+chr11_5` | 0 bp | 1,415,923 | 100.0 | `fusion;unsupported(ref10+ref11:0f/8s)` |

Both were already `fwd` and did **not** flip under the patch, so 100%-inverted projection is
not an orientation-call failure. `0f/8s` = fused in 0 voters, split in 8: a chimeric join,
which no orientation logic can fix because a scaffold whose halves need opposite orientations
has no correct single value. Both flag spellings exist — `unsupported(...)` for no voter
support, `chimera_suspect(...:1f/7s)` for a fused minority.

### 0.3 The reference "private-sequence excess" was a clipping artifact — RESOLVED
| | clip | full |
|---|---|---|
| `Sde-CMat_203_hap2` (reference) | 122,550,615 = **15.08%, highest of 10** | 122,010,508 = **8.08%, lowest of 10** |
| `Sde-CTlk_104#2` | 83,521,123 | 203,673,145 (+144%) |
| `Sde-CBau_104#1` | 88,779,495 | 163,550,261 (+84%) |

The reference is unchanged between flavors (540 kb of 122 Mb) while every other haplotype
gains 70–144%. The reference is the graph backbone and is never clipped; everything else loses
its unaligned sequence. Below even share in the full graph.

**Consequence:** the 5,913 fixed-alternative `SV_DEL` alleles read as reference-specific
insertions are clip-derived and suspect for the same reason. Re-derive on the full graph.

### 0.4 The real private-sequence signal is chromosome-level
21 contigs exceed 5 Mb and 20% private. **19 of 21 are plain `chromosome` class**, not
composites — the chimera hypothesis is refuted. The 20–27% band is almost entirely `chr8_1`
and `chr9_1`, appearing across seven and four haplotypes respectively, and the clip graph
agrees from the reference's own side (chr8_1 19.9%, chr7_1 17.4% against a ~9.5% baseline).
Ten haplotypes, two flavors, same chromosomes elevated: genome biology, not assembly artifact.

Separately, three 5–8.5 Mb **secondary** pieces run 60%+ private — `Sde-CTlk_104_hap2 chr3_2`
63.3%, `chr7_2` 61.9%, `Sde-CPla_115_hap1 chr8_5` 60.3% — with no fusion flags. Those are the
QC candidates.

**Nothing currently emits a per-chromosome private table.** It should — this is the strongest
signal in the data and is invisible in the per-haplotype and per-contig views.

### 0.5 Requirements on REARRANGE (implemented)
Artifact discrimination joins `fusion` / `unsupported(...)` / `chimera_suspect(...)`, with
per-query %inverted secondary. Candidate table sorted **unflagged first, then fewest
carriers**, so a clean 2-of-N locus outranks a flagged 1-carrier artifact. Both footprint
measures reported (union and run span; fill is a quality signal). `pct_inv` and
`pct_inv_unfiltered` both emitted, and `ORIENTATION_SUSPECT` needs both high, so a `-j` filter
alone cannot raise the flag.

---

## 1. Measured facts

### The vcfbub cap is resolution recovery, not a ceiling
| invocation | records |
|---|---|
| `--max-level 0 --max-ref-length 100000` | 31,365,160 |
| `--max-level 0` (no cap) | 23,877,797 |
| `--max-ref-length 0` | **0** |

The cap **adds** 7.49M records by popping oversized parents into children. `0` is a literal
zero limit that empties the catalog. 542 LV=0 records carry 234 Mb above it. Stays at 100000.

### The parent tier is exactly LV == 0
23,877,797 — identical to `vcfbub --max-level 0`. No span floor OR-ed in; vcfbub's hierarchy
handles nesting and OR-ing one would double-count a nested bubble inside its parent.

### `AT` is universal and only valid pre-decomposition
41,920,948 / 41,920,948 raw records, 0 unparseable. vcfwave and `bcftools norm -m` rewrite
REF/ALT while `AT` is inherited, so allele *i* stops matching traversal *i+1*. Enforced by a
hard guard, not by trusting the wiring.

### Topology classification, measured on the real parent view
`records=23,877,797  alt_alleles=31,314,795 (1.31/record)`

| primary_class | alleles | per-allele bp | pangenome node bp | novel node bp | merged ref footprint |
|---|---|---|---|---|---|
| SNP | 17,741,162 | 17,741,162 | 527,452,057 | 17,741,162 | 15,686,744 |
| INDEL | 11,722,380 | 128,706,035 | 251,726,496 | 27,553,309 | 27,206,020 |
| SUBST | 1,133,177 | 2,507,120,835 | 911,915,764 | 485,046,813 | 89,538,307 |
| INS | 469,008 | 321,423,391 | 261,710,853 | **246,632,865** | **5,196,150** |
| DEL | 249,021 | 105,896,522 | 100,215,970 | **0** | 57,809,651 |
| INV_PATH_EXPLICIT | 27 | 3,102,759 | 3,243,938 | 281,980 | 2,960,523 |
| DUP | 18 | 547,316 | 417,644 | 171,832 | 244,622 |
| INV_DUP | 2 | 8,330,920 | 4,934,127 | 768,402 | 4,165,460 |

Topology reclassifies 233,015 length-called `SV_INS` and 178,583 `SV_DEL` as SUBST.

### Three bp measures, because one number cannot answer the question
- **per-allele SV sum: 2,946,421,743 bp on a ~1 Gb reference** — 2.9× the genome. SUBST alone
  is 28× its merged footprint. Multiallelic sites reuse one reference span, and
  `max(REF,ALT)` counts ALT length when ALT is longer. Correct as a size-spectrum y axis,
  meaningless as a total.
- **merged reference footprint** is reference-*biased*: INS 5.2 Mb against DEL 57.8 Mb,
  because an insertion barely touches the reference.
- **novel node bp** is graph-native: INS 246.6 Mb, **DEL exactly 0** — the correctness proof,
  since a deletion traverses only reference nodes and cannot contribute novel sequence.

Graph totals reconcile: `graph_total_bp = 1,926,892,905` matching `odgi stats -S` exactly;
`pangenome_node_bp` 1,788,647,791 (92.8% of the graph in some bubble traversal);
`novel_node_bp` 770,851,175 against ~912 Mb of non-reference sequence (85%, remainder in
nested bubbles only) and 813 Mb private of which 123 Mb is reference-private → 690 Mb
non-reference private, leaving ~222 Mb shared among 2+ non-reference haplotypes. Three
independent derivations in a consistent frame.

### Inversions are overwhelmingly not bubble-representable
| detector | result |
|---|---|
| `AT` path-explicit | 29 alleles at LV=0 (32 across all levels), ~11 Mb |
| revcomp alignment on SUBST | 411 alleles, 2.0 Mb minus-strand, **all <100 kb, zero above 500 kb** |
| `odgi untangle`, chr10 alone | ~21.5 Mb inverted, 5–7 blocks |

Rescue rate is **flat at 0.11–0.20% from 1 kb to 100 kb**, then 0.02% at 100–500 kb and zero
above — so the calls follow SUBST's own size distribution rather than being size-skewed, and
large SUBST is definitively not inversion. Overall 411/309,879 = **0.13%**, so 99.87% of SUBST
≥1 kb is genuine allele replacement.

**Positive control: 30 of 32 path-explicit inversions recovered by the alignment test**, so
`--min-frac 0.8` is calibrated. A threshold sweep adds 425 → 719 calls from 0.8 → 0.3 but only
2.13 → 2.95 Mb, so nothing large is hiding below the threshold.

Romain et al. (bioRxiv 2025.03.14.643331) supply the framework: path-explicit inversions
traverse shared nodes in opposite directions; alignment-rescued ones are disjoint paths.
`INV_PATH_EXPLICIT` is therefore a **floor, not an inversion count**, and must be labelled so.

### SUBST subdivided by measured homology — replaces the "unresolved tier"
Of SUBST alleles ≥100 kb: **4,015 of 4,329 align to their reference allele at ≥0.5 with no
inversion.** Large SUBST is homologous sequence the *graph* fragmented, not an alignment hole.
And "no homology either way" is a **small**-allele phenomenon — 22.9% at 1–10 kb, 11.3% at
10–100 kb, 3.3% at 100–500 kb, 1.0% at 500 kb–1 Mb, 3.5% above 1 Mb.

So a span/ratio threshold would have routed almost exactly the wrong records.
`pangenome_unresolved_span` and `pangenome_unresolved_ratio` were **deleted** in favour of
`SUBST_HOMOLOGOUS` / `SUBST_PARTIAL` / `SUBST_UNRELATED` from `max(frac_fwd, frac_rev)`.
Nothing is excluded from the catalog.

### Allele frequencies from AC/AN, not sample columns
`AN = 9` (four diploid columns + one haploid; reference absent). **`SV_DEL` is the only class
with a real frequency spectrum** — 62% singleton, then 45,317 / 20,552 / 12,227 / 8,658 /
6,935 / 6,772 / 6,815 / 5,913 across AC 2–9. Every other class collapses: `SV_INS` 92%
singleton, `SV_COMPLEX` 90%, `SV_BLOCKSUB` 96%.

**AC=9 is fixed alternative: 5,913 loci where the reference carries sequence no other
haplotype has.** That revises down from the 30,224 I got by counting sample columns, which
conflated diploid with haploid. That only the deletion class shows a spectrum is the
reference-bias signature — but see §0.3, it is clip-derived.

### Q8: the two private-sequence derivations agree
Normalising by path count (diploid ÷2, haploid ÷1): CBau 2.35 M alleles / 86.5 Mb ·
CMat_203_hap1 2.20 / 80 · CLim 2.16 / 78 · CTlk 2.10 / 82 · CPla 1.64 / 59. Same first, same
last, middle three within ~5%. Caveat: CPla's 7.74% missingness depresses its allele count
independently.

### Private segments
| | segments | bp | mean |
|---|---|---|---|
| clip, <1 kb | 33,573,074 | 137,878,974 | **4.1 bp** |
| clip, ≥1 kb | 163,418 | 675,062,257 | 4.1 kb |
| full, <1 kb | 33,920,277 | 140,397,398 | 4.1 bp |
| full, ≥1 kb | 170,311 | 1,371,472,664 | **8.1 kb** |

**99.5% of segments carry 17% of the bp**, so a count histogram is a single spike and the
bp-weighted one is the answer — exactly Chris's point. The sub-1 kb population is identical
across flavors; all 698 Mb of the clip/full difference lands in ≥1 kb, gaining only 6,893
segments, so the recovered sequence is in a small number of very large segments.

Regression: clip `private_bp = 812,601,533`, byte-identical to the pre-patch run, and
`reconstructed_graph_bp` matches the S-line total, so the fifth pass did not perturb the four
verified passes. `segment_bp_all − private_bp = repeat_traversed_bp` exactly (339,698 clip /
908,093 full).

### Clipping removes almost nothing but private sequence
704,499,041 bp removed, of which **698,360,436 (99.1%) was private**. The clip graph
understates private content by 46%. Clipping also cuts paths into subpaths, so the clip graph
carries 3,077 *more* paths despite less sequence — which is why its `.og` is 118 GB against
the full graph's 32 GB, and why untangle runs primarily on full.

### Label co-occurrence is nearly nested — upset plot dropped
Only **7** combinations exist; three cover 99.995%: `DEL,NOVEL_INS` 1,133,177 · `NOVEL_INS`
469,008 · `DEL` 249,021. Every `INV_PATH_EXPLICIT` co-occurs with `DEL,NOVEL_INS`; `DUP` never
appears without `DUP_NOVEL`. `primary_class` discards essentially nothing. Keep the
co-occurrence table (7 rows, proves the point), drop the figure.

### odgi v0.9.2 in cactus v3.1.4 has a broken internal step-index builder
`sdsl::int_vector::operator[]` assertion during "Collecting Steps". Not ID space (`odgi sort
-O` does not fix it); `odgi validate` passes. Standalone `odgi stepindex` succeeds and
`odgi untangle -a` then returns 0. Hence `PANGENOME_STEPINDEX` as its own process — which is
also worth 26% of the cost each time an untangle parameter is retuned. `odgi optimize` does
not exist; it is `odgi sort -O`. This `vg` build has **no `-e` flag** (relevant to Swave,
whose README uses `vg deconstruct -e -a`).

### Untangle: measured, and all queries beat filtered queries
chr1_1 (94.7 Mb, 4 threads): stepindex 3m42s / 6.9 GB; untangle 10m26s / 7.7 GB at 293% CPU
(serial phases inside untangle, so >8 threads buys little).

Cost is **flat in `-e`**: 6m29s / 6m45s / 6m19s at 1 Mb / 100 kb / 10 kb, memory constant.
Finer is also more *accurate*, so `-e 10000`: 98.8% of inverted bp survives a 1 kb floor,
median segment 3.8 kb, breakpoints localised to ±10.8 kb rather than ±603 kb.

**All 394 queries ran in 3m35s against 6m19s for 43** — unrestricted is faster, because cost
is dominated by chromosome-scale paths. And the 346 unplaced scaffolds project **coherently**:
bp/span 0.81–1.00 across 50.7 Mb, making untangle a scaffold-placement instrument for free.
Four queries project nowhere above `-j 0.1` — sequence with a chromosome assignment and no
reference homology, i.e. private sequence with a location, and worth reporting explicitly.

Per-chromosome untangle **cannot see inter-chromosomal translocations** (cactus splits by
refContig, so `ref.name` has one value). Composites *are* visible, giving a free cross-check
against fusion flags. Translocations need one whole-graph untangle with all reference paths —
one bigmem job, costed separately.

### Batch 1 harmonization: correct, and near-null on this cohort
| check | baseline | patched |
|---|---|---|
| rows | 7687 | 7687 |
| class counts | 540 / 139 / 7008 | **identical** |
| `inflated_aln` | 645 | **645** |
| `overlap` | 371 | 340 (75 lost, 44 gained) |
| **orientation flips** | — | **0** |

Both controls held. Site 1071 (orientation) is a null result here. Overlap churn is
bidirectional as designed — the 44 gains are non-adjacent overlaps the old running-envelope
walk structurally could not see — and 95 of 96 changes are in the passengers. Total voter
effect: one Tlk overlap flag, six Tlk part renames, one Bau `chr10_1 ↔ chr10_2` swap.

38 chromosome rows are non-monotone in part index versus reported `ref_span` (36 passengers).
Expected, not a bug: the sort key is `footprint_start` (largest merged block) while `ref_span`
reports `fp[0][0]` (leftmost). **The sort key is not emitted, so ordering is not auditable** —
still outstanding.

---

## 2. Architecture as built

| Layer | Instrument | Delivers |
|---|---|---|
| 1 | `odgi untangle`, per chromosome, full primary | inversions, dispersed duplications, rearrangement, scaffold placement. `self.cov > 1` is the **only** duplication signal anywhere (AT found 27 node re-visits in 3,268,312 SV alleles). |
| 2 | `AT` traversals, parent tier only | SUBST / INS / DEL partition replacing `COMPLEX`/`BLOCKSUB`; `INV_PATH_EXPLICIT` as a named floor |
| 2b | revcomp + minimap2 on SUBST ≥1 kb | `INV_ALN_RESCUED` + the SUBST homology subdivision |
| 3 | `AC`/`AN` | AF spectrum, per-individual private variants |
| 4 | `gfa_hap_coverage.py` pass 5, **both flavours** | private-segment spectra and BED |
| 4b | `PRIVATE_FASTA` → `MAP` / `KMER` → `JOIN` | per-segment evidence: does the sequence exist elsewhere, and is it repeat-derived |

Exclusive `primary_class` for all totals; non-exclusive labels retained per variant.

**No chromosome scatter for `CLASSIFY` or `INVERSION_RESCUE`** — 23.9M records classify in
well under an hour single-threaded, and rescue completed as one task. `HAP_COVERAGE` is
whole-graph because 139 composites span multiple chromosome subgraphs and a private run cut at
a boundary would corrupt the segment histogram.

### Layer 4b: the private-sequence chain, and why it is shaped this way

| process | tasks | scatter |
|---|---|---|
| `PRIVATE_FASTA` | 2 | per flavour — GFA read **once**; emits one private and one control FASTA per haplotype |
| `PRIVATE_INDEX` | 1 | flavour-independent; tagged multi-FASTA of all assemblies + **one-part** minimap2 index |
| `PRIVATE_MAP` | 20 | per haplotype × flavour; **both sets in one task** |
| `PRIVATE_KMER` | 20 | per haplotype × flavour; **both sets in one task, control first** |
| `PRIVATE_JOIN` | 2 | per flavour — joined table, cross-tab, and the R-ready CSV |

**Two independent verdicts, by design.** `PRIVATE_MAP` aligns each segment against the other
assemblies; `PRIVATE_KMER` counts its k-mers in the sample's own reads. Neither sees the
other's evidence, so their cross-tabulation is informative rather than two views of one
measurement:

| combined | reading |
|---|---|
| `NOT_PRIVATE` + `REPEAT_LIKE` | present elsewhere AND high copy → graph collapse |
| `PRIVATE_CONFIRMED` + `UNIQUE_LIKE` | absent elsewhere AND single copy → novel sequence |
| `NOT_PRIVATE` + `UNIQUE_LIKE` | present elsewhere, single copy → the graph failed to merge homologous sequence; an alignment failure, not a repeat problem |
| `PRIVATE_CONFIRMED` + `REPEAT_LIKE` | absent elsewhere, high copy WITHIN this sample → haplotype-specific expansion |
| `NO_ALIGNMENT` + * | aligned nowhere, not even to its own assembly — low-complexity, enumerated from the FASTA rather than the PAF so it is distinguishable from a failed task |

### The control set is the load-bearing idea

Private measurement alone cannot say whether private sequence is repeat-enriched. That needs
sequence from the **same haplotype**, measured the **same way**, differing only in privateness.
`PRIVATE_FASTA` therefore emits size-matched **non-private control windows** per
(haplotype, contig), and every threshold that could be hand-set is derived from them instead.

Three requirements on a control window, each from a measured failure:

1. **No private content** (`max_private_frac = 0.05`). ~14% of a haplotype is private, so
   unfiltered windows would be ~14% contaminated toward the private value. Rejected rather
   than masked: masking creates junction k-mers that exist nowhere in the genome.
2. **Size-matched per (haplotype, contig)**, because multiplicity correlates with length and
   private sequence is not uniformly distributed across chromosomes. A genome-wide bp match
   would confound both with privateness.
3. **Cross-individual** (`min_cross_frac = 0.95`): 95% of bp on nodes walked by ≥2
   INDIVIDUALS, not merely ≥2 haplotypes. `cov ≥ 2` is satisfied by a window's own sister
   haplotype, which is not the contrast `PRIVATE_MAP` performs.

### Derived, not guessed

| quantity | source | measured |
|---|---|---|
| k-mer single-copy reference | **control** run, k-mer-weighted median of per-segment medians | 22–65× per haplotype, tracking read depth |
| `min_identity` | **control** p5 of per-segment best identity | 0.634–0.674 across ten haplotypes |
| `min_frac` | **constant 0.5** — deliberately not derived | see below |

Each of these was a hand-set constant that was wrong:

- k-mer threshold `3.0` as an **absolute** multiplicity, against ~14× single-copy coverage,
  called everything `REPEAT_LIKE`. It is now a **multiple** of the control's derived level.
- `min_identity = 0.90` discarded nearly all real homology (control identities run 0.75–0.90),
  so 70% of control windows reported no homologue anywhere. `0.80` fixed this cohort but was
  still a number from one clupeid.
- `min_frac` derived as a control p5 gave **0.0818 → 0.9414** across ten haplotypes of one
  species — noise in the long left tail of coverage, not biology, and it made the criterion
  differ per haplotype. Reverted to a constant.

**The distinction that matters:** identity measures cross-haplotype divergence, a species
property that does not transfer between taxa, so it must be derived. Coverage is a definitional
choice — how much of a segment must be found elsewhere before it stops counting as private —
and is the same choice for any taxon. `aligned_frac_merged` is strongly bimodal (segments pile
near 0.00 or 1.00), so anything from ~0.3 to ~0.7 gives nearly identical verdicts; 0.5 sits in
the empty middle.

**Every derivation carries its own tripwire.** Computing the identity percentile requires
counting control segments with no non-self hit, and `max_control_no_hit = 0.5` fails the task
when that is implausible. Measured 0.000–0.010 once the index was single-part — it is the check
that would have caught both the split minimap2 index and the 0.90 floor on the first run.

---

## 3. Built and applied

**Batch 1.** `harmonize_names.py` four envelope → merged-footprint swaps.

**Batch 2.** `cactus_pangenome.nf` (`--vcf full clip` + seven full-graph emits) ·
`pangenome_variants.nf` (parent tier + tier audit; awk classifier removed) ·
`pangenome_hap_coverage.nf` + `gfa_hap_coverage.py` (pass 5) · `pangenome_plots.R` ·
five new modules · `classify_variants.py`, `rescue_inversions.py`,
`rearrange_from_untangle.py` · wiring across `main.nf` / `pangenome.nf` /
`harmonize_scaffolds.nf` / `nextflow.config`.

**Batch 3.** `pangenome_hap_coverage.nf` flavour-parallelised · five private modules ·
`extract_private_fasta.py` (4-pass + control windows), `summarise_private_map.py`,
`summarise_private_kmer.py`, `join_private_evidence.py` · `main.nf` passes
`BUILD_MERYL_DB.out.meryl_db` · config params and five resource labels.

**Batch 4.** `rearrange_from_untangle.py` (carrier split into chromosome-scale vs unplaced,
re-sorted on chromosome-scale support, `n_chrom_individuals` added) ·
`PANGENOME_PRIVATE_PLOTS` + `pangenome_private_plots.R` (8 figures, BOTH arms in one task) ·
`PANGENOME_REARRANGE_PLOTS` + `pangenome_rearrange_plots.R` (5 figures, flavour-parallel) ·
`pangenome_plots.R` trimmed of the two private figures and extended with
`ref_footprint_by_chrom` + `length_class_summary` · wiring, and the `ch_hap_priv` double-read
resolved by REMOVAL rather than by forking.

**Batch 5.** `pangenome_report.R` gains the "Which view says what" matrix and per-section view
labels · `pangenome_report.nf` inputs 9 → 15 · wiring for the six matrix channels · figure
resource labels tightened to measured values.

### Verified on real data

| claim | evidence |
|---|---|
| chromosome-blind footprint fixed | `ref_footprint_by_chrom.tsv` has **15** chromosomes; `merged_ref_footprint_all_classes` 493,058,268 bp vs the buggy 90,307,914; no per-chromosome row exceeds its chromosome |
| SUBST was the worst-suppressed class | 89.5 Mb → **420.2 Mb** (4.7×), because it has the most intervals spread across all fifteen chromosomes |
| the reference-bias case, stated cleanly | INS: 246,819,137 novel node bp / 5,366,604 ref footprint. DEL: **exactly 0** novel / 94,471,353 ref. 46× one way, zero the other |
| `graph_total_bp` matches odgi | 1,926,884,214 exactly, independently derived from S lines |
| decomposition guard works | fine tier `topology_enabled False`, all four markers detected (`ID=ORIGIN`, `bcftools_normCommand`, `bcftools_normVersion`, `vcfwave`) |
| the two tiers are what they claim | parent 23,877,797 records at 1.31 alleles each; fine 60,453,457 at exactly 1.0000 |
| inversion rescue is calibrated | control 27/29 recovered (93.1%), against 30/32 pre-rebuild. `control_path_explicit_tested = 29` matches `classify`'s LV=0 path-explicit count exactly — two code paths, same population |
| SUBST subdivision reproduces its motivation | `SUBST_HOMOLOGOUS` 153,498 (75.5%) / `PARTIAL` 13,195 / `UNRELATED` 36,403 |
| compaction fixed the odgi crash | `UNTANGLE` 30/30 exit 0 including all fifteen clip tasks |
| `REARRANGE` runs both arms | 1,516,080 / 1,515,536 rows kept, 7,687 harmonization rows joined |
| private extraction is self-consistent | 170,322 segments / 1,371,790,691 bp; `segment_bp == fasta_bp_total`; 0 nodes missing sequence |
| the k-mer contrast | private `repeat_like` **0.9306** vs control **0.0669** on the full arm |
| flat across chromosomes | private `repeat_like` 0.92–0.94 on all fifteen; chr8 0.9438 and chr9 0.9198 unremarkable against chr1 0.9424 |
| measured resources | `PRIVATE_INDEX` peak RSS 28.4 GB (guess was 96 GB), 115 s wall for a one-part index over 10.4 Gb |
| the reference-rank inversion, as a number | `reference_rank_clip` **1 of 10**, `reference_rank_full` **10 of 10**. Most private haplotype on one arm, least on the other. This is why the private figures had to leave `PANGENOME_PLOTS`, which is pinned to clip |
| rearrangement is arm-INSENSITIVE | 665,401,148 bp inverted on clip vs 665,572,453 on full — **0.026%**. Two distinct audit files, so the wiring is right. Subpath fragmentation from clipping is real (556 vs 394 on chr10) but costs no measurable inverted sequence |
| the control is measured, not assumed | derived `min_identity` 0.634–0.674 across ten haplotypes; control no-hit fraction 0.000–0.010; `control_bp_ratio_achieved` 1.009–1.019 |
| the k-mer contrast survives every threshold change | private `repeat_like` 0.9306 vs control 0.0669, unchanged across three rounds of threshold work, because it does not depend on cross-assembly alignment |
| 18 figures, zero skips | 8 private + 5 rearrange × 2 arms; `candidate_rows_unplottable` and `duplication_rows_unplottable` both 0 |
| figure cost, measured | `PRIVATE_PLOTS` 18.5 s / 730.3 MB reading both arms incl. two ~400k-row CSVs; `REARRANGE_PLOTS` 4.7 s / 225.2 MB |

### Bug classes that cost the most time — record these as conventions

1. **A plain channel read more than once starves all but one consumer, silently.** Process
   outputs (`X.out.y`) are broadcast and safe; anything built by `map`/`flatMap`/`combine` is
   not. Three instances: `ch_cactus_in` (five consumers), `CACTUS_PANGENOME.out.gfa` into
   `CLASSIFY` (lockstep against a two-item channel → the parent tier never ran), and the
   `parents_vcf` double read. **The symptom is always a missing task, never an error.**
2. **For a missing Nextflow task, instrument before theorising.** `.view()` shows items
   flowing; `count()` only emits on channel close. Together they distinguish empty / flowing
   but never closing / fine in one run. Five wrong diagnoses preceded two minutes of
   instrumentation.
3. **Read the log, not the progress display.** `[-] process > X -` in the terminal summary does
   **not** mean a process did not run — cached tasks completing during DAG resolution can
   appear as never-started. `grep "Cached process\|Submitted process" .nextflow.log` is
   authoritative. This single misreading drove roughly half a day.
4. **SIGPIPE under `set -o pipefail`.** Any `cmd | head`/`| tail` over a data file is exit 141
   waiting to happen. Rewrite as a single `awk`. `--version | head -n1` is safe. Documented
   once in `pangenome_untangle.nf`, then reintroduced into three other modules.
5. **`$` in a `"""…"""` script block.** Shell needs `\$`; Groovy interpolation does not.
   Post-condition checks should scan the script block for unescaped `$`.
6. **Post-condition strings must be fully qualified.** An unqualified `min_frac     = 0.5`
   matched the unrelated `harmonize_dropoff_min_frac` and blocked a correct patch.
7. **Destructive commands ship with their variable definitions inline.** An unset `$W` turned
   a delete loop into 80 refusals; twice.
8. **Rewrite a file rather than patch it a fourth time.** Three successive anchored patches to
   `pangenome_private_join.nf` stacked into overlapping copies of the same guard block.
9. **`comment.char = "#"` corrupts PanSN names.** `#` is legitimate inside a haplotype key
   (`Sde-CBau_104#1`), so `read.delim(comment.char = "#")` truncates every key at the
   separator and shifts every subsequent column left. It made a per-chromosome figure report
   "no chr* contigs" from a file containing 602 of them. Strip `#` LINES explicitly and keep
   `comment.char = ""`. `pangenome_report.R` already documented this; the knowledge existed in
   the repo and was reintroduced anyway — **twice**, in two new scripts.
   The mirror-image trap: `read_kv` left `comment.char` at `""` and so treated a leading `#`
   line as the header, returned `ncol < 2`, and yielded nothing. Both failure modes are the
   same root cause — comment handling and `#`-in-data cannot be the same mechanism.
10. **`read.delim` does not guarantee numeric types.** One stray non-numeric token turns a
    whole column character and `cut()` dies with `'x' must be numeric` at PLOT time, not read
    time. Coerce the columns you do arithmetic on, explicitly, on read.
11. **A `log10` scale silently discards values ≤ 0.** ggplot warns ("Removed 505 rows") rather
    than failing. Filter explicitly and record the count, or a data problem hides as a warning.
12. **`facet_wrap` on an empty variable aborts from inside `ggsave`** with a stack trace rather
    than a diagnosis. Guard every faceted plot on the variable actually having values, AFTER
    all filtering.
13. **Post-condition FORBIDDEN strings need BLOCK scope, not substring distinctiveness.** Three
    correct patches were blocked by loose anchors: `min_frac     = 0.5` matched
    `harmonize_dropoff_min_frac`, `NO_HAP_PRIVATE` matched `PANGENOME_REPORT`'s own legitimate
    placeholder, and a bare `memory = { ... 16.GB : 48.GB }` matched four labels.
14. **A spanning anchor fails wholesale on any local divergence.** One anchor covering two
    config label blocks could not match once one of them had been edited by hand. Per-block
    edits tolerate that; and a label three patches have touched should be anchored narrowly.
15. **Assert `count(old) == 1` in scratch edits too, not just in delivered patches.** Four
    successive edits to `pangenome_private_plots.R` silently no-opped before I noticed, because
    the scratch `str.replace` calls had no assertion — the exact discipline the delivered
    patches enforce.

---

## 4. Remaining

| Batch | Contents |
|---|---|
| ~~4~~ | **DONE** — see §3 and §4a. 18 figures, zero skips. |
| ~~5~~ | **DONE** — the "Which view says what" matrix, three views × fifteen measures, every unavailable cell carrying its reason. |
| 6 | Pangenome construction: `--lastTrain` (v3.1.4, available now) then `--gref` (needs v3.2.1). See §4b — two runs, not one, so the scoring change is not confounded with a version bump. |
| 6b | **`PRIVATE_ENRICHMENT`**: the GLMM, as its own module so the model can be re-fit without redoing the k-mer lookups. See §4a2. Waiting on the R code. |
| 7 | Swave as a locus-level direction annotation, EXPLORATORY. See §4c — it merges where we decompose, so the integration is an annotation layer, not a replacement classifier. |
| 8 | **Does the full arm earn its keep?** Evaluate the clip/full differences now that they are measured, and reconsider why the pipeline carries both when most published minigraph-cactus work uses clip alone. See §4d. Deliberately AFTER everything works. |
| later | Whole-graph untangle for translocations. `ref_span` emitting its sort key. GraffiTE post-annotation (the meryl k-mer proxy is a permanent self-contained feature, **not** a placeholder for it). `svim-asm` as the non-graph check. |

**Parked, unevaluated:** `INVPG-annot` (2025 preprint). PGGE / `peanut`.

---

## 4a. Batch 4 — DONE. What the real REARRANGE output showed

`REARRANGE` ran and its candidate table is dominated by artifacts of how it sorts and pools:

**Unplaced scaffolds swamp the carrier lists.** The chr10 locus at 58,005,273–68,935,373
(span 10.93 Mb, union 9.78 Mb, fill 0.8945) has **80 carriers, of which 74 are unplaced
scaffolds** — every one flagged `ORIENTATION_SUSPECT`, because a small unplaced contig
projecting inverted has no forward component to compare against. 2,615 queries carry that flag,
which makes `any_artifact_flag` nearly useless at high carrier counts. The table should
separate chromosome-scale from unplaced carriers rather than pooling them.

**The sort is inverted for the interesting case.** 160 single-carrier loci sort first — all
unflagged, all small — while the 10.93 Mb locus with six chromosome-scale carriers sorts last.
Unflagged-first was right when artifacts were assumed to be the noise; the actual noise is
single-carrier unplaced projections.

**§0.1 needs revising.** On the rebuilt graph the chr10 locus is carried by **six
chromosome-scale haplotypes across three individuals** (`CBau_104` ×2, `CLim_110` ×2,
`CTlk_104` ×2), both haplotypes each — not the 2-of-10 heterozygous pattern recorded earlier.
That is closer to fixed divergence between the reference and those three individuals, with
`CMat` and `CPla` as the exceptions. The earlier reading came from the pre-patch graph and from
looking at 43 queries.

---

## 4a2. Batch 6b — `PRIVATE_ENRICHMENT`: the GLMM

`PRIVATE_JOIN` emits `<taxid>.<flavor>.private_evidence.csv`, one row per segment, carrying
both response forms so a binomial and a continuous model come off the same file:

`haplotype, individual, sample, flavor, chromosome, segment, set, start, end, span_bp,
log_span, is_private, n_other_assemblies, best_identity, aligned_frac_merged, map_verdict,
n_kmers_observed, n_kmers_expected, frac_absent, mean_copy, median_copy, max_copy,
single_copy_ref, copy_ratio, kmer_verdict, repeat_like, combined`

`individual` is derived so the nested random effect works. `log_span` is precomputed.

### The model as specified

```
repeat_like ~ is_private + log(span_bp) + (1 | haplotype) + (1 | chromosome)
                                        + (1 | haplotype:chromosome)
```

- **Binomial as the headline** (`repeat_like` 0/1) for an interpretable odds ratio; **lognormal
  on `log(median_copy)`** as the effect-size check, since the response spans 14 to 360,000.
- `log(span_bp)` is **not optional** — longer segments span more repeat classes, and private
  and control size distributions will not match exactly even after size-matched sampling.
  Without it the privateness coefficient absorbs length.
- `(1 | haplotype:chromosome)` is the term that answers the chr8/chr9 question: large variance
  means the effect is haplotype-specific rather than a property of the chromosome. Given
  `repeat_like` is flat at 0.92–0.94 across all fifteen, expect this near zero.
- With 10 haplotypes and 5 individuals there is almost no d.f. to estimate an individual-level
  variance separately, so `(1 | individual/haplotype)` will likely fail to converge or return
  zero. Fit `(1 | haplotype)` and note the within-individual correlation, or check empirically
  whether the two haplotypes of a sample behave alike.
- ~380,000 rows on the full arm; the binomial will want `nAGQ=0`.

### What could invalidate it

Segments within a chromosome are **not independent** — adjacent private segments often flank
the same repeat array, so their multiplicities are correlated beyond what the chromosome random
effect absorbs. Effective sample size is much smaller than nominal and standard errors will be
too small. Options: cluster-robust variance by locus, or thin to one segment per window. Worth
reporting rather than quietly emitting a *p*-value that is too good.

Also: this tests **association, not causation**. Repeat-enrichment of private sequence is
consistent with graph collapse producing spurious private sequence, but equally with genuine
repeat expansion being genuinely haplotype-specific. `PRIVATE_MAP` separates those — collapse
means the sequence exists in other assemblies, a real expansion means it does not. The two
together are the argument; neither alone.

### OPEN QUESTION for the model: graph-node sharing does not imply block alignability

**117,492 of 315,303 control windows have zero other-assembly hits** — 79% of control
failures — despite being ≥95% cross-individual by graph coverage. The graph says two
individuals walk those nodes; minimap2 says the sequence is not in their assemblies at ≥65%
identity over ≥50% of its length. The remaining 29,145 failures have 6–9 hits but sit at
0.1–0.3 merged coverage. Windows that succeed do so decisively: **162,908 at ≥0.9**.

Both statements can be true simultaneously: a window can be 95% covered by cross-individual
nodes while every individual node is only a few hundred bp, so the *window* has no contiguous
homologue even though its parts are shared. That is a real difference between graph-node
sharing and assembly-level alignability, not a defect in either measure — and it **caps the
control at ~54% `NOT_PRIVATE`**, which is the honest ceiling this measure can reach.

Three rounds of threshold tuning moved control from 0.19 → 0.5312 → 0.5360; the last change
bought 0.5%, so this is not a tuning problem. It is a question the CSV can answer directly,
and it bears on the model:

- Does the zero-hit control population differ from the ≥0.9 population in `span_bp`? If the
  failures are the short windows, node-level sharing at small scale is the explanation.
- Is `n_other_assemblies` bimodal *within* the control, and does that structure track
  chromosome or haplotype?
- **Should `aligned_frac_merged` enter the model as a covariate rather than only as a
  threshold?** The bimodality suggests two distinct populations, and collapsing them to a
  binary `map_verdict` may be discarding the informative axis.
- If graph sharing and alignability are measuring genuinely different things, the
  `map_verdict × kmer_verdict` cross-tab has a third dimension and the four-cell reading in
  §2 is an approximation.

Current contrast, for reference when fitting:

| | mapping `not_private` | k-mer `repeat_like` |
|---|---|---|
| private | 0.2372 | 0.9306 |
| control | 0.5360 | 0.0669 |

The k-mer arm does not depend on cross-assembly alignment and is unaffected by any of this.

---

## 4b. Batch 6 — pangenome construction: `--lastTrain` and `--gref`

The most expensive batch (a rebuild), so the comparison is specified in advance rather than
explored afterwards.

### Why it is two runs, not one

`--lastTrain` is **already present in v3.1.4** (added v2.9.3; bug-fixed in v3.0.1 for
`--mgSplit` interaction and v3.1.1 for a crash on unplaced contigs). `--gref` needs
**v3.2.1** — v3.2.0 shipped it as a prototype whose option name was incompatible with vg, and
3.2.1 fixed that and added `--grefL`.

Bundling them would confound the scoring change with a version bump that also moves vg
(1.71 → 1.74), odgi, gfaffix, abPOA, bcftools and taffy. So:

| run | image | change | isolates |
|---|---|---|---|
| **B** (in flight) | v3.1.4 | patched harmonization, `--vcf full clip`, topological classifier | the current baseline |
| **C1** | **v3.1.4, unchanged** | `+ --lastTrain` | the alignment scoring, single variable |
| **C2** | v3.2.1 | `+ --gref` (± lastTrain if C1 wins) | reference-bias representation |

Run A (pre-batch-1: original harmonization, clip-only VCF, awk classifier) is preserved in
`$T/pre_rebuild_baseline` but confounds three changes at once and is a provenance record
rather than a comparison arm.

### C1 — RUN AND FALSIFIED. `--lastTrain` is not the explanation for SUBST.

**Result: the hypothesis is wrong, cleanly, and the pre-registered predictions are what make
that a usable answer rather than an ambiguous one.**

| prediction | B | C1 | verdict |
|---|---|---:|---|
| SUBST falls **substantially** | 1,133,130 | 1,119,141 | **−1.2% — no** |
| INS rises | 469,002 | 496,889 | +5.9%, right direction, trivial magnitude |
| DEL rises | 249,032 | 273,134 | +9.7%, right direction, trivial magnitude |
| megabase bubbles decompose (`alt_alleles_per_record` falls from 1.3114) | 1.3114 | 1.3120 | **rose — no** |
| `INV_PATH_EXPLICIT` rises from 27 | 27 | 28 | +1, noise |
| `INV_ALN_RESCUED` falls from 294 | 294 | 298 | **rose — wrong direction** |

**The decisive evidence is the SUBST homology subdivision, not the counts.** If HOXD70 scoring
were producing substitutions where indels belong, `SUBST_HOMOLOGOUS` — the near-perfectly
aligned alleles that should have become INS/DEL — would have collapsed specifically. It did
not:

| | B | C1 | Δ | share |
|---|---:|---:|---:|---|
| `SUBST_HOMOLOGOUS` | 153,498 | 150,632 | −1.9% | 75.5% → 76.3% |
| `SUBST_PARTIAL` | 13,195 | 12,766 | −3.3% | |
| `SUBST_UNRELATED` | 36,403 | 34,065 | −6.4% | |

A flat 2–6% reduction across every category with the **proportions unchanged**. The 1–5 kb
unresolved-ratio profile shifts by the same uniform ~3%. Training the matrix made the
alignment marginally tighter everywhere and changed nothing structural.

**So SUBST at ~1.12M alleles, 76% of it homologous, is a property of the graph decomposition
or of the assemblies — not of the scoring matrix.** That answers open question 5 in the
negative.

**The under-alignment risk did not materialise.** Graph total 1,927,069,595 vs 1,926,884,214
(+0.0096%); clip private 812,774,152 vs 812,642,040 (+0.016%). Essentially the same graph.

**Unexplained and worth noting:** SNP fell 2.45% (17,741,390 → 17,306,442) and INDEL 2.35%
(11,720,570 → 11,445,277) — 710,241 fewer small variants against only ~38,000 more
INS/DEL/SUBST, total alleles down 2.15%. Something did change in the alignments; it was not
the thing predicted.

**DECISION: `--lastTrain` is KEPT** (`pangenome_cactus_lasttrain = true`). INS/DEL typing
improves slightly, nothing degrades, reverting would cost another 17 h, and C2 then tests
`--gref` against trained rather than borrowed scoring.

**What C1 changes for C2.** `--gref` targets the reference-bias asymmetry, and C1 establishes
that the asymmetry is NOT a scoring artifact — which strengthens the case for C2. The
asymmetry is untouched: `DEL` still has ~273k alleles with **exactly 0 novel node bp** against
a 94 Mb reference footprint, while `INS` has ~497k alleles with 246 Mb novel and a 5 Mb
footprint.

---

### C1 — the original hypothesis, and what would falsify it (retained for provenance)

The default scoring is derived from the HOXD70 matrix, which the cactus docs describe as
appropriate for **very diverged** genomes, noting that for pangenomes it can produce "long
runs of transitions that really should be gaps". We are aligning ten haplotypes of one
species with it. Measured consequences that fit that failure shape:

- `SUBST` is 1,133,177 alleles, 60% of the SV catalog
- 273 of 290 SUBST alleles ≥500 kb sit at **0.9–1.0 forward** homology — near-perfect
  alignment, still called allele replacement
- **13 loci carry 7 alt alleles each**: one megabase bubble per locus, no internal
  decomposition, essentially every haplotype contributing its own allele
- per-allele SV bp sums to 2.95 Gb against a ~1 Gb reference

**Predictions if the hypothesis holds.** SUBST count falls substantially; INS and DEL rise as
those events are correctly typed; the 7-allele loci decompose into smaller bubbles;
`alt_alleles / records` falls from 1.31; large SUBST at high forward homology largely
disappears.

**If SUBST does not move, the hypothesis is wrong** and SUBST is real divergence or an
assembly property — which is itself the answer to open question 5.

**A second, independent prediction worth checking.** Better alignment should convert some
alignment-rescued inversions into path-explicit ones, because path-explicit representation
requires the aligner to have recognised inverted homology. So `INV_PATH_EXPLICIT` should
**rise** from 29 while `INV_ALN_RESCUED` falls from 411. That is a direct test of the Romain
et al. framework on our own data, and it is a prediction in the opposite direction from the
SUBST one, so a confounding artifact is unlikely to satisfy both.

**Risk.** last-train infers its model from the reference against the *most diverged input*.
With within-species haplotypes that is barely diverged, so the fitted model could be tight
enough to under-align. Watch graph total bp and private-sequence totals for a collapse.

### C2 — `--gref`, and a hazard that must be handled first

`vg paths -u` computes a reference path cover: it finds the parts of the graph the reference
does not walk and promotes fragments of haplotype paths covering them into synthetic
reference paths, forming a new sample `gref_<reference>` with fragments suffixed `_<N>_alt`.
Deconstructing against that sample puts non-reference material into the VCF.

**Predictions.** The 5,913 fixed-alternative (`AC=9`) `SV_DEL` alleles — sequence present in
the reference and absent from every other haplotype — should shrink or vanish, since that
material becomes representable. `INS`'s merged reference footprint (5.2 Mb against 246.6 Mb of
novel sequence) should rise sharply. The AF-spectrum asymmetry, where `SV_DEL` is the only
class with a real frequency spectrum, should reduce.

**THE HAZARD: `gref_<reference>` is a new SAMPLE in the graph.** Every haplotype-based
analysis we have groups paths by PanSN `sample#haplotype`:

- `gfa_hap_coverage.py` would count `gref_*` as an 11th haplotype. Because gref paths are
  *copies* of existing sequence, every node they cover gains +1 coverage — so nodes that were
  private (`cov == 1`) become `cov == 2` and **stop being private**. The 813 Mb / 1,511 Mb
  private figures would be silently destroyed.
- panacus `--groupby-haplotype` would do the same to the growth and openness curves.
- `PANGENOME_UNTANGLE` would project gref paths as queries against the reference — mostly
  self-alignment noise.

So C2 is **not** just adding a flag. It requires `gref_*` path exclusion in
`gfa_hap_coverage.py`, `PANGENOME_GROWTH`, and untangle's query list, plus a check that the
private total is unchanged by enabling it. Build and verify that exclusion *before* running
C2, or the comparison is worthless.

### Possible shortcut: `--gref` without a rebuild

v3.2.0 added bypass options to `cactus-graphmap-join` to rerun indexing and VCF generation
**without re-clipping graphs**. Our pipeline uses the all-in-one `cactus-pangenome`, but the
`.gbz`, `.full.gbz` and per-chromosome graphs are all on disk. If graphmap-join accepts an
existing graph, `--gref` becomes a cheap VCF regeneration rather than a 17 h rebuild, and C1
and C2 decouple entirely.

Cheap to test once run B lands: pull the v3.2.1 image and read
`cactus-graphmap-join --help` for the bypass options and whether `--gref` is accepted there.

### Comparison harness

The comparison is only as good as the tables it runs on, so all three arms must be classified
by the **same** `classify_variants.py` — including re-running it on run A's VCF, which was
never classified topologically. Metrics, in priority order:

1. `primary_class` counts and bp (SUBST / INS / DEL / INV_*)
2. `alt_alleles / records`, and the count of loci with ≥5 alt alleles
3. `INV_PATH_EXPLICIT` vs `INV_ALN_RESCUED` split
4. SUBST homology profile at ≥500 kb (the 273-of-290 figure)
5. graph total bp, node count, private-sequence total per flavour
6. untangle inverted bp on chr10, and whether the §0.1 candidate survives
7. AF spectrum by class, and the `AC = AN` population


---

## 4c. Batch 7 — Swave as a locus-level annotation layer (EXPLORATORY)

Swave (Nat Genet 2026) was run to completion on the **pre-rebuild** graph and compared
against our catalog. Verdict: **evaluated, promising, NOT adopted** -- and the candidate
integration is not the one originally imagined.

### What was compared
Swave's `swave.sample_level.split.vcf` (2,124,066 records) against
`373251.clip.parent.all.sv_sizes.tsv` (31,314,795 alt alleles). Both derive from the same
`raw.vcf.gz`; Swave carried `LV=0` on every record, so the source population matches our
parent tier. Sample columns are identical (5, four diploid + one haploid), and `AN` varies
1--9 with a mode at 9, so its allele counts are real.

### THE FINDING: Swave MERGES, we DECOMPOSE
This is the whole result and it invalidates any one-to-one comparison.

- A Swave call spanning >10 kb overlaps a mean of **10.5 of our alleles**.
- Swave `INS` calls routinely span 4--87 kb of REFERENCE (`END - POS`), i.e. loci where
  reference sequence is replaced by longer alternative sequence, netted to "insertion" on
  balance of bp. We call the same loci `SUBST`, because the traversal both gains and loses
  node content. Neither is wrong; they are different summaries of one event.
- Under 50% RECIPROCAL overlap, 1,394,646 of 2,124,066 Swave records (65%) matched nothing --
  which looked alarming and was **an artifact of the size mismatch**, since reciprocal overlap
  penalises exactly a merged-vs-decomposed comparison. Under asymmetric containment only
  203,320 of our ~747k SV alleles have no overlapping Swave call at all.

**So Swave supplies something our catalog lacks: a DIRECTION for the locus.** 381,529 of our
`SUBST` alleles sit inside a Swave `DEL` and 90,452 inside a Swave `INS` -- where our
allele-level topology says "gains and loses node content," Swave says "this locus is net a
deletion of N bp."

### Candidate integration
Not a replacement for Layer 2 and not a competing classification. A **locus-level direction
annotation** joined onto the catalog: "these N alleles constitute a net INS/DEL/DUP of M bp."
That is additive, it is what the per-allele view cannot produce, and it requires trusting
neither `SVLEN`, the genotypes, nor the absent `QUAL`.

The join MUST be Swave-locus to our-allele-SET, by interval containment. A positional join on
`CHROM:POS` is wrong twice over: Swave emits ~1.37 records per position (1,889,887 LV values
across 1,383,522 positions), and the intervals differ in scale by ~10x.

### Two findings that survive every keying scheme
1. **`13 INV` matched `INV_PATH_EXPLICIT` under 50% reciprocal overlap**, against our 29
   path-explicit alleles. Small sets, strict criterion, two methods with nothing in common --
   a trained RNN over realignment dotplots versus signed node traversal. The strongest
   cross-validation in the comparison.
2. **Swave emits ~15,500 duplication-class calls (`DUP` 12,712, `invDUP` 1,964,
   `hyperCPX_DUP` 816, `hyperCPX_invDUP` 69) where `AT` can express none** -- only 27 alt
   alleles in 3,268,312 revisit a node. The capability gap is real even though the specific
   loci are not established (see below).

### What must be resolved before any integration
- **`SVLEN` is uninterpretable on path-valued REF/ALT.** Swave's REF/ALT hold graph PATHS, not
  sequence, unless `convert_seq` is run -- so `length($5)` counts path notation. Observed:
  a record with `END == POS` and `SVLEN=1671` but `|ALT|`=12; another spanning 87,201 bp with
  `SVLEN=90,603`. **Any size analysis must use `END - POS`, not `SVLEN`.** This retires the
  earlier "all 317 `INV` under 100 kb" result until redone.
- **Containment inflates small-variant cells.** `59,816 SNP <-> DUP`, `43,808 SNP <-> INS`,
  `22,131 INDEL <-> DUP`: a 1 bp SNP clears "50% of ours" against any overlapping call
  trivially. These rows are geometry, not agreement -- and they undercut the specific
  duplication loci, though not the existence of the gap.
- **No confidence stratification.** `QUAL` is uniformly `.`, so Swave's own guidance to drop
  LowQual cannot be applied to this output.
- **Genotypes are ~10x more conservative than ours** (CPla 26.5% vs 7.74% missing, CBau 7.2%
  vs 0.62%) but the RANKING is identical and tracks assembly quality, and hap1/hap2 rates
  match within half a percent. Not broken, just stricter. Use our `AC`/`AN` regardless.

### Method note
This comparison produced four claims that were built and then retracted -- asymmetric phasing
failure, a size-threshold explanation, an anchoring offset, and the SUBST reclassification
volume. Each came from a summary statistic accepted before it was interrogated; the
`1,889,887 > 1,383,522` line that exposed the keying error was sitting in output already
requested. **Any future tool comparison should start by establishing record granularity on
both sides before computing a single cross-tab.**

### Batch 7 scope
1. Re-run Swave on the REBUILT graph (post-batch-2), so the comparison is against the current
   catalog rather than the superseded one.
2. Redo all size analysis on `END - POS`.
3. Build the locus-level annotation join (Swave locus -> our allele set, by containment) as a
   standalone script first, validated before any module.
4. Decide adoption on whether the direction annotation changes any conclusion, not on
   agreement counts.

---

## 4d. Batch 8 — does the full arm earn its keep?

Deferred deliberately until the pipeline works end to end, then answered from measurements
rather than from the assumptions that motivated carrying both arms in the first place.

### What is now measured

| measure | clip vs full |
|---|---|
| private sequence bp | **812,983,199 vs 1,511,870,062** — +46% on full, and the reference's rank INVERTS (1 of 10 → 10 of 10) |
| inverted bp (untangle) | **665,401,148 vs 665,572,453** — 0.026%, arm-insensitive |
| variant catalog | clip only — `PANGENOME_VARIANTS` and `CLASSIFY` do not run on full |
| growth / partition | clip only — panacus runs on the clip GFA |
| graph totals | 1,926,884,214 vs 2,631,391,946 bp |

**So the full arm currently earns its keep on exactly one measure.** And the justification for
it that appears in several module headers — "clipping cuts paths into subpaths (556 vs 394 on
chr10) and a rearrangement straddling a boundary is lost to path projection" — is now
measured FALSE. The fragmentation is real; `odgi untangle` projects each subpath independently
and recovers the inversions regardless.

### The question that decides it

Is the 46% private difference a FINDING, or an artifact of what clipping is FOR?

Clipping removes sequence that had no reference alignment. If full-arm "private" sequence is
largely that unaligned material, then:

- the clip figure is the correct one to publish, and the full-arm number is measuring
  "sequence minigraph-cactus could not place" rather than "sequence unique to a haplotype";
- the reference-rank inversion has a mundane explanation — the reference is never clipped, so
  on clip it retains material every other haplotype loses;
- and the private-sequence result should be framed as a property of the clip graph, with the
  full arm cited as the bound on what clipping discards rather than as a competing estimate.

If instead the extra 698 Mb is genuine haplotype-specific sequence that merely failed to
align to THIS reference, the full arm is the honest denominator and the clip figure understates
a real biological quantity.

### How to tell them apart

The evidence needed already exists in `private_evidence.csv`, both arms:

1. **Do the full-arm-only private segments align to other assemblies?** They are in the
   evidence table with `map_verdict`. If they are overwhelmingly `PRIVATE_CONFIRMED`, they are
   genuinely absent elsewhere and clipping is discarding real sequence. If `NOT_PRIVATE`, they
   exist in other assemblies and only failed to align to the REFERENCE, which is a different
   claim.
2. **Are they repeat-like?** `repeat_like` was 0.9306 on the full arm. If the clip-only subset
   is markedly less repeat-heavy, clipping is preferentially removing repeat.
3. **The control cap.** Control tops out at ~54% `NOT_PRIVATE` because graph-node sharing does
   not imply block alignability (§4a2). Whether that ceiling differs between arms bears
   directly on whether the arms are measuring the same thing at all.

### Publication framing, which is the real deliverable

Most published minigraph-cactus work reports the clip graph exclusively. Two consequences:

- a result quoted from the full arm needs justifying against that convention, not merely
  stating;
- and if the answer is "clip is what to report", the pipeline can stop carrying the full arm
  for most measures — which would remove the flavour-parallel duplication in `UNTANGLE`,
  `REARRANGE`, `REARRANGE_PLOTS`, `HAP_COVERAGE` and the whole private chain, roughly halving
  that part of the run.

**Keep both arms until this is settled.** The cost is a few minutes per run and the comparison
is now documented rather than assumed, which is the only reason the question is answerable.

---

## 5. Open questions

### Answered by the rebuild

| question | answer |
|---|---|
| Does `.full.raw.vcf.gz` appear? | **Yes** — 4.2 GB raw, 1.6 GB filtered. The naming inference was right, and the full arm has a variant catalog. |
| Does the chr10 candidate re-derive? | **Yes**, same locus and terminal position, but with **six chromosome-scale carriers across three individuals**, not the 2-of-10 heterozygous pattern originally recorded. See §4a. |
| Does clip `hap_private.tsv` still read 812,601,533? | **No** — 812,983,199 including `repeat_traversed_bp`, and per-haplotype values moved ±0.1% bidirectionally. Harmonization moved the graph slightly; the pattern (largest gains in the renamed assemblies) is consistent with cause. |
| Do output-block additions bust a task hash? | **Still unverified.** `CLASSIFY` re-ran when its output block changed, but its script also changed, so it is not a clean test. Worth answering when a cactus rerun is cheap. |

### Live

1. ~~**Are the 288 SUBST alleles ≥500 kb a scoring artifact?**~~ **ANSWERED: NO.** Batch 6 C1
   ran `--lastTrain` as a single-variable change and SUBST fell only 1.2%, with the homology
   subdivision proportions unchanged (75.5% → 76.3% `SUBST_HOMOLOGOUS`). A scoring cause would
   have collapsed `SUBST_HOMOLOGOUS` specifically. SUBST is a property of the graph
   decomposition or of the assemblies. **The follow-on question is which** — decomposition or
   assembly — and the 13 seven-allele loci are where to look.
2. **`--gref` WILL corrupt the private-sequence analysis unless excluded first.** This is no
   longer an open question but a blocking prerequisite for C2, and batch 3 made it bigger.

   `vg paths -u` computes a reference path cover: it finds graph regions the reference does not
   walk and promotes fragments of HAPLOTYPE paths covering them into synthetic reference paths,
   forming a new sample `gref_<reference>` with fragments suffixed `_<N>_alt`. Those paths are
   **copies of sequence already in the graph**, so every node they touch gains a walker.

   | file | what breaks | fix |
   |---|---|---|
   | `gfa_hap_coverage.py` | `cov` rises on every covered node, so private sequence (cov == 1) is under-counted and the private column is silently deflated | skip paths whose sample matches `gref_*` in the coverage pass |
   | `extract_private_fasta.py` | same coverage problem, **plus** `min_cross_frac` counts INDIVIDUALS — `gref_<ref>` reads as an eleventh individual, so windows shared only with it would qualify as cross-individual when they are shared with nothing | exclude in both the coverage pass and `indiv_key` |
   | `classify_variants.py` | the synthetic sample enters `AC`/`AN`, so the AF spectrum gains a phantom haplotype and every frequency shifts | exclude from the AC/AN denominator |
   | `pangenome_popstruct.R` | a `gref_*` row in the odgi similarity matrix becomes a phantom tip in the NJ tree and a phantom point in the PCoA | filter before the distance matrix |

   **The exclusion must be a shared, named pattern, not four copies of a regex.** Four
   independent implementations of "is this a gref path" is exactly the shape that produced the
   `comment.char` bug twice.
3. **Why does graph-node sharing not imply block alignability?** 117,492 of 315,303
   cross-individual control windows have zero other-assembly hits, capping the control at ~54%
   `NOT_PRIVATE`. Framed for the model in §4a2.
4. **Singleton skew unexplained.** 96% of `SV_BLOCKSUB` private to one of five samples is not
   what segregating variation looks like.

---

## 6. Conventions

Content-anchored patch scripts, idempotent, dry-run default, `--apply` to write, `.bak`
backups, post-condition sentinels, `count(old) == 1` and `new not in text` asserted before
emitting, anchors built **programmatically from actual source** (hand-transcribing them cost
two failed patches). Delimiter-balance checks must assert the delta is **unchanged**, not
zero: `pangenome_variants.nf` is natively +1 brace / −1 paren because of embedded awk.
`nextflow.config` **and `main.nf`** are CRLF; every patch detects per file. Delivered via
`present_files`. Jason applies all edits.

`.bak` backups go to a **`deprecated/` subfolder beside the file** (`modules/deprecated/`,
`py_scripts/deprecated/`, `workflows/deprecated/`), not alongside it.

Post-condition strings must be **fully qualified**: an unqualified `min_frac     = 0.5` matched
the unrelated `harmonize_dropoff_min_frac` and blocked a correct patch.

**Rewrite rather than patch a fourth time.** Three successive anchored patches to the same
region of `pangenome_private_join.nf` stacked into overlapping copies of one guard block, and
the line-based cleanup then cut through `versions.tsv` and `stub:`.

**Every destructive command ships with its variable definitions inline.** An unset `$W` turned
a delete loop into 80 refusals, twice. The `case "$d" in "$W"/??/*)` guard did its job both
times, which is why it stays.

### Nextflow-specific, learned the hard way

- **A plain channel read more than once starves all but one consumer, silently.** Process
  outputs are broadcast; `map`/`flatMap`/`combine` results are not. Fork with `multiMap`,
  `.first()`, or an explicit split. **The symptom is a missing task, never an error.**
- **Lockstep consumption:** a process with two input channels takes one item from each per
  task. A two-item channel against a one-item channel makes **one** task and silently discards
  the rest. Pair by key with `combine` + `multiMap` rather than passing a singleton alongside a
  fanned-out channel.
- **Instrument before theorising.** `.view()` shows items flowing; `count()` only emits on
  channel close. Both together distinguish empty / never-closing / fine in one run.
- **Read `.nextflow.log`, not the progress display.** `[-] process > X -` does not mean a
  process did not run.
- **SIGPIPE under `set -o pipefail`:** no `| head` / `| tail` over a data file inside a script
  block. Single `awk` instead.
- **Resource directives (`time`, `memory`, `queue`, `errorStrategy`) are outside the task
  hash** and can be changed freely. Input tuple **shape** changes do re-run the task.
- **`0` means "derive"** for any threshold that should be data-driven, with a non-zero value
  overriding, and the derived value written into the audit. Used for
  `pangenome_private_map_min_identity` and `--single-copy`; keep the convention consistent.

### R-specific

- **Never `comment.char = "#"`** on any table whose fields can contain `#` — which is every
  table keyed on a PanSN haplotype name. Strip `#` LINES explicitly, keep `comment.char = ""`,
  and set `quote = ""` while you are there.
- **Comment placement is not standardised across the writers.** `classify_variants.py` puts
  its notes BEFORE the header; the newer scripts put them after. A reader that only tolerates
  one of those will silently return nothing for the other — which is exactly how the matrix's
  audit rows came out blank.
- **Coerce numeric columns on read.** `read.delim` returns character for any column with a
  stray token, and the failure surfaces at plot time.
- **Guard every `facet_wrap` and every log scale** after all filtering, not before.
- **Comment every figure with what it is FOR**, not what it draws. The per-chromosome
  footprint figure exists because it makes a specific past bug visually impossible; that is
  worth more in the file than a description of the axes.

---

## 7. Risks

1. **`INV_PATH_EXPLICIT` undercounts by construction.** State it with the citation wherever
   it appears, or 29 gets read as an inversion count.
2. **Non-exclusive labels break naive summing.** `primary_class` for totals only.
3. **Per-allele bp is not a total** — 2.9× the genome. Three measures exist for this reason.
4. **Singleton skew unexplained.** 96% of `SV_BLOCKSUB` private to one of five samples is not
   what segregating variation looks like.
5. **Clip-derived conclusions inherit a 46% private-sequence understatement** and the
   reference-backbone bias. Prefer the full arm for anything about private sequence.
6. **odgi workaround.** An image bump could fix the internal builder; the `-a` workaround is
   harmless either way, but the comment must survive.
7. **Tool churn.** Swave, `INVPG-annot`, GraffiTE all postdate this design. Keep the
   classification module thin enough to swap.
8. **The control caps at ~54% `NOT_PRIVATE`,** because graph-node sharing does not imply
   block alignability (§4a2). Read the private figure against that ceiling, not against 1.0.
9. **Derived thresholds are circular if the control is contaminated.**
   `max_private_frac = 0.05` bounds it, but if `control_bp_ratio_achieved` ever falls well
   below 1.0 with high rejection counts, the derived values shift toward the private
   distribution and the contrast weakens without failing.
10. **`meryl-lookup` report-type names have changed between releases.** There is no `-dump`;
    the module captures `-wig-count -help` into the task log so a future rename is diagnosable
    from the output rather than from guesswork.
11. **One guessed resource figure remains.** `pangenome_private_map` at 64 GB is sized for
    loading a ~10 Gb `.mmi` and has not been profiled. `PRIVATE_INDEX` (28.4 GB peak),
    `PRIVATE_PLOTS` (730 MB) and `REARRANGE_PLOTS` (225 MB) are now measured.
12. **The R scripts are validated statically, not executed.** There is no R in the authoring
    environment, so bracket balance, dependency lists and use-before-assignment are checked by
    script and everything else is checked by running it. Expect one round of real errors per
    new figure script; three of the four bug classes 9–12 were found that way.
13. **`pangenome_untangle_flavors` still runs both arms** even though inverted bp is
    arm-insensitive to 0.026%. Kept deliberately so batch 8 can settle the question from
    evidence, at a cost of ~5 s and a duplicate figure set per run.
