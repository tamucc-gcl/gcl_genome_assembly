# Pangenome variant classification + SV/private-sequence rework

**Status:** batches 1 and 2 built and applied. Rebuild running (`CACTUS_PANGENOME`, ~17 h,
launched with `--vcf full clip`). Batches 3–5 outstanding.

**Origin:** Chris Bird, 2026-08-26 — bp-weighted SV spectra, private-haplotype size spectra,
independent mapping of private haplotypes, transposon drivers.

**Cohort:** 5 individuals × 2 haplotypes = 10 paths. Reference `Sde-CMat_203_hap2`.
Passengers `Sde-CPla_115_hap1/hap2` retained. Clip graph 1,926,892,905 bp / 151,697,494
nodes; full graph 2,631,391,946 bp / 153,003,433 nodes.

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
| — | `gfa_hap_coverage.py` pass 5 | private-segment spectra and BED |

Exclusive `primary_class` for all totals; non-exclusive labels retained per variant.

**No chromosome scatter for `CLASSIFY` or `INVERSION_RESCUE`** — 23.9M records classify in
well under an hour single-threaded, and rescue completed as one task. `HAP_COVERAGE` is
whole-graph because 139 composites span multiple chromosome subgraphs and a private run cut at
a boundary would corrupt the segment histogram. Net effect: **no `groupTuple` barriers anywhere
in the new work** except untangle's collection, which deliberately has no `size:` because
`PANGENOME_UNTANGLE` carries `errorStrategy 'ignore'` and a fixed size would hang forever on a
failed chromosome.

---

## 3. Built and applied

**Batch 1.** `harmonize_names.py` four envelope → merged-footprint swaps; `publish_dir_mode =
'copy'`; `stageInMode = 'copy'` on cactus.

**Batch 2.** `cactus_pangenome.nf` (`--vcf full clip` + seven full-graph emits) ·
`pangenome_variants.nf` (parent tier + tier audit; awk classifier removed) ·
`pangenome_hap_coverage.nf` + `gfa_hap_coverage.py` (pass 5) · `pangenome_plots.R` (new
schemas, four figures) · `harmonize_scaffolds.nf` + `main.nf` + `pangenome.nf` (wiring) ·
`nextflow.config` (params + five labels) · five new modules · three new `py_scripts`.

The harmonization report is keyed by **filename**, not by a tupled emit, so
`harmonize_species.nf` stays byte-identical and no task hash can move.

---

## 4. Remaining

| Batch | Contents |
|---|---|
| 3 | `PRIVATE_FASTA`, `PRIVATE_MAP`, `PRIVATE_KMER`. One joined per-segment table. **Motivating question has changed:** the reference excess is resolved as a clipping artifact (§0.3), so the live question is whether chr7/chr8/chr9's elevated private fraction (§0.4) is real divergence or collapsed repeat. |
| 4 | Private-segment spectrum plots on the SV bins; per-chromosome private table (§0.4); Layer 1 rearrangement figures. |
| 5 | Report matrix, test × {clip, full}, with explicit "n/a — clip only" for the fine view. |
| later | Whole-graph untangle for translocations. `ref_span` emitting its sort key. GraffiTE post-annotation. `svim-asm` as the non-graph check (`PAIRWISE_ALIGNMENT` already runs all 45 pairs; still needs `-c`, secondaries, and a lower length filter). |

**Parked pending evaluation:** Swave (Nat Genet 2026 — supersedes Layer 2 in principle;
py3.7/EOL dependency friction and a human-trained model), `INVPG-annot` (2025 preprint).

---

## 5. Open questions

1. **Does `.full.raw.vcf.gz` appear?** The single naming inference in batch 2. All full-graph
   emits are `optional: true`, so a wrong guess costs a missing channel, not a failed task.
2. **Does the chr10 candidate re-derive?** Carrier scaffold was renamed.
3. **Does clip `hap_private.tsv` still read 812,601,533 / 122,550,615?** A change means
   harmonization moved the graph more than the null result implied.
4. **Do output-block additions bust a task hash?** Asserted twice in this project, never
   verified. Sidestepped rather than answered. Worth testing when a cactus rerun is cheap.
5. **Are the 288 SUBST alleles ≥500 kb real large-scale divergence?** Not inversions (rescue
   found zero), mostly homologous (§1). What are they?

---

## 6. Conventions

Content-anchored patch scripts, idempotent, dry-run default, `--apply` to write, `.bak`
backups, post-condition sentinels, `count(old) == 1` and `new not in text` asserted before
emitting, anchors built **programmatically from actual source** (hand-transcribing them cost
two failed patches). Delimiter-balance checks must assert the delta is **unchanged**, not
zero: `pangenome_variants.nf` is natively +1 brace / −1 paren because of embedded awk.
`nextflow.config` **and `main.nf`** are CRLF; every patch detects per file. Delivered via
`present_files`. Jason applies all edits.

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
