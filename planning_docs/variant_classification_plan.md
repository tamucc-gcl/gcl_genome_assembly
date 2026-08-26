# Pangenome variant classification + SV/private-sequence rework — implementation plan

All design questions resolved. Swave and `INVPG-annot` are **parked pending evaluation**, not
adopted; everything here is ours, with a join hook so a later Swave comparison is a diff.

**Origin:** Chris Bird's feedback (2026-08-26) on bp-weighted SV spectra, private-haplotype
size spectra, independent mapping of private haplotypes, and transposon drivers.

**Cohort:** 5 individuals × 2 haplotypes = 10 paths. Reference `Sde-CMat_203_hap2`.
Passengers `Sde-CPla_115_hap1/hap2` stay in. Graph: 151,697,494 segments; `reference.fa` is
1015 MB, so ~1 Gb per haplotype.

---

## 1. Decisions

| # | Decision |
|---|---|
| 1 | bp-weighted SV histogram uses the **longer allele** per alt allele. y-axis labelled "summed allele bp"; **no total printed under it** — 41.9M records carry 55.3M alts, so a 7-alt record's reference span is counted 7 times. |
| 2 | Totals are two separate numbers: **merged non-overlapping reference footprint** (denominator = reference length, since it is literally a footprint on the reference) and **summed novel-sequence bp** (denominator = total graph bp, additive across haplotypes because different insertions are different sequence). The existing private-sequence plots already normalise against graph content correctly; follow them. |
| 3 | All existing plots kept; new plots additive. |
| 4 | Private-segment spectra in count and bp on the same coarse bins as the SV spectrum. |
| 5 | Private-segment FASTA floor `pangenome_private_min_bp = 1000`. |
| 6 | `PRIVATE_MAP` emits data + plot now; statistical tests deferred, so all columns land in **one joined per-segment table** and the later test is a regression on existing columns. |
| 7 | TE question deferred to GraffiTE once the annotation pipeline delivers a RepeatModeler library. Interim: meryl k-mer copy number. **Everything emits BED** so the later intersection is possible. |
| 8 | Clip and full in parallel; report becomes a test × {clip, full} matrix — **except the fine view**, see 9. |
| 9 | **Fine view (vcfwave) on clip only.** The full graph's value is large-event and private-sequence content, not SNP/indel resolution, and vcfwave is the single most expensive step (~6 h at 16 cores on the clip catalog; worse unclipped). Parent + unresolved tiers and all private-sequence work run on both flavors. |
| 10 | `INV` is the **union** of two detectors, with provenance as a column: `INV_PATH_EXPLICIT` (topology) and `INV_ALN_RESCUED` (revcomp alignment). Both reported; the decomposition stays auditable. |
| 11 | `INV_DUP` folds into `INV` — 2 records. `DUP` from `AT` dropped as a design goal: 27 node re-visits graph-wide means it is structurally unrepresentable here. |
| 12 | Unresolved tier is **classify-then-route**, sequenced after the revcomp detector. Span and ratio are the mechanism, never the sole criterion. Routed to a published BED with bp stated, never deleted. |
| 13 | Revcomp detector: 1 kb floor, k-mer strand-bias prefilter (numpy) → **minimap2** confirmation. No `edlib`; one tool covers short and long alleles and it is already in the pipeline. 80% aligned fraction against the **shorter** allele, matching how `SUBST` is defined. |
| 14 | Revcomp detector is its **own process**, not folded into `CLASSIFY` — it is the expensive new thing and its thresholds will be iterated. |
| 15 | `HAP_COVERAGE` is **not** split by chromosome. 139 composite scaffolds appear in more than one chromosome subgraph; a private run cut at a boundary inflates segment counts and deflates lengths, corrupting the exact histogram Chris asked for. |
| 16 | `UNTANGLE` gets its own resource label, not `cactus_tools` (cpus=4 / 32 GB / 4h, sized for a stats step). Profile set from the probe. |

---

## 2. Settled by measurement

| Fact | Value |
|---|---|
| vcfbub cap **adds** records by popping oversized parents | 31,365,160 with cap vs 23,877,797 without; `--max-ref-length 0` yields 0 |
| Parent view **is** `LV == 0` | 23,877,797 — exact match to `vcfbub --max-level 0`. No span floor: OR-ing one would double-count nested content inside its parent. |
| `AT` present, parseable | 41,920,948 / 41,920,948, 0 unparseable |
| `AT` invalid post-decomposition | vcfwave + `norm -m` rewrite REF/ALT; alt index *i* ≠ traversal *i+1* |
| Path-explicit INV / DUP / INV_DUP | 30 / 25 / 2 of 3,268,312 SV alleles |
| Node re-visits | 27 alleles; 3,268,285 have max visits of 1 |
| Topology partition | SUBST 1,956,621 · INS 853,248 · DEL 458,386 |
| Reclassified by topology | 233,015 length-called `SV_INS` and 178,583 `SV_DEL` are SUBST |
| Published summary stale | understates SV by 943,760, overstates INDEL by the same; deltas reconcile exactly to two awk branches added since |
| Genotypes | 5 columns (4 diploid + `Sde-CMat_203_hap1` haploid), reference absent; `AC`/`AN` on 100% of records — **AF from AC/AN, never from counting sample columns** |
| Missingness | 3.14% (CBau) to 12.84% (CPla) |
| Records ≥ 1 Mb | 28; 22 at ratio 0.9–1.3, then a clean gap through 3.1, 13.5, 78.8, 352, 4021+ |
| SUBST ≥ 500 kb | 288 alleles — the unresolved candidate set |
| Private sequence | 813 Mb; reference 123 Mb (15.1% vs 10% even); passengers **lowest** at 58–60 Mb |
| Reference bias, second signal | 30,224 `SV_DEL` alleles carried by all 5 samples = reference-specific insertions |
| Singleton skew | 97.5% of `SV_BLOCKSUB`, 94% of `SV_INS` private to one sample; `SV_DEL` only 65.7% |
| `.og` sizes | clip 118 GB, full 32 GB — backwards, probe step 2 will say whether it is content or serialisation. Whole-graph odgi work needs bigmem+; per-chromosome files are ~1.8 GB each. |

### Why path-explicit INV is 30
Romain et al. (bioRxiv 2025.03.14.643331) define two inversion motifs: **path-explicit**, where
ancestral and inverted alleles traverse shared nodes in opposite directions, and
**alignment-rescued**, where the alleles are disjoint paths and sequence alignment is needed.
`AT` sees only the former, and for Minigraph-Cactus the path-explicit set skews large (boundary
around 500 kb). Alignment-rescued inversions land in `SUBST` — which is what the revcomp
detector is for.

---

## 3. Layers

| Layer | Instrument | Delivers |
|---|---|---|
| 1 | `odgi untangle`, per chromosome | inversions, dispersed duplications, **translocations** (structurally impossible in a per-`refContig` VCF), rearrangement without clean anchors |
| 2 | `AT` traversals, parent view only | SUBST / INS / DEL partition replacing `COMPLEX`/`BLOCKSUB`; `INV_PATH_EXPLICIT` |
| 2b | revcomp + minimap2 on SUBST ≥ 1 kb | `INV_ALN_RESCUED`; also the routing input for the unresolved tier |
| 3 | `AC`/`AN` | AF spectrum, per-individual private variants, core/shell/cloud partition of the variant catalog |
| — | `gfa_hap_coverage.py`, `PRIVATE_MAP`, `PRIVATE_KMER` | private-sequence spectra and validation |

Layer 2/2b output carries `chrom`, `pos`, `primary_class`, `topology_labels`, `inv_provenance`
in a flat table so a Swave comparison is a join.

**Exclusive primary class** for anything that sums: `INV` · `DUP` · `SUBST` · `INS` · `DEL` ·
`REORDER` · `SNP` · `INDEL`. Non-exclusive label set retained per variant for the upset plot.

---

## 4. Batch 1 — harmonization, pangenome off

1. **`harmonize_names.py`: four envelope → merged-footprint swaps** (patch script delivered).
   - **1071 `prim`** — load-bearing. Picks the piece deciding a scaffold's *orientation* by
     envelope span. A wrong pick reverse-complements a chromosome-scale scaffold; the comment
     above it records that the pre-Step-4 version left 22–26% of some assemblies unnormalised
     and inflated apparent SV in the pangenome.
   - **1009/1018** overlap detection — 371 flags. `chr4_1` (42.3 Mb, 79.2 Mb envelope) and
     `chr4_2` (45.7 Mb, 73.5 Mb envelope) both report `cov≈0.98`, so their merged footprints
     tile. The flags are artifacts.
   - **1033** part-index key — order by the start of the largest merged block.
   - **1055** `ref_span` — report the merged footprint.
2. **Before/after diff** of the name map and orientation column against existing PAFs. No
   change means the envelope was adequate and we stop thinking about it; any orientation flip
   is the SV-inflation mechanism caught in the act.
3. `publish_dir_mode = 'copy'`, `stageInMode = 'copy'` on cactus.
4. Popstruct intermediates: individual `Di`, PCoA coordinates and eigenvalues, NJ Newick.
5. Doc fixes: the `pangenome_variants.nf` docstring pairs a current-code numerator (401,642)
   with an old-code denominator (1,567,141); note the published report's SV/INDEL skew.

---

## 5. Batch 2 — the rebuild

### 5.1 `cactus_pangenome.nf`
`--vcf full clip`. Add emits `gfa_full`, `og_full`, `gbz_full`, `snarls_full`, `vcf_full`,
`raw_vcf_full`, `chrom_og_full`. Invert the `.full.og` exclusion filter in
`workflows/pangenome.nf` for the full arm.

### 5.2 `pangenome_variants.nf` — decomposition only, three tiers
Classification stripped out. Per flavor:
- **unresolved** → `unresolved_regions.bed` + record table (both flavors)
- **parent** — `LV == 0` minus unresolved (both flavors)
- **fine** — vcfbub (`--max-ref-length` stays **100000**) → chunked vcfwave → norm → sort,
  plus `blocks.vcf.gz` and `bcftools stats` (**clip only**, decision 9)

Split by chromosome via tabix. Re-profile resources: 15 chromosome tasks each handle ~1/15 of
the records, so per-task `cpus` drops from 16 — otherwise 15 × 16 cores is queue-bound.

### 5.3 `pangenome_classify.nf` + `py_scripts/classify_variants.py`
Layers 2 and 3, per `(flavor, chromosome)`. Promote the validated parts of
`audit_variant_layers.py` v2 with both fixes enforced **in the module**:
- `DUP` restricted to nodes the reference also visits
- hard provenance guard refusing topology on any VCF whose header shows `vcfwave`,
  `ID=ORIGIN` or `bcftools_normCommand`

### 5.4 `pangenome_inversion_rescue.nf` + `py_scripts/rescue_inversions.py`
SUBST ≥ 1 kb: k-mer strand-bias prefilter, minimap2 confirmation, 80% against the shorter
allele. Emits `INV_ALN_RESCUED` calls plus the unresolved routing decision.

### 5.5 `pangenome_untangle.nf` + `pangenome_rearrange.nf`
Scatter over per-chromosome `.og`, copying the `PANGENOME_2D_VIZ` pattern. New resource label.
`REARRANGE` emits inversion / translocation / duplication BED, size tables, bp totals, and an
intersection against `unresolved_regions.bed` — the 25.6 Mb chr7 record at 16,349,874 is the
first test.

### 5.6 `gfa_hap_coverage.py` + `pangenome_hap_coverage.nf`
Fourth pass: collapse maximal runs of consecutive private (`cov == 1`) nodes along each
haplotype's walk into segments. Emit `private_segments.tsv` and `.bed`. Add
`--min-private-bp` (default 1000) and flavor plumbing. Whole-graph per decision 15. Leave the
three verified passes untouched.

### 5.7 `nextflow.config`
```
pangenome_vcfbub_max_ref   = 100000   // UNCHANGED; 0 empties the catalog
pangenome_parent_view      = true     // LV==0, no span floor
pangenome_fine_view_flavors = 'clip'  // decision 9
pangenome_unresolved_span  = 500000
pangenome_unresolved_ratio = 100
pangenome_graph_flavors    = 'clip,full'
pangenome_private_min_bp   = 1000
pangenome_inv_rescue       = true
pangenome_inv_rescue_min_bp = 1000
pangenome_inv_rescue_min_frac = 0.8
pangenome_untangle         = true
pangenome_sv_bins          = '50,100,500,1000,5000,10000,50000,100000,500000,1000000'
publish_dir_mode           = 'copy'
```
Only CRLF file in the repo — patches preserve line endings per file.

### 5.8 Parallelization
| Step | Axes | Fan-out |
|---|---|---|
| `VARIANTS` tiers | flavor, then chrom | 2 → 15 |
| `VARIANTS` vcfwave | chrom × chunks (clip only) | 15 × chunks |
| `CLASSIFY` | flavor × chrom × tier | 2 × 15 × 3 |
| `INVERSION_RESCUE` | flavor × chrom | 2 × 15 |
| `UNTANGLE` | flavor × chrom | 2 × 15 |
| `HAP_COVERAGE` | flavor only | 2 |
| `PRIVATE_*` | flavor × haplotype | 2 × 10 |

All merges use explicit `size:` — `groupTuple`/`.collect()` without one are scheduling barriers
gating everything on the slowest task. Counts known: 15, 10, 2, 3.

---

## 6. Batches 3–5

| Batch | Contents |
|---|---|
| 3 | `PRIVATE_FASTA`, `PRIVATE_MAP` (minimap2 each haplotype's private FASTA against every other assembly — the discriminator for the reference's 123 Mb), `PRIVATE_KMER`. One joined per-segment table. |
| 4 | `pangenome_plots.R`: bp-weighted histogram by `primary_class` faceted by tier on both bp axes; upset label co-occurrence; private-segment count and bp histograms; AF spectra; Layer 1 rearrangement bp. **Inversions get their own panel** — an inversion adds and removes no sequence, so stacking orientation-difference bp with presence/absence bp makes the column total meaningless. |
| 5 | `pangenome_report.R`: test × {clip, full} matrix, "n/a — removed by clipping" where a cell cannot exist, "n/a — clip only" for the fine view, explicit unresolved-bp exclusion line. |
| later | GraffiTE after the annotation pipeline. `svim-asm` as the non-graph check, needing the `PAIRWISE_ALIGNMENT` fixes (no `-c`, `--secondary=no`, `min_aln_bp=150000`, `within_sample` — fix by emitting two PAFs from one alignment). |

---

## 7. Risks

1. **`odgi untangle` is load-bearing and unmeasured.** Probe gates 5.5's scope. The 118 GB clip
   `.og` means whole-graph untangle may be unaffordable regardless, which makes the
   per-chromosome scatter a necessity rather than a preference.
2. **Singleton skew unexplained.** 97.5% of `SV_BLOCKSUB` private to one of five samples is not
   what segregating variation looks like. Report carries the flag until understood.
3. **Reference private-sequence excess.** 123 Mb / 15.1%, corroborated by the 30,224 all-sample
   `SV_DEL` alleles. `PRIVATE_MAP` is the discriminator; private plots carry the caveat until it
   runs.
4. **Non-exclusive labels break naive summing.** `primary_class` for totals, label set for the
   upset plot. Every output states which it uses.
5. **`INV_PATH_EXPLICIT` undercounts by construction**, stated wherever it appears with the
   Romain et al. citation, or 30 gets read as an inversion count.
6. **Tool churn.** Swave (Nat Genet 2026), `INVPG-annot` (2025), GraffiTE all postdate this
   design. Building our own is right for control now; keep the classification module thin
   enough to swap.
