# Pangenome arm: refocus plan

**Handoff document.** Written to close out the chimera/gref session and open the next one.
Everything below is either measured in this session or flagged explicitly as untested.

---

## 0. Where things stand

The chimera arm is finished and validated. The gref comparison is done. What remains is
deciding what the pangenome arm should *be* — which analyses belong in a general-purpose
pipeline, which graph each uses, and what the output directory should look like.

### Confirmed this session

| | |
|---|---|
| chimera detection → breaking → evidence | built, run, `chr9` recovered from 1.9% to 100.8% of cohort median |
| chimera QC checkpoint + report section | wired through `qc_phase.nf`, `compile_qc.R`, `generate_summary_report.R` |
| `--gref` on cactus v3.2.1 | run, compared, **keep** |
| pangenome labels off node-local scratch | `cactus_pangenome`, `pangenome_variants`, `untangle`, `stepindex`, `hap_coverage`, `classify` |

### Open, carried into the next session

| | |
|---|---|
| **`PANGENOME_PCA_NJ` segfaults** | `odgi similarity` on the 112 GB `.og`, exit 139, 32 GB label. Worked on odgi 0.9.2, fails on 0.9.4. Test with more memory first, then against the 3.1.4 image — if the old image works it is a regression and that one process gets pinned. |
| INS/DEL asymmetry unexplained | 1.62 in clip, 1.60 in gref. `--gref` did **not** resolve it, so it is not the reference bias the plan assumed. Either real biology or a `vg deconstruct` direction artefact. |
| `Sde-CPla_115` | 755–771 Mb against a cohort of 925–1050, 135 sub-threshold composites. Whether it belongs in the cohort. |
| `CTlk_104#1` chr11 (0.759) and chr6 (0.812) | low coverage fraction, **no** composite flag, so not chimeras. Gaps or divergence. |

---

## 1. The graph rule

Three flavours come out of one cactus run. The rule is not a lookup table — it follows from
what a question counts.

```
FULL      everything, including what clipping removes
CLIP      the primary biological graph
GREF      the primary coordinate system for variants
```

### The principle

**A coverage-based question must use CLIP or FULL. A coordinate-based question should use
GREF.**

`gref_Sde-CMat_203_hap2` is a synthetic 11th sample whose 656,022 paths are *copies* of
sequence already in the graph. So:

- every node a gref path covers gains +1 coverage
- a node that was private (`cov == 1`) becomes `cov == 2` and stops counting
- the 821,682,929 bp private figure would collapse to near zero

That is not a comparability argument. It is that "covered by exactly one haplotype" is
meaningless in a graph containing duplicate paths. Excluding `gref_*` from the gref graph
reproduces the clip graph exactly, so there is nothing to gain by switching.

Conversely a variant needs reference coordinates, and 800,779,642 bp of this graph has none
in clip space.

### Measured

| | clip | gref |
|---|---|---|
| paths | 9,039 | **665,061** |
| samples | 6 real | 6 real + `gref_Sde-CMat_203_hap2` |
| VCF records | 30,765,568 | 32,441,762 (**+5.4%**) |
| **fixed-alternative sites (AC==AN)** | **2,916,967 (9.48%)** | **2,339 (0.007%)** |
| alt alleles / record | 1.321 | 1.314 |
| synthetic segments | — | 656,006 spanning 800,779,642 bp |

**The fixed-alternative collapse is the finding.** 2.9 M sites — 9.5% of the clip catalog —
were positions where every haplotype carries the alt and the *reference* is the outlier.
Those are artefacts of the reference lacking sequence, not biology, and they contaminate any
AF spectrum, SV count or classification built on the clip VCF. A 1,247-fold reduction is not
a marginal improvement.

### The assignment

| flavour | use for | why |
|---|---|---|
| **FULL** | QC; what clipping removed; accessory-sequence sensitivity; odgi viz when needed | carries the 704 Mb clipping discards, 99.1% of it private |
| **CLIP** | **primary biological graph** — graph stats, haplotype sharing, core/accessory, presence/absence, panacus, private sequence, giraffe mapping | every one of these counts coverage; gref would inflate all of them |
| **GREF** | **primary graph-to-VCF** — deconstruction, variants inside non-reference insertions, non-reference SV alleles, anything needing coordinates for sequence the reference lacks | removes 2.9 M reference-bias sites |

**Giraffe mapping is the one genuinely open case.** Listed under CLIP as the safe default, but
gref gives reads more places to align — which could improve recall or create spurious
mappings to synthetic paths. Untested.

---

## 2. Which analyses stay

The pipeline currently runs ~20 pangenome processes. Several were built to answer questions
about *this* dataset and do not generalise. Arguments both ways for each.

### Keep — infrastructure

`PANGENOME_STATS`, `PANGENOME_QC`, `PANGENOME_MANIFEST`, `PANGENOME_REPORT`,
`MULTIQC_PANGENOME`, `PANGENOME_REF_FASTA`

Cheap, universally useful, and the manifest/report are how anyone finds anything.

### Keep — the core pangenome questions

**`PANGENOME_HAP_COVERAGE` (clip + full)**
- **For:** core/shell/private is *the* pangenome question. Everything else is commentary.
- **Against:** none.
- **Verdict:** keep, clip primary, full for the sensitivity comparison.

**`PANGENOME_GROWTH` (panacus)**
- **For:** openness/saturation is the standard published pangenome figure, and it directly
  answers "have we sequenced enough individuals".
- **Against:** with 5 individuals the growth curve is short and the fit is weak.
- **Verdict:** keep. It gets better with every added individual, and the weakness is a
  property of the cohort, not the analysis.

**`PANGENOME_VARIANTS` + `PANGENOME_CLASSIFY`**
- **For:** the variant catalog is the main quantitative output.
- **Against:** vcfwave decomposition is slow (30.8 M block records) and the classification
  scheme is bespoke.
- **Verdict:** keep, **moved to the gref VCF**. The classification scheme is worth keeping
  because it separates SV classes that `vg deconstruct` lumps together.
- **TODO** are there any non-bespoke classifications schemes we can replace the bespoke method with. Ideally published implementable tools that we can replace the bespoke method with

**`PANGENOME_INPUT_COVERAGE`**
- **For:** this is what caught the chr9 chimera. It compares each chromosome against the
  cohort median and flags anything far below with a composite explanation. Cheap, and it
  catches a class of error nothing else in the pipeline sees.
- **Against:** none.
- **Verdict:** keep, and arguably promote — it is the highest-value-per-cpu-second process
  in the arm.

### Arguable — decide in the next session

**`PANGENOME_PCA_NJ` / `PANGENOME_POPSTRUCT`**
- **For:** population structure straight from the graph, no variant calling needed. Cheap
  when it works.
- **Against:** currently **segfaults**. With 5 individuals a PCA has 4 degrees of freedom and
  an NJ tree of 10 haplotypes says little that is not already obvious. Population structure
  is better answered from the variant catalog with proper population genetics.
- **Verdict:** lean **cut**, or keep as opt-in. It is a nice demonstration rather than an
  analysis anyone would publish from a 5-individual cohort.
- **TODO** should this be replaced with population structure from the variant catalog that is produced already??
- **Extra** I like having a PCA/tree showing the similarity of each haplotype. The individual level can/probably should be based on the variant catalog though

**`PANGENOME_STEPINDEX` + `UNTANGLE` + `REARRANGE` + `REARRANGE_PLOTS`**
- **For:** detects inversions and rearrangements that the VCF represents poorly. Genuinely
  hard to get any other way.
- **Against:** four processes, 139 GB of staged graphs, the most fragile part of the arm
  (SIGPIPE, scratch), and the output has not yet driven a conclusion.
- **Verdict:** **keep but make opt-in**. The capability is real; the cost is high and most
  runs will not need it.

**`PANGENOME_PRIVATE_*` (fasta, index, map, kmer, join, plots)**
- **For:** six processes that establish whether private sequence is real or an artefact —
  read-depth support, k-mer support, mapping support. That question is unavoidable whenever
  private sequence is reported, and it *is* reported by hap_coverage.
- **Against:** built specifically to validate this dataset's 813 Mb private figure. Heavy —
  needs meryl DBs, minimap2 indexes, per-haplotype mapping.
- **Verdict:** **keep the evidence chain, make it opt-in.** A pangenome that reports 42%
  private sequence without evidence it is real is not trustworthy, but not every run needs
  the full chain. Default off, on for any run whose numbers will be published.

**`PANGENOME_INVERSION_RESCUE`**
- **For:** recovers inversions that deconstruct misses.
- **Against:** very specialised; minimap2 over ~310,000 allele pairs.
- **Verdict:** **cut or opt-in.** Overlaps with untangle. Keeping both needs a reason.

**`PANGENOME_PROGRESSIVE` + `PROGRESSIVE_PLOT`**
- **For:** shows how the graph grows as haplotypes are added.
- **Against:** largely duplicates panacus growth curves.
- **Verdict:** **cut** unless it shows something panacus does not.

**`PANGENOME_2D_VIZ`**
- **For:** odgi 2D layouts are the recognisable pangenome picture.
- **Against:** already off by default, expensive, and rarely interpretable at 1.9 Gb.
- **Verdict:** keep off by default; it has no label of its own, which should be fixed if
  it is retained.

### Summary of the proposed default

```
always      stats, qc, manifest, report, multiqc, ref_fasta,
            hap_coverage (clip + full), growth, input_coverage,
            variants + classify  [ON THE GREF VCF]

opt-in      private evidence chain, untangle/rearrange,
            inversion_rescue, 2d_viz, progressive

cut         pca_nj / popstruct   (or opt-in if the segfault is fixed cheaply)
```

---

## 3. Output directory and reports

The pangenome directory is **234 GB** with ~100 top-level files, and the previous run's
outputs sit alongside the current one because publishDir does not remove what a rerun no
longer produces.

### Problems

- **Superseded figures persist.** `private_by_haplotype.png` and four others are dated Sep 11
  and were replaced by the `size_by_sharing` / `tier_by_chromosome` set in a later batch.
  Nothing removed them.
- **No flavour in the filename** for some outputs, so clip and gref results will collide.
- **The three VCF flavours are undistinguished** to a reader: `.vcf.gz`, `.full.vcf.gz`,
  `.gref.vcf.gz` with no indication which one to use.
- **112 GB `.og` and 40 GB variants VCF are published**, which is right for reproducibility
  and wrong for anyone browsing the directory.

### Proposed structure

```
pangenome/<species>/
  graphs/        the .og/.gbz/.gfa per flavour, with a README naming the primary
  variants/      the VCFs, gref primary and clearly labelled
  coverage/      hap_coverage tables and figures  (clip + full)
  growth/        panacus
  qc/            stats, input_coverage, the QC tables
  evidence/      the private chain, when run
  synteny/       untangle/rearrange, when run
  report/        pangenome_report.md, multiqc, the manifest
```

And the report should state, once and near the top, **which graph each number came from**.
A figure mixing clip coverage with gref variant counts is wrong in a way no reader would
catch.

---

## 4. Suggested order for the next session

1. **Fix or cut `PANGENOME_PCA_NJ`** — one memory bump or one container pin, and it decides
   whether the process stays.
2. **Move `PANGENOME_VARIANTS` + `CLASSIFY` onto the gref VCF.** The single highest-value
   change: it removes 2.9 M artefactual sites from every downstream variant number.
3. **Make the opt-in set opt-in** — params and gates, no logic changes.
4. **Restructure the output directory** and update the report to name its graph per section.
5. **Re-run** and compare against `tst/pangenome_pregref_20260923/`.

Items 2 and 3 are independent and can go in either order. Item 4 is cosmetic but touches
every `publishDir`, so it is worth doing in one pass rather than incrementally.

---

## 5. Conventions that survived this session

Worth carrying forward, all of them learned the hard way here.

- **Never edit inside a `script:` block** — including comments and whitespace — on a process
  that is expensive to re-run. The task hash is the literal script text. A two-word comment
  correction started a 17-hour rebuild.
- **A process moved into a subworkflow re-runs.** The fully-qualified name
  (`CHIMERA:BREAK_CHIMERAS`) is part of the hash.
- **Anchor patches on complete statements.** Four separate bugs this session came from an
  anchor matching the first line of a multi-line statement, leaving an orphaned `.mix()` or
  `.set{}` that Groovy accepts as a no-op — so the file compiles and silently loses the
  statement.
- **Post-conditions must strip comments** before searching for forbidden tokens, or they
  match the patch's own explanatory text. This happened three times.
- **Check for bare `$`, not just escaped-`$` counts.** Comparing escape densities passes
  happily while an unescaped `$0` sits in the replacement, and Groovy then fails to compile.
- **`sort | head` under `set -o pipefail` is a race, not a size limit.** It will work for
  months and then fail once.
- **Read the code before forming a hypothesis.** Several rounds this session were spent on
  theories that one `grep` would have eliminated.
