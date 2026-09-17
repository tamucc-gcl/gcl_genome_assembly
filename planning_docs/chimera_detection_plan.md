# Chimeric scaffold detection and breaking

**Companion to** `variant_classification_plan.md`. That document covers the pangenome
variant/private-sequence work; this one covers detecting and repairing mis-joined scaffolds
upstream of it.

**Status:** detection built and running. Evidence enrichment (`CHIMERA_EVIDENCE`) not yet
built. Two break candidates identified and independently confirmed; neither has been cut.

---

## 0. Origin

`PANGENOME_INPUT_COVERAGE` was written to join cactus's `sample-stats.tsv` to harmonization's
composite flags. The first thing it surfaced:

```
chr9_1  Sde-CTlk_104#1  1,281,066 bp     1.9% of the 68 Mb cohort median
```

`Sde-CTlk_104_hap1` has **no chromosome 9 scaffold at all**. Its chr9 is fused to chr5 in a
111.6 Mb scaffold, and `cactus-graphmap-split` assigns each contig to a SINGLE chromosome —
so the whole thing went into chr5's subgraph and **~73 Mb of genuine chromosome 9 is absent
from every graph built so far.**

It also explains why chr9 has almost no core sequence in the per-chromosome composition
figure: core requires all ten haplotypes, and one has no chr9.

**Nothing else in the pipeline catches this class of error.**

| tool | why it misses it |
|---|---|
| cactus | does not break chimeric contigs; assigns each to one chromosome and clips the rest |
| Inspector | read-to-contig alignment. A Hi-C join is an N-gap with no reads spanning it, so there is nothing to disagree with. It found two structural errors on that haplotype, neither on the fused scaffold. |
| gap filling | the absence of spanning reads that makes a join chimeric also stops TGSGapCloser filling it — the junction gaps survive as 100 bp N-runs |
| harmonization | DID detect it, flagged `chimera_suspect(ref5+ref9:1f/7s)`, and nothing acted on it |

---

## 1. Why detection belongs in `harmonize_names.py`

It already holds everything detection needs, and nothing else does:

- the PAF of every assembly against the chosen reference
- composite naming (`chr5_1+chr9_1`) identifying a scaffold spanning two chromosomes
- the **cross-haplotype concordance vote** — `1f/7s` means one other haplotype carries this
  junction and seven keep those chromosomes separate
- merged per-member footprints

A separate detector would re-derive all of it from a second alignment. `PANGENOME_INPUT_
COVERAGE` was only ever re-reading harmonization's output, and it runs inside `PANGENOME` —
after harmonization, and only when the pangenome is enabled. Breaking must happen BEFORE
harmonization renames anything, so detection belongs at the source.

---

## 2. The evidence hierarchy, and why it is ordered this way

| signal | role | why |
|---|---|---|
| **cross-haplotype concordance** | **the gate** | the only signal independent of how the scaffold was BUILT |
| interstitial telomeres | corroborates | independent of Hi-C; diagnostic of end-to-end fusion specifically |
| Hi-C cross-contact depletion | locates precisely | more precise than the PAF, but it is the signal that MADE the join |
| N-gap proximity | picks the cut point | cutting in a gap loses nothing; cutting in sequence severs real bases |

**Telomere ABSENCE must never veto a break.** A mid-arm fusion leaves none — `chr6_3+chr12_1`
is exactly that case.

### The Hi-C circularity, and a measurement error worth recording

Hi-C made the join, so it cannot independently justify breaking it. That argument stands on
its own.

But the empirical support originally given for it was **wrong**. The first measurement of
cross-junction contact on `chr5_1+chr9_1` returned **1.209** — elevated — and was explained as
subtelomeric repeat attracting spurious contacts. Three things were wrong with it:

1. the junction was taken as 35 Mb, estimated from 5 Mb PAF bins; the real minimum is 37.4 Mb,
   so contact was measured 2.4 Mb inside the chr9 arm
2. the null excluded ±6 Mb around that wrong centre
3. windows near either end averaged truncated blocks (`max(0, a-w)`), which also produced a
   spurious minimum at 2.0 Mb on the other scaffold

Corrected — right position, full windows only — it is **0.748**. The junction is *depleted*.

**Lesson: a surprising measurement was explained rather than re-checked.** The explanation was
plausible enough to stop the investigation, and a whole design decision was built on an
artifact. Requiring full windows and scanning for the minimum rather than assuming its
position are both now non-negotiable.

Hi-C is still not the gate — but it is informative, and it was informative all along.

---

## 3. What is built

### Applied

| patch | does |
|---|---|
| `patch_harmonize_chimera_detect.py` | keeps QUERY intervals per (scaffold, ref chrom) and adds `chimera_junction()` |
| `patch_harmonize_chimera_write.py` | writes `<species>.chimera_candidates.tsv` with verdicts |
| `patch_harmonize_chimera_scope.py` | carries `qivals` through `classify()` into `info`; the naming loop had referenced `ev`, which is local to the placement loop |
| `patch_harmonize_chimera_rid.py` | uses the `rid` loop variable; `a.assembly_id` does not exist — the script runs once per species off a manifest |
| `patch_harmonize_chimera_gate.py` | replaces the `n_switches` gate with an arm-coverage test |
| `break_chimeras.py` + `.nf` + wiring | the split and name-map rewrite, gated, between harmonize and finalize |
| the four `patch_main_*_gates.py` | `qc_mode 'none'` and `run_post_assembly` |

### `chimera_junction()`

Bins the scaffold at 1 Mb, credits each bin to whichever member contributes the most aligned
bp, compresses to runs, and takes the two longest runs of DIFFERENT members. Binned rather
than per-alignment because a chimeric scaffold has thousands of alignments and spurious
cross-mappings anywhere along it — the measured chr5+chr9 scaffold has ref9 alignments inside
the ref5 arm and vice versa.

### The gate

```
span >= chimera_min_span                    20 Mb
both member footprints >= min_member_bp      5 Mb
(left_bp + right_bp) / span >= min_arm_frac  0.80
n_f <= chimera_max_nf                        1
n_s >= chimera_min_ns                        3
```

**`n_switches` is reported but NOT gated on.** It counts runs of any size, so a 200 kb
spurious `chr10` alignment at bin 0 of the 111.6 Mb scaffold made a clean chr5→chr9 transition
read as "interdigitated (2 switches)" and rejected the best-evidenced chimera in the cohort.
The arm fraction separates the cases cleanly: **1.10 and 1.16** for the two real chimeras,
**0.20** for genuine interdigitation. (The ratio exceeds 1 where alignments overlap, which is
why the floor is 0.80 rather than near 1.)

### Three modes

```
chimera_break = false      detect only; candidates written, nothing cut     DEFAULT
chimera_break = '<path>'   break every row in that file -- the file IS the instruction
chimera_break = 'auto'     break what the gate marked BREAK_CANDIDATE, same run
```

The candidates file is written on **every** run regardless of mode.

---

## 4. Measured results

**139 composites**, matching the manual `fai` survey exactly:

| assembly | composites |
|---|---|
| `Sde-CPla_115_hap1` | 76 |
| `Sde-CPla_115_hap2` | 59 |
| `Sde-CTlk_104_hap2` | 3 |
| `Sde-CTlk_104_hap1` | 1 |

137 are `NOT_A_CANDIDATE` on the span floor — median 1.9 Mb, max 16.2 Mb. All `CPla_115`.
That is the right outcome: a 2 Mb chimeric fragment is a fragmented assembly, not a mis-joined
chromosome, and splitting one leaves two unplaced halves.

### The two candidates

| | `CTlk_104_hap1` chr5_1+chr9_1 | `CTlk_104_hap2` chr6_3+chr12_1 |
|---|---|---|
| scaffold | `scaffold_1`, 111,642,300 bp | `scaffold_3`, 76,662,708 bp |
| footprints | ref5 36,036,670 / ref9 73,589,204 | ref12 50,749,347 / ref6 38,365,629 |
| concordance | **1f/7s** | **1f/7s** |
| `inflated_aln` | 1.15 | 1.17 |
| arms / span | 1.10 | 1.16 |
| PAF junction | 37.0 Mb | 43.0 Mb |
| **Hi-C minimum** | **37.4 Mb, ratio 0.748** | **44.4 Mb, ratio 0.740** |
| Hi-C contiguity | five lowest windows consecutive, 37.2–37.6 Mb | three of five lowest, 44.2–44.5 Mb |
| interstitial telomere | **arrays both orientations 37.6–43.7 Mb**, peak rev=148 at 43.67 | none: junction 6–16 against background 4–8 |
| terminal telomere | neither end | **fwd=353 at 0.01 Mb**, nothing at 76.66 |
| nearest N-gap | 37.286 Mb (−114 kb) | **44.408 Mb (+8 kb)** |
| reading | **end-to-end fusion**, three signals agreeing | **mid-arm join**, two signals; intact chr12 arm from its own telomere joined to a truncated chr6 |
| **cut at** | **37,400,000** (or the gap at 37.286) | **44,408,000** — the gap |

Both are scaffolding errors, not biology. The decisive argument is that **the two haplotypes
of one individual fuse DIFFERENT chromosome pairs** — hap1 chr5+chr9, hap2 chr2+chr3 and
chr6+chr12. A real karyotypic fusion would be shared or at least heterozygous for the same
pair. Six of ten haplotypes have no composites at all and the reference has fifteen clean
chromosomes.

---

## 5. `CHIMERA_EVIDENCE` — to build

### The coordinate problem, and the solution

Evidence must be generated from **pre-split** sequence — the split has not happened yet and
the evidence is what justifies it. But there is no usable contact map at detection time:

- `QC_PHASE` (which builds the scaffold-stage maps) is at main.nf ~1160; `HARMONIZE_SCAFFOLDS`
  is at ~953
- `ch_hic_pairs_scaffold_round2_space` exists earlier and is in round-2 scaffold coordinates,
  but `GAP_FILLING` (−45 to −47 kb) and especially `TELOCLIP_EXTEND` (**+3.7 to +5.0 Mb**) run
  after. Teloclip prepends to scaffold ends, which shifts every internal coordinate.
- the `_final` mcool is produced by `FINAL_HIC_MAPS`, downstream of harmonization

**Solution: per-candidate mini-reference re-mapping.** The old BAM is used only to SELECT
reads by scaffold name, not for coordinates — and names are preserved through gap filling and
teloclip (`>scaffold_1` identical at yahs, gap_filling and teloclip). Re-mapping to the
current sequence makes coordinates correct by construction.

```
for each BREAK_CANDIDATE (assembly, scaffold):
  1. samtools faidx <post-teloclip asm> <scaffold>       -> mini reference (~100 Mb)
  2. samtools view <scaffold-stage Hi-C BAM> <scaffold>  -> reads that hit it
                                                            (name-based; positions ignored)
  3. -> fastq -> align to the mini reference
  4. pairtools parse/sort/dedup -> cooler cload -> 100 kb cool
  5. full-window cross-contact scan
```

Cheap because it is scoped to candidates: 2 scaffolds rather than 10 assemblies, ~100 Mb
reference rather than 1 Gb, ~1/10 of the reads. And it generalises — no assumption about what
ran between scaffolding and harmonization.

### Per candidate, emit

1. **Hi-C profile** — minimum position, ratio to the scaffold median, and **how many of the
   five lowest windows are contiguous with it**. That last is what made `chr5/chr9`
   convincing: five consecutive windows is a boundary, one isolated window is noise. FULL
   WINDOWS ONLY — truncated edges produced a spurious minimum at 2.0 Mb.
2. **Refined `cut_bp`** — the Hi-C minimum, with the PAF estimate kept alongside. Both
   scaffolds needed 0.4–1.4 Mb of refinement.
3. **PAF/Hi-C agreement** — `REVIEW` if they disagree by more than **2 Mb**. Not zero
   tolerance: the PAF junction comes from 1 Mb bins, so ±1 bin is expected, and the measured
   disagreements were 0.4 and 1.4 Mb.
4. **N-gap snap** — search ±500 kb of the minimum for a gap ≥100 bp; cut at its midpoint and
   record `snapped_to_gap=<pos>`, else record `no_gap_within_500kb`. A reviewer must be able
   to tell whether the cut severed sequence.
5. **Telomere profile** — max interstitial within ±2 Mb against the scaffold background, plus
   **terminal** signal at both ends. The terminal reading is what distinguished "intact chr12
   arm plus truncated chr6" from "neither end is a real chromosome end".
6. **One figure** — the contact matrix at 100 kb with the PAF estimate, Hi-C minimum, chosen
   cut and the telomere track marked. That single image is what a reviewer actually needs.

### Gating, mode by mode

`auto` cannot see telomere or Hi-C evidence if the enrichment runs after the cut — so in
`auto` the enrichment runs **before** `BREAK_CHIMERAS` within the same run, off the
mini-reference maps, and the refined verdict is what `auto` acts on.

On a run where a file was supplied, the enrichment is the record of what justified it.

**`CHIMERA_EVIDENCE` must be gated on not-yet-broken assemblies.** Run against post-split
sequence it would look for an interstitial array on a scaffold that no longer exists, find the
array at the end of one piece and the start of another — which is what a correct break
produces — and that is no longer evidence the break was justified.

---

## 6. Deferred: moving gap filling

Considered, then dropped. The proposal was to move `GAP_FILLING` and `TELOCLIP_EXTEND` after
harmonization so that gaps exist at detection time, filling cannot bridge a false join, and
Hi-C coordinates match.

**Two of the three reasons were falsified by measurement.** The junction gaps *survive* gap
filling — they are still 100 bp N-runs in the final assembly, because the absence of spanning
reads that makes a join chimeric also prevents the filler closing it. So snapping works now,
and no false join was bridged.

The third reason (coordinates) is better solved by mini-reference re-mapping, which needs
nothing moved.

**Residual risks if it is ever revisited:** the name map's `length` column would go stale
(`FINALIZE_ASSEMBLY` applies the map by NAME, so extraction still works — but the recorded
lengths would be pre-fill, and they should be recomputed in finalize regardless);
`QC_PHASE` stages a `'gap_filled'` checkpoint the report references; and harmonization's
`inflated_aln` and coverage fractions are computed against scaffold length, which changes when
gaps are filled.

---

## 7. Gates

| gate | default | turns off |
|---|---|---|
| `chimera_break` | `false` | `BREAK_CHIMERAS` |
| `qc_mode = 'none'` | `all_stages` | `QC_PHASE`, read QC, `SNAIL_PLOT_FINAL` |
| `run_post_assembly` | `true` | `FINAL_VIZ`, `FINAL_HIC_MAPS`, `COVERAGE_BOOK` |
| either of the last two off | | `REPORTING` |

Three processes are ungated on purpose: `FINALIZE_ASSEMBLY`, `PANGENOME` (gates itself), and
`BUILD_MERYL_DB` — both `QC_PHASE` and the pangenome private-kmer chain read its output
directly, so gating the producer without rewiring those reads reintroduces the failure below.

A `'chimera_broken'` QC checkpoint sits between `'teloclip'` and `'final'`, so the QC table
shows contiguity before and after the break.

**Chimera work requires ≥2 long-read assemblies per species and chromosome-scale scaffolds** —
both conditions harmonization already enforces for the concordance vote to exist at all.

---

## 8. Bug classes from this arm

1. **Groovy over-escaping in a module script block.** Writing the `.nf` as though it were a
   Python string gave `\\$(` where `\$(` was needed, and `\\\\t` for `\\t`. The lexer error
   points at one line while the fault is uniform. `cat -A` does not disambiguate backslash
   counts — count the raw bytes and compare against a module known to work.
2. **Workflow variable scoping.** A variable assigned only inside `if`/`else` is not visible
   after the block: "No such variable: ch_pre_finalize". Assign the default BEFORE the branch
   and override inside.
3. **Python scope in `harmonize_names.py`.** `ev` is local to the placement loop;
   `a.assembly_id` does not exist because the script runs once per species off a manifest.
   Both were reaching for the name expected rather than the one bound at that point. An AST
   check — are all referenced names bound in this function? — catches both, where
   `ast.parse` catches neither.
4. **Gating a process needs TWO checks, not one.** Ungated READS of gated processes ("Access
   to `.out` is undefined") and processes never gated at all are different failures, and the
   first check says nothing about the second. Four rounds were spent because only one existed,
   and because the reads check used a hardcoded process-name list rather than deriving it from
   the file.
5. **Anchors created by an earlier edit in the same patch cannot be counted.** All anchor
   counts run against the original text before any replacement, so an edit whose anchor a
   previous edit produces always reports `old=0`. Merge them.
6. **A post-condition that scans for a token will match the patch's own comments.** Strip
   comment tails before checking.
7. **Explaining a surprising measurement instead of re-checking it.** §2. Cost a design
   decision built on an artifact.

---

## 9. Open questions

1. **Should the Hi-C minimum override the PAF junction in `auto`?** Decided yes — it is being
   used to LOCATE, not to justify, and concordance already did the justifying. With the 2 Mb
   disagreement check as the safeguard.
2. **`CPla_115`'s 135 sub-threshold composites.** Excluded correctly, but they mean that
   assembly is fragmented in a way nothing currently addresses. It has the lowest read depth
   (11 vs 18–19 single-copy k-mer multiplicity) and the lowest graph content (777/760 Mb).
   Whether it belongs in the cohort is a separate question from chimera breaking.
3. **`CTlk_104_hap2`'s other two composites** — `chr2_1+chr3_2` and one more — are below the
   span floor. Worth confirming they are fragments rather than under-detected mis-joins.
4. **After breaking, does chr9 gain core sequence?** The falsifiable prediction: chr9 core
   goes from ~0, `CTlk_104_hap1` graph content rises from 822 Mb, and the full arm's ≥1 Mb
   private bin shrinks.
