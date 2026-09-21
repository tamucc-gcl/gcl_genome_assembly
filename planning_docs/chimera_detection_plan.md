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

### RESOLVED

| | before | after the break |
|---|---|---|
| chr9 in `Sde-CTlk_104#1` | 1,281,066 bp, **1.9%** of median | **69,095,920 bp, 100.8%** |
| chr9 core (all ten haplotypes) | ~0 | **17,551,310 bp, 14.0%** |
| chr9 subgraph | 118.7 Mb | 125.4 Mb |
| `Sde-CTlk_104#1` graph content | 822 Mb, lowest in cohort | **925.2 Mb**, mid-pack |
| clip graph total | 1,926,884,214 bp | **1,939,208,633 bp** (+12.3 Mb) |
| `cells_low` | 1 | **0** — every chromosome reads `ok` |

The graph GREW by 12.3 Mb, which is the under-alignment tripwire passing: sequence that could
not be placed before is now aligned and contributing, rather than the graph shrinking because
the aligner was given less to work with.

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

## 2. The evidence hierarchy

| signal | role | why |
|---|---|---|
| **cross-haplotype concordance** | **the gate** | the only signal independent of how the scaffold was BUILT |
| **AGP joins** | **where, exactly** | records the join scaffolding made, with a 100 bp gap at it |
| **PAF per component** | **which join** | most joins are same-chromosome and drop out here |
| interstitial telomeres | corroborates | independent of Hi-C; diagnostic of end-to-end fusion |
| Hi-C cross-contact | confirms | it made the join, so it cannot justify breaking it |

**Telomere ABSENCE must never veto a break.** A mid-arm fusion leaves none.

### The AGP replaced inference as the locator, and inference was closer than I claimed

Positions were first estimated by binning the reference PAF at 1 Mb and finding where the
dominant chromosome changed. Against the AGP:

| scaffold | PAF binning | Hi-C minimum | AGP join (truth) |
|---|---|---|---|
| `chr5_1+chr9_1` | 37.0 Mb (−144 kb) | 37.4 Mb (+256 kb) | **37,144,322** |
| `chr6_3+chr12_1` | 43.0 Mb (−695 kb) | 44.4 Mb (+705 kb) | **43,694,678** |

Both inference methods were within ~0.7 Mb. Neither was exact, and both land in *sequence*
rather than in the 100 bp gap the join actually is — so they cannot be cut at. The AGP gives
the position; inference cannot.

**But the AGP alone is not enough either.** `scaffold_3` has 173 joins and one is chimeric;
`Sde-CPla_115_hap1` has 919 across the assembly. Breaking at every join would undo
scaffolding. The PAF says *which* join separates two chromosomes, the concordance vote says
whether to act, and the AGP says where.

### A wrong turn worth recording

For one round I argued the junction was the **round-2** join at 41,919,050 — 4.77 Mb from the
truth — and built an explanation about repeat-rich regions misleading the PAF to account for a
discrepancy I had manufactured. Every component from 37.14 to 41.92 Mb is chr9, so that join
is chr9→chr9 and not chimeric at all.

Two things caused it. I took the largest/most-recent join as the relevant one instead of
testing every join, and I explained a surprising 4.5 Mb gap rather than checking it. The same
failure mode as the Hi-C measurement below. Hence `chimera_joins.py` tests **every** join.

### The Hi-C measurement error

The first measurement of cross-junction contact on `chr5_1+chr9_1` returned **1.209** —
elevated — and was explained as subtelomeric repeat attracting spurious contacts. Three things
were wrong: the junction was taken as 35 Mb from 5 Mb PAF bins, the null excluded ±6 Mb around
that wrong centre, and windows near either end averaged truncated blocks.

Corrected — scanned position, full windows only — it is **0.748**, and `chr6_3+chr12_1` is
**0.740**. Both depleted.

The circularity argument survives on its own: Hi-C made the join, so depleted contact across a
junction it created is not independent evidence. But the empirical support was an artifact, and
a design decision was built on it. **Full windows only, and scan for the minimum rather than
assuming its position.**

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
| `patch_harmonize_sister_vote.py` | splits the vote into sister / other-individual carriers |
| `patch_harmonize_header_fix.py` | corrects two header claims the AGP work superseded |
| `agp_joins.py` | both AGPs chained -> every join, exact, in final coordinates |
| `chimera_joins.py` | components + PAF -> which joins separate two chromosomes |
| `break_chimeras.py` (reworked) | N cuts -> N+1 pieces, repeated chromosomes named distinctly |
| `chimera_evidence.py` | written; Hi-C path untested, not wired |

### RUN AND APPLIED

`chimera_break = 'auto'`, one cut per candidate, confirmed by the audits:

| assembly | mode | broke | cuts | pieces |
|---|---|---|---|---|
| `Sde-CTlk_104_hap1` | auto | 1 | 1 | 2 |
| `Sde-CTlk_104_hap2` | auto | 1 | 1 | 2 |
| the other eight | auto | 0 | 0 | 0 |

And the finalized assembly proves the whole chain landed: `Sde-CTlk_104_hap1` now has
**separate `chr5_1` and `chr9_1` records** where it previously had one `chr5_1+chr9_1`
composite. `chr9_1` existing at all in that assembly is what was impossible before.

Evidence at the applied cuts:

| | `hap1` | `hap2` |
|---|---|---|
| cut | **37,144,208** | **43,694,648** |
| gap snap offset | −114 bp | −30 bp |
| Hi-C ratio at the cut | **0.167** | 0.740 |
| low windows contiguous | 5 | 2 |
| interstitial telomere | **12×** background, both orientations | 7× |
| telomere source | tidk | tidk |

`hap1`'s 0.167 is contact at a sixth of the scaffold median. Both cuts land within 114 bp of a
real 100 bp scaffolding gap, so neither severs sequence.

### The Hi-C pairs translation, validated on real data

`chimera_hic_pairs.py` on `Sde-CTlk_104_hap2`: 94,610,761 pairs read, 3,365,110 written for
`scaffold_3`, 1,601,553 re-canonicalised to the upper triangle, and **0 distance mismatches
across 2,717,813 intra-contig pairs**.

That last number is the validation that mattered. Orientation composition across two AGP
rounds could not be self-checked statically -- §8 records why -- and 2.7 million pairs
preserving their distances exactly is a stronger test than any fixture. The 1.6M swaps
confirm the re-canonicalisation was load-bearing: without it roughly half the contacts would
have sat on the wrong side of the diagonal, in a matrix that loads without complaint.

### A false alarm worth recording

Midway through, the evidence module flagged a Hi-C minimum at 2.0 Mb on `hap2 scaffold_3` --
ratio 0.534, MORE depleted than the cut at 0.740 -- and the region aligned to reference
`scaffold_12` rather than `scaffold_6`. I concluded a second chimeric junction had been
missed and built a block-assignment fix for it.

Reference `scaffold_12` is **`chr12_1`**. The same chromosome as the rest of that arm, so the
join at 2,974,101 is chr12→chr12 and `chimera_joins.py` was right with one transition all
along. I had inferred "different chromosome" from a different reference SCAFFOLD name --
exactly the reference-piece-versus-consensus-chromosome confusion the translation work in §5
exists to prevent, made by the person who wrote it.

The block-assignment change is harmless and verified not to create false transitions, so it
stays: it would genuinely help a fragmented terminal component. It just was not needed here.

**The 2 Mb depletion is still a real finding.** Contact is depleted across a join WITHIN
chr12, so `scaffold_36 | scaffold_6` looks like a genuine mis-join in order or orientation.
This arm detects inter-chromosomal fusions by construction and cannot see that class.

### `chimera_junction()` — SUPERSEDED

The original 1 Mb binning estimator. Replaced by `agp_joins.py` + `chimera_joins.py`, which
take the position from the AGP rather than estimating it. Retained in `harmonize_names.py`
because its output (`junction_bp`) is a useful cross-check: measured against the AGP it was
−144 kb and +695 kb on the two candidates.

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
| components | 213 | 174 |
| AGP joins on it | 212 | 173 |
| **chimeric transitions** | **1** | **1** |
| **cut at** | **37,144,322** | **43,694,678** |
| both round-1 joins, lifted | `round1_via_round2` | `round1_via_round2_rev` |
| concordance | **1f/7s**, `sis0,oth0` | **1f/7s**, `sis0,oth0` |
| transition | chr5 → chr9 | chr12 → chr6 |
| Hi-C at the junction | 0.748, five lowest windows consecutive | 0.740, three of five |
| interstitial telomere | arrays both orientations 37.6–43.7 Mb, peak 148 | none; junction 6–16 vs background 4–8 |
| terminal telomere | neither end | fwd=353 at 0.01 Mb, nothing at 76.66 |
| reading | end-to-end fusion, three signals | mid-arm join, two signals |

Both are scaffolding errors. The decisive argument remains that **the two haplotypes of one
individual fuse different chromosome pairs** — hap1 chr5+chr9, hap2 chr2+chr3, chr6+chr12,
chr7+chr11. Six of ten haplotypes have no composites and the reference has fifteen clean
chromosomes.

### Joins that are NOT chimeric, and why that matters

| join | verdict |
|---|---|
| `hap1` round-2 at 41,919,050 | chr9 → chr9. I argued for cutting here; it is 4.77 Mb from the truth. |
| `hap2` round-2 at 2,974,101 | chr12 → chr12 |
| `hap2` round-2 at 46,454,171 | chr6 → chr6. Also argued for; also wrong. |

All three would have been cut by an AGP-only rule. This is why `chimera_joins.py` tests every
join against the PAF rather than trusting the AGP alone.

### `Sde-CPla_115_hap1`: 79 callable chimeric joins

Across ~80 scaffolds, on an assembly with 919 joins total. Every one excluded by the 20 Mb
span floor, so none is a `BREAK_CANDIDATE` — but the transitions look real, and 79 chimeric
scaffolds is an assembly-quality finding in its own right. Detection reports it; the gate
declines to act. Breaking them would take ~80 scaffolds to ~160 pieces of 0.3–3 Mb, which is
dismantling an assembly rather than repairing one.

`hap2`'s other two composites also transition — `scaffold_18` chr3→chr2 at 9,749,995 (16.2 Mb)
and `scaffold_22` chr11→chr7 at 2,152,658 (7.1 Mb) — and are excluded by the same floor.

---

## 5. The detection chain as built

```
agp_joins.py       both AGPs, chained -> every join, exact, in final coordinates
chimera_joins.py   components + PAF -> which joins separate two chromosomes
harmonize_names.py the concordance vote -> whether to act        (the gate)
break_chimeras.py  N cuts -> N+1 pieces, name map rewritten
chimera_evidence.py Hi-C + telomere confirmation per join         (to wire)
```

### `agp_joins.py`

Chaining both rounds is **required, not thorough**: round 1 makes 1,830–1,891 joins per
haplotype here against round 2's 32–54, and *both real junctions are round-1 joins*. Without
lifting them ~97% of the search space is invisible.

Three cases the naive lift gets wrong, all present in real data:

- **subrange** — `scaffold_6` enters as 4,306,001–47,785,970, so round 2 broke it and used the
  middle. A round-1 position outside the used range has to be dropped, not mapped.
- **reverse orientation** — `scaffold_6` again. The offset is measured from the far end, AND
  `left_cid`/`right_cid` must be swapped, because reverse-complementing puts the round-1 left
  neighbour at the higher final coordinate. Caught by a neighbour-consistency check: before
  the fix one contig was listed as both the left of one join and the right of the next.
- **absent** — a round-1 scaffold that never entered round 2.

Validated: 1,945 joins on `Sde-CTlk_104_hap2`, **0 neighbour mismatches** across
`scaffold_3`'s 173 joins, 0 dropped.

### `chimera_joins.py`

Components, not windows. The AGP partitions a scaffold into components at contig granularity,
and a contig cannot straddle a scaffolding join by construction — a fixed window around a join
can and does.

Three rules that look like details and are not:

- **`unplaced` is not a chromosome.** It is the reference's own unassigned sequence, so it says
  nothing about which chromosome a component belongs to. Counting it as one produced a spurious
  out-and-back transition pair on an 86 kb component, taking `hap1` from 1 transition to 3.
  Only `^chr` targets vote. Measured, `unplaced` is the only non-chr target: 0.2 Mb per
  assembly.
- **Unassigned components are filled from neighbours.** `h2tg000137l_1` aligns over 0.83 Mb of
  its 4.85 Mb span and components near a junction have less — unfilled, every repeat-rich
  contig reads as a transition.
- **A component needs 100 kb aligned AND a 2× dominance margin** before it votes, so one split
  near-evenly between two chromosomes abstains rather than choosing arbitrarily.

A transition more than 250 kb from any AGP join is reported `callable=no`: the chromosome
changes **inside a contig**, which is a contig-level mis-assembly, cannot be cut at a gap, and
belongs to Inspector.

### The polymorphic-fusion safeguard

`n_f/n_s` separates an artifact from a **fixed** fusion — real biology would be carried by most
haplotypes. It cannot separate an artifact from a **polymorphic** fusion: present in one
individual, that scores `n_f=1, n_s=8`, identical to an artifact.

The discriminator is *which* haplotype carries it, because each haplotype is scaffolded
independently:

| | verdict |
|---|---|
| `n_f_other > 0` | `NOT_A_CANDIDATE` — real, possibly fixed |
| `n_f_sister > 0, n_f_other = 0` | **`REVIEW`** — possible polymorphic fusion, never auto-broken |
| both 0 | `BREAK_CANDIDATE`, "sister clean" |

Measured: both candidates are `sis0,oth0`. Supporting evidence — the two haplotypes of
`Sde-CTlk_104` fuse **different** chromosome pairs (hap1 chr5+chr9; hap2 chr2+chr3, chr6+chr12,
chr7+chr11), which a real fusion would not do.

**Inversions cannot be affected at all.** An inversion is intra-chromosomal, so it produces no
chromosome transition and this detector never fires on one. That is structural, not a threshold.

### `break_chimeras.py`

N joins → N+1 pieces. Cuts are validated as a **set** — `min_piece_bp` against adjacent cuts,
not the scaffold ends — and a single too-short piece rejects the whole scaffold, because a
partial break leaves a known chimera half-fixed and harder to reason about than an untouched
one.

Repeated chromosomes get distinct names: `scaffold_3` is chr6 / chr12 / chr6, so the pieces
become `chr6_3`, `chr12_1`, `chr6_3_2`.

### `chimera_evidence.py` — written, not wired

Hi-C and telomere **confirmation** of a known cut, no longer localisation. Per join: the
cross-contact profile with full windows only, `n_low_contiguous` (the longest run among the
five lowest — five consecutive is a boundary, one isolated window is noise), interstitial and
terminal telomere signal from tidk, and one figure.

Coordinates are the open problem. There is no usable contact map at detection time —
`QC_PHASE` is downstream of harmonization, `ch_hic_pairs_scaffold_round2_space` predates
`GAP_FILLING` (−45 kb) and `TELOCLIP_EXTEND` (**+3.7 to +5.0 Mb**), and the `_final` mcool comes
after. The plan is per-candidate **mini-reference re-mapping**: select reads by scaffold name
from the scaffold-stage BAM (names survive gap filling and teloclip — `>scaffold_1` identical
at all three stages), re-map to a single-scaffold reference, `pairtools` → `cooler cload`.
Scoped to candidates it is minutes, not hours.

**tidk is the source of record for telomeres**, run on the mini reference. Reimplementing the
count was rejected: two sources for one number means no way to adjudicate a disagreement, and
tidk normalises the canonical repeat. Note tidk reports the window END, converted to START on
read.

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
7. **Explaining a surprising measurement instead of re-checking it.** §2. Twice: the Hi-C
   1.209 (wrong position, truncated windows) and the 4.77 Mb AGP discrepancy (wrong join).
   Both times the explanation was plausible enough to stop the investigation.
8. **Taking the largest or most recent candidate as the relevant one.** The round-2 join was
   assumed to be the junction because it was the one I had just read. Testing every join costs
   nothing and would have caught it immediately.
9. **A statistic that is nondeterministic under ties.** `np.argsort` defaults to quicksort, so
   with a flat low region *which* window counted as the minimum varied between runs on
   identical data. Fixed by `kind="stable"` plus taking the longest run rather than the run
   through whichever tie led.
10. **Reverse-complement coordinate mapping.** Both the position (measured from the far end)
    and the left/right labels (swapped) change. Caught by a neighbour-consistency check --
    one contig appearing as both the left of one join and the right of the next -- not by
    reading the code.
11. **`str.index()` on a list to find position.** It returns the FIRST occurrence, so with
    duplicate values it silently pointed at the wrong element.
12. **A COMMENT inside a `script:` block is part of the task hash.** Two comment-only edits to
    `cactus_pangenome.nf` -- one of them a correction I suggested -- busted the cache and
    started a 17-hour rebuild of a graph that was already built. The rule was already recorded
    in this project and I handed over the change without applying it.

    Documentation about a process belongs in the module's HEADER docstring or in the
    `nextflow.config` label, neither of which hashes. Never touch anything between the `"""`
    markers of an expensive process for documentation reasons.

    Recovery: diff the module's script block against the cached task's `.command.sh` and
    remove every literal difference. Interpolations (`${taxid}` and the like) and shell
    escaping differ legitimately.
13. **Asserting what a cached task contained instead of reading it.** I claimed the successful
    cactus run did not include the cleanup line and reasoned a whole trade-off from it. The
    diff showed the opposite -- it DID have it, and the work directory listing had already
    said so (no `js`, no `cactus_work`).

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
4. ~~After breaking, does chr9 gain core sequence?~~ **RESOLVED, all three predictions
   held.** chr9 core ~0 -> 17.55 Mb; `CTlk_104#1` 822 -> 925.2 Mb; and chr9 itself 1.9% ->
   100.8% of the cohort median. See §0.

   Two residual low fracs on that haplotype -- `chr11` at 0.759 and `chr6` at 0.812 -- carry
   NO composite flag, so they are not chimeras. Genuine assembly gaps or divergent regions,
   and a different problem.

   The composite names (`chr5_1+chr9_1`) persist as annotations on the now-correct rows.
   That is provenance, not a live flag: status `ok`, frac 1.008.
5. **Should the span floor apply to detection or only to the gate?** Currently only the gate,
   which is why `Sde-CPla_115_hap1`'s 79 joins are reported. That seems right — suppressing
   them at detection would hide a real finding — but it means the output is dominated by rows
   nothing will act on.
6. **Is the 250 kb `callable` window right?** A transition further than that from any join is
   called contig-internal. Neither real candidate came close to the limit (51 bp and within a
   component boundary), so the threshold is untested against a genuine contig-level chimera.
7. ~~`chimera_evidence.py`'s Hi-C path is unexecuted.~~ RESOLVED: ran, and the pairs
   translation validated on 2.7M pairs. `tidk search` and `cooler cload pairs` both worked as
   written.
8. **The `chimera_broken` QC checkpoint does not instantiate.** DEFERRED, chased far enough
   to scope. Gate and scope verified correct at main.nf line 1282, `run_all_qc` demonstrably
   true (other stages stage), `BREAK_CHIMERAS` ran 10 tasks, and `chimera_broken` appears
   ZERO times in the run log.

   Most likely `BREAK_CHIMERAS.out.assemblies` read twice -- once into `ch_pre_finalize` for
   FINALIZE_ASSEMBLY at 1091, once for the QC stage at 1282 -- with the first consumer
   exhausting it. Process outputs are normally broadcast, so this needs verifying rather than
   assuming; if confirmed, the fix is a `multiMap` fork, which is the same remedy the
   channel-read-twice pattern has needed five times in this project.

   Costs only the before/after contiguity row in the QC table. Everything it would show is
   already in the break audits, so it does not block. Best diagnosed with
   `--qc_mode all_stages --run_pangenome false`, where the channel behaviour is visible in
   minutes.
9. **The assembly QC report and plots do not mention chimera detection or breaking.**
   `generate_summary_report.R` / `compile_qc.R` / the assembly report have no section for it,
   so a reader of the published report cannot tell that two scaffolds were cut, where, or on
   what evidence. Needs:
   - a summary table: candidates found, break candidates, cuts applied, per assembly
   - the per-cut evidence row (cut position, Hi-C ratio, telomere ratio, gap offset)
   - the per-candidate figure surfaced rather than left in the chimeras directory
   - the `chimera_broken` contiguity checkpoint once item 8 is fixed, since before/after
     scaffold N50 is the most legible summary of what a break did

   Scope it with batch 9 (output consolidation): both are about what the report tells a
   reader, and both are cheap once the outputs stop changing.
10. **Node-local scratch is not big enough for two pangenome processes.** `CACTUS_PANGENOME`
   failed twice at `odgi_squeeze` with ENOSPC on a 447 GB node disk -- the second time with
   the node entirely to itself, which ruled out the contention theory. The toil job store
   reaches 337 GB and squeeze needs >110 GB on top of it. `PANGENOME_VARIANTS` then failed
   the same way: vcfwave decomposes 30.8M block records from a 39 GB VCF into 64 uncompressed
   chunks.

   Fix: `scratch = false` on both labels, putting the task dir on `/work` (NFS, measured
   ~1.1 GB/s, `rsize/wsize=32768`). Every OTHER pangenome process runs fine on scratch with
   the same 112 GB graph -- including `stepindex`, `untangle` and `rearrange`, which was
   checked rather than assumed.

   Consequence: with `scratch = false` nothing discards the task dir, so `cactus_pangenome.nf`
   removes `js` and `cactus_work` explicitly at the end of a successful run. Deliberately not
   added to `pangenome_variants.nf` yet -- it is untested, and a failed run's chunks are worth
   keeping.
11. **The broken name maps are not published.** `BREAK_CHIMERAS` emits them but `publishDir`
   has a pattern only for the audit, so the provenance of each renamed piece --
   `scaffold_1_sub_0_37144208 -> chr5_1` -- lives only in a work directory. One-line fix.
