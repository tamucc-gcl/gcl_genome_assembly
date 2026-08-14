# Generalized chromosome identification & name harmonization — working plan

**Repo:** `gcl_genome_assembly` · **Owner:** jselwyn · **Opened:** 2026-08-14

**Goal.** A chromosome-identification and name-harmonization method that works in any
system, not just *S. delicatulus*; tolerates genuine between-haplotype differences; and
flows into the pangenome without losing chromosomes or over-splitting them.

**Non-goal.** Resolving *S. delicatulus* sex chromosomes. Sde is the test case that exposed
the bugs, not the target.

**Design constraint.** Every evidence channel beyond scaffold size is **optional**. The
method must still run on a species with nothing but assemblies, and must degrade to
size-only behaviour rather than failing or silently changing its answer.

---

## Status

| # | Step | State |
|---|---|---|
| 0a | Composite-aware batch consensus; voter/passenger roles; per-assembly chromosome-set metrics | **applied** (`patch_harmonize_consensus.py`) |
| 0b | Size-invariant composite-member test; join graph / presence matrix / provisional consensus map (report-only) | **patch ready** (`patch2_harmonize_graph.py`) |
| 1 | Move Hi-C mapping to post-teloclip; rename/revcomp pairs lift | not started |
| 2 | Hi-C junction evidence module | not started |
| 3 | Telomere evidence module | not started |
| 4 | Frame decoupling (consensus frame becomes the naming frame) | blocked on 2, 3 |
| 5 | BUSCO overlap evidence | not started |
| 6 | Read-depth evidence | not started |
| 7 | Merqury phase-switch evidence | not started |
| 8 | Reporting: interval, not point estimate | blocked on 4 |

---

## Architecture

The current harmonizer conflates three questions. Separating them is what makes the method
portable.

| | Question | Scope | Currently | Target |
|---|---|---|---|---|
| **A** | Is scaffold S of assembly A chromosome-scale? | within-assembly | size only (dropoff + `min_chrom_frac`) | size nominates; Hi-C / telomere / depth / BUSCO promote or demote |
| **B** | Which scaffolds across assemblies are the same chromosome? | cross-assembly | minimap2 star topology to one reference | unchanged mechanism |
| **C** | How many chromosomes are there, and which reference pieces belong together? | frame | whatever the reference happens to have | join graph over alignment + Hi-C evidence |

### Evidence-channel contract

Evidence producers are separate modules emitting two standard TSVs. `harmonize_names.py`
takes them via optional arguments; absent → current behaviour exactly. This keeps the
harmonizer free of dependencies on cooler, tidk, BUSCO and meryl, and lets each channel be
developed, cached and validated independently.

```
scaffold_evidence.tsv
assembly  scaffold  source  metric               value  call            confidence

junction_evidence.tsv
assembly  scaf_a  scaf_b  source  metric              value  call        orientation
```

`call` is the channel's own verdict in its own vocabulary; the harmonizer maps calls to
graph actions. `confidence` lets a channel self-report low-quality regions (N-rich terminal
windows, low-mappability bins) so the harmonizer can discount rather than trust blindly.

### The asymmetry rule

**Presence is informative; absence is not.** Every channel is written to this. It is what
makes non-T2T assemblies, gene-poor chromosomes and low-coverage regions safe to work with.

- An internal telomere array at a candidate junction **refutes** the join.
- No telomere anywhere **says nothing** — never a demotion criterion.
- Shared BUSCOs between two candidates **refute** them being neighbours (they're haplotigs).
- Disjoint BUSCOs **say little** — a real chromosome can be gene-poor.

### Two numbering systems, never conflated

- `chrN` — reference chromosome, descending-size rank within the chosen reference. What the
  name maps use today.
- `consN` — consensus chromosome, a connected component of the join graph, id ordered by
  smallest member. What frame decoupling (step 4) will make the naming frame.

The provisional map from step 0b already emits both side by side so the difference is
auditable before anything is renamed.

### Pangenome guarantees

| Failure | Cause | Fix |
|---|---|---|
| chromosome lost | batch-consensus demotion | step 0a — flag, don't demote |
| chromosome lost | small chromosome absorbed below `secondary_frac` | step 0b — size-invariant member test |
| over-split | reference frame has more units than the consensus | step 4 — frame decoupling |
| per-chromosome build won't group | `chrA_i+chrB_j` composites vs separate `chrA_i`/`chrB_j` | step 4 — `consN` ids |
| `+` in PanSN contig names | composite naming | step 4 — disappears with `consN` |
| real haplotype differences hidden | silent absorption | step 0b — `presence_matrix.tsv` |

---

## Step 0b — apply now

`patch2_harmonize_graph.py`, applied after `patch_harmonize_consensus.py`.

### Size-invariant composite-member test

`classify()` decided membership purely by share of a scaffold's total aligned bp
(`--secondary-frac`, default 0.2). A small chromosome joined to a large one contributes a
small share **by construction**:

```
chr20 (8.71 Mb) inside a chr4-scale scaffold (84.6 Mb)
    share = 8.7 / 93.3 = 0.093  <  0.2   -> not a member, no flag
```

Six of eight good assemblies were therefore recorded as holding a plain `chr4_1`, and chr20
vanished from the frame. The smaller the chromosome, the more invisible it became. This is
also what made chr20 look sex-restricted.

New second criterion: a chromosome is a member if ≥ `--member-cover-frac` (0.5) of **its own
length** lies inside the scaffold. Size-invariant. Uses **merged interval coverage**, not
summed block length, because block-length sums double-count overlapping alignments — observed
at 120–130% of query length on this data. Rescued members are flagged
`size_asymmetric(chrK:share=..,cover=..)`.

Also emits `inflated_aln(bp=..,cov=..)` where summed block length exceeds
`--inflated-aln-tol` × scaffold length. **Diagnostic only.** Known latent issue this
surfaces but does not change: `--min-aligned-frac` is still evaluated against the inflatable
summed value, so a repeat-rich scaffold can clear the placement threshold on double-counted
alignment. Decide separately whether to switch that to merged coverage.

### Join graph, presence matrix, provisional consensus map

Report-only — **names are unchanged**. New outputs:

| File | Contents |
|---|---|
| `{taxid}.chromosome_graph.tsv` | one row per candidate junction: `n_fused`, `n_split`, and whether each of three rules keeps it |
| `{taxid}.chromosome_components.tsv` | connected components under each rule, so the supported count range is explicit |
| `{taxid}.presence_matrix.tsv` | chromosome × assembly: `chromosome` / `composite` / `absent`, plus per-individual counts |
| `{taxid}.consensus_chromosome_map.provisional.tsv` | the schema a per-chromosome pangenome build would consume |

Three rules, all reported:

| rule | edge when | behaviour |
|---|---|---|
| `permissive` | `n_f > 1` or `n_s == 0` | current in-flag rule; passes 2f/10s, over-joins |
| `majority` | `n_f > n_s` | conservative; leaves ties open |
| `plurality` | `n_f >= n_s` | treats a tie as a join |

Reporting all three is the point: it makes the two 4f/4s ties visible as *open* rather than
silently resolved by whichever rule happens to be wired into naming.

### Verification

Synthetic case with a small chromosome joined to a large one (share 0.167), an assembly
keeping them separate, and a repeat-inflation trap: a scaffold with **9 Mb of block length**
to an 8 Mb chromosome but only **1.5 Mb of merged coverage**.

```
before:  Ind-B q_big  chr1_1         chromosome   -
after:   Ind-B q_big  chr1_1+chr4_1  composite    fusion;concordant(chr1+chr4:2f/2s);
                                                  size_asymmetric(chr4:share=0.167,cover=0.95)
trap:    Ind-F q_rep  chr2_1+chr3_1  composite    (chr4 correctly NOT a member)

join graph rule=permissive edges=2 components=2
join graph rule=majority   edges=1 components=3
join graph rule=plurality  edges=2 components=2
```

The trap is the load-bearing test: with summed block length the cover test would read
9/8 = 1.1 and wrongly promote. Merged coverage reads 1.5/8 = 0.19 and correctly declines.

### Expected on Sde

CPla_115 stays passenger. `chr4_1+chr20_1` composites appear where a plain `chr4_1` was named
in ~6 assemblies, each flagged `size_asymmetric`. The join graph gains the chr4↔chr20 edge.
Component counts land near **15 (plurality) / 17 (majority) / 16 or fewer (permissive)** —
the two 4f/4s ties are the difference. Names change, so HARMONIZE_SPECIES and everything
downstream re-runs.

**First thing to read after the run:** `chromosome_components.tsv`. If majority gives 17 and
plurality gives 15, the two ties are the whole remaining question and step 2 is the answer.

---

## Step 1 — move the Hi-C mapping (zero net compute)

### Why it's free

`FINAL_HIC_MAPS` already re-maps Hi-C after everything and emits `.pairs.gz`, `.cool`, and a
balanced `.mcool` (`run_hic_balance = true`) at 2500k/1000k/500k/250k/100k/50k/10k. The Hi-C
the harmonizer needs already exists — it's just computed one step too late:

```
906  GAP_FILLING
931  TELOCLIP_EXTEND       -> ch_final_assembly
955  HARMONIZE_SCAFFOLDS(ch_final_assembly)
965  FINALIZE_ASSEMBLY
997  FINAL_HIC_MAPS         <- maps to the FINALIZED assembly
```

Consuming it in HARMONIZE would be a dependency cycle.

### Why the reorder is exact

`finalize_assembly.nf` on the harmonized path does three things: rename `old_name → new_name`,
reverse-complement where `orient == rev`, reorder. It does **not** filter (the
`min_scaffold_bp` threshold applies only on the `NO_HARMONIZE` size-rank path) and does not
alter sequence. So post-teloclip → finalized is a pure **rename + revcomp**, with scaffold
set and lengths invariant.

Pairs therefore lift by arithmetic alone:

- `chrom`: `old_name → new_name`
- `pos`: if `rev`, `pos → L - pos + 1`
- `strand`: flip if `rev`
- `pairtools sort` against the finalized `chrom.sizes`

This is strictly simpler than what `hic_liftover_pairs.nf` already does for AGP-based lifts —
that module handles the reverse-orientation strand flip already.

### Change

Retarget `MAP_HIC_TO_FINAL` at `ch_final_assembly` (post-teloclip); add a
`HIC_LIFT_HARMONIZED_PAIRS` step consuming the name map; feed the lifted pairs to
`CONTACT_MAP_FINAL`. Mapping count unchanged at one. Final plots identical.

### Rejected alternatives

- **Add a second, earlier Hi-C mapping** — doubles the most expensive Hi-C step.
- **Move harmonization before gap-filling** — name maps would be computed on assemblies that
  then change length, forcing finalize to re-derive them; and post-teloclip is the most
  complete assembly, which is the right thing to harmonize.
- **Lift the YaHS-stage contig-coordinate Hi-C forward** — TGSGapCloser changes internal
  offsets, teloclip shifts everything by the 5′ extension, Inspector can break scaffolds.
  A lift chain through all four needs offset records from each step and is fragile.

### Build in

Teloclip extensions add terminal sequence with no Hi-C support, so extended ends appear as
low-coverage terminal bins. Given the documented over-aggressiveness (126 extensions vs ~22
tidk-confirmed), the corner test must **mask** low-coverage terminal bins rather than read
them as evidence of absence.

---

## Step 2 — Hi-C junction evidence

The only channel that can break a 4f/4s tie, and it works **within a single assembly**, so it
needs no cohort consensus.

### Test

For each candidate junction (A, B) in an assembly where they are separate scaffolds:
terminal-corner contact density, normalized against the distribution of terminal-corner
densities across all non-candidate scaffold pairs in the same assembly.

**Test all four end-pairings** (A-start↔B-start, A-start↔B-end, A-end↔B-start,
A-end↔B-end) and take the max. Orientation isn't known a priori — and *which pairing wins
gives the orientation for free*, which step 4 needs in order to reorient composites.

Corner size scales: `min(10 Mb, 25% of the shorter piece)`, floored at ~8 bins.

### Decision rule, fixed before looking

Supported if normalized enrichment exceeds the 95th percentile of the non-candidate
background **and** the signal is corner-localized. A whole-pair elevation is a
repeat/mapping artifact, not a junction. Mask low-mappability bins — the maps show strong
dark stripes throughout that would otherwise skew the background.

### Resolution: auto-select, don't hardcode

From the Sde reference maps, judging where the inter-chromosomal background sits on each
colourbar:

| resolution | inter-chromosomal background | usable |
|---|---|---|
| 100 kb | at the floor of the scale, single digits | no — can't estimate a background |
| 250 kb | ~10 | marginal; good for localization |
| 500 kb | ~10–50, crisp diagonal blocks | **yes — operating point** |
| 1 Mb | ~30–100, best SNR | yes, but smallest chromosome only 9 bins |

Portable rule: **finest resolution in the mcool where the median inter-chromosomal bin count
is ≥ ~10 raw contacts and the shorter piece spans ≥ 12 bins.** Both computable from the
`.cool` before running. The module records its choice.

### What the Sde reference maps already show

Read at the confidence a downsampled PNG supports:

- **chr8_1 ↔ chr17_1** — unambiguous off-diagonal block, strongest inter-chromosomal signal
  in the matrix, visible at 1 Mb / 500 kb / 250 kb. Matches 7f/1s alignment support.
- **chr18_1 ↔ chr19_1** — present, lower-right region.
- **chr4_1 ↔ chr20_1** — bright streak at the far right of the chr4 row. Independent
  corroboration of the H1 exclusive-footprint result.
- **chr10↔chr15, chr14↔chr16** — not resolvable from the image; those chromosomes occupy a
  compressed strip. These are exactly the two ties, so this is the part that needs the matrix.

A missed join visible in the reference's **own** Hi-C is independent of every
alignment-based conclusion. Three lines now agree the reference frame is split.

---

## Step 3 — telomere evidence (tidk)

Three-valued per scaffold end: `telomere` / `no_telomere` / `low_confidence` (N-rich terminal
window, teloclip-extended end without independent tidk support, low complexity).

| observation | strength | action |
|---|---|---|
| internal array above threshold at a candidate junction | strong | **refuse the join** — both pieces are complete chromosome ends |
| telomeres at outer ends of both pieces, none at junction ends | weak | mild support for the join |
| no telomere anywhere | none | **no action** — never demote for absence |

Thresholds are array length (monomer count), not mere presence: terminal call requires the
array within a terminal window; internal call requires it above threshold and beyond that
window. Use tidk on the assembly, **not** teloclip's extension calls.

The refuting direction is a genuine safety check on step 2's conclusions, which is why this
runs before step 4 rather than after.

---

## Step 4 — frame decoupling

The reference becomes an **alignment coordinate system only**. The naming frame becomes the
join graph.

1. Align all assemblies to the reference (unchanged).
2. Build the join graph: alignment majority **or** a Hi-C join call adds an edge; any channel
   refuting vetoes it.
3. Connected components = consensus chromosomes → `consN`.
4. Name in the consensus frame.

A split reference then just accumulates more `_2` parts in itself while every other assembly
gets a clean single name — reference choice stops propagating into names. No second
alignment pass; components come out of the pass-1 PAFs.

Also in scope:

- move reference selection into `harmonize_names.py --select-reference` so it uses the dropoff
  estimator instead of the bash `nchrom()` approximation (they disagree by 0–2 on good
  assemblies, ~240 on failed ones)
- plurality/mode over voters instead of nearest-median, N50 tie-break, ambiguity flag when the
  modal share is thin
- **reorient composites** — 22–26% of CMat/CLim sits in composites that get no strand
  normalization today, which will inflate apparent SV in the pangenome. Step 2 supplies the
  orientation.
- tighten the in-flag `concordant` rule to `n_f > n_s` with an `equivocal` tag

**Deliberately last.** Its premise is that the reference is split; establish that with Hi-C
and telomeres before renaming ten assemblies on it.

---

## Steps 5–7 — haplotig and phase-switch channels

Different failure mode from the one currently biting, hence lower priority. Real targets
exist: CTlk_104_hap1's 45 Mb `chr5_2` and hap2's six `_2` copies, which `contained_frac = 0.9`
did not catch.

**BUSCO overlap (5).** Duplicate detection only, never completeness gating. Metric: fraction
of the smaller candidate's complete-single-copy BUSCOs also on the larger. High → haplotig.
Answers a question Hi-C cannot: Hi-C shows inter-scaffold contact for *both* "two arms of one
chromosome" and "two haplotypes of the same chromosome". Needs `full_table.tsv` with
coordinates retained. Advisory only: gene-density-vs-length residual flags degenerate or
repeat-dominated chromosomes as annotation candidates, never as demotions.

**Read depth (6).** The unifying channel — haplotig ≈0.5×, collapsed repeat ≈2×, sex-limited
≈0× in one sex, normal ≈1×. One computation, four failure modes. Check first whether
`mapping_qc.nf` / `hic_coverage.nf` already emit per-scaffold HiFi depth.

**Merqury (7).** Copy-number spectra are weak for haplotig detection *within* a phased
assembly — hap1 and hap2 are each haploid representations, so nearly everything sits at 1×
regardless. Real value is **hap-specific k-mer switch detection**: alternating hap1/hap2
k-mer blocks are a phase switch, which manufactures false SVs in the pangenome graph. Worth
having, off the critical path for chromosome counting.

---

## Step 8 — report an interval

Replace the point estimate. `§Chromosome number` carries:

- per-voter chromosome-set count distribution
- component count under all three graph rules
- per-junction evidence: concordance f/s, Hi-C enrichment, telomere status, BUSCO overlap
- the presence matrix, so genuine between-haplotype differences are stated not hidden
- an explicit uncertainty statement, and the falsification conditions tested and survived

Machine-readable inputs already exist: `chromosome_sets.tsv`, `chromosome_graph.tsv`,
`chromosome_components.tsv`, `presence_matrix.tsv`, and the `# assembly_chromosome_set` /
`# chromosome_graph` report headers.

---

## Sde interval, current

Not the goal, but the running test case and the thing that will validate each step.

| basis | estimate |
|---|---|
| per-voter chromosome-set counts (8 clean) | 15, 15, 15, 15, 17, 19, 20, 26 — mode **15** |
| join graph, majority (chr8/17, chr18/19, chr4/20) | **17** |
| join graph, plurality (adds chr10/15, chr14/16) | **15** |
| four male haplotypes independently | **15** each, same junctions, same breakpoints |

**15–17 per haplotype**, and the two 4f/4s ties are the entire remaining width. chr20 is
resolved: a distinct ~7.5–7.9 Mb segment of chr4 (exclusive footprint 0.69–0.90 in five
assemblies), not sex-limited, not a haplotig. The ZW reading of chr20 is dead — it was an
artifact of the `secondary_frac` blind spot, and the one apparently informative female
observation was a single coin flip. ZW as a *system* is untested, not refuted; step 6 is
where that gets looked at, without a nominated candidate.

## Lessons to keep

- **Silent absorption is worse than a wrong answer.** Both the batch-consensus demotion and
  the `secondary_frac` blind spot produced confident, unflagged output. Every promotion and
  demotion now carries a reason string.
- **A failed assembly can carry real signal.** CPla_115 was the only non-reference male that
  revealed chr20 — precisely because fragmenting into ~1.2 Mb pieces raised chr20's share
  above the threshold that hid it everywhere else. Passengers are unreliable for **topology**
  and can be uniquely informative for **presence**. Emit passenger coverage as advisory rather
  than filtering it out of view.
- **Test the mechanism before building on a correlation.** The chr20 sex pattern had exactly
  one informative observation and was presented as replicated.
- **Fix order matters.** Removing CPla before step 0a would have demoted chr8 (51.5 Mb) and
  chr17 (26.5 Mb) to unplaced *and* silently renumbered every chromosome after them.
