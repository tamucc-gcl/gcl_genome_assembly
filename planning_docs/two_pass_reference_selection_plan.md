# Two-pass reference selection — deferred extension

**Status:** Python complete and tested. Nextflow wiring not started. Not required for
*S. delicatulus*, where the answer is already known and can be pinned.

---

## Why this exists

Reference selection currently runs in the bash block of `harmonize_species.nf`: count
scaffolds ≥ `min_scaffold_bp`, take the assembly nearest the batch median count, break ties
on N50. Four problems:

- **Wrong estimator.** `nchrom()` counts scaffolds above a size threshold; the harmonizer
  uses the drop-off caller. They disagree by 0–2 on good assemblies and ~240 on failed
  ones, so selection optimises a number the harmonizer never uses.
- **Wrong statistic.** Nearest-median assumes the batch is centred on truth. Two failed
  assemblies moved the median from 17.5 to 19.5 and chose a reference split at five
  junctions.
- **No access to roles.** Selection runs before Python starts, so passengers vote.
- **No quality dimension.** It cannot see `genome_fraction`, `cut_ratio`, or whether the
  chromosome-set call succeeded or fell back to a threshold.

A candidate can only really be judged by running the alignment and the join graph against
it. Hence two passes: score every candidate that way, then harmonize once with the winner.

## What patch 16 changed about the stakes

Frame decoupling made reference choice matter **less**, and the experiment proved it. Three
references with 15, 17 and 20 pieces produced **identical** groupings of every shared
scaffold (ARI 1.000). Only the shattered 26-piece candidate diverged (ARI 0.86–0.89).

So this is a quality-of-life improvement, not a correctness fix. It is deferred for that
reason.

---

## The metric, validated before anything was built

```
                mult   |dev|   edges   gfrac    ARI vs others
CMat_203_hap2   0.981  0.019     0     0.994    1.000   <- selected
CBau_104_hap2   1.105  0.105     5     0.989    1.000
CBau_104_hap1   1.124  0.124     2     0.991    1.000
CTlk_104_hap2   1.163  0.163    18     0.966    0.86-0.89
```

```
multiplicity = placed_scaffolds(other voters) / (n_consensus x n_other_voters)
```

1.0 is a clean 1:1 chromosome mapping. Above 1, assemblies contribute several scaffolds per
chromosome — a fragmented frame, or redundancy it failed to absorb. Below 1, chromosomes
are missing.

**The reference's own scaffolds are excluded.** They define the frame and are trivially
placed, so counting them measures which assembly was picked rather than how good the frame
is. An earlier version of this analysis included them and concluded "maximise placement",
which ranked the frames backwards.

**Selection rule:** voter → minimise `|multiplicity − 1|` → fewest edges → highest
`genome_fraction` → highest N50 → id.

**Edges is the tie-break, not the primary.** A reference that wrongly *fuses* two
chromosomes needs no edges to correct and would score best on edges alone — but queries
that keep them apart pile onto one piece and multiplicity rises. Edges catch
over-splitting; multiplicity catches over-fusion. On this data the two disagree on the
middle pair, which is why the order matters.

**Rejected:** `cut_ratio` (a property of the assembly, not of the frame it induces — kept
as a reported diagnostic); nearest-median; maximise-placement.

---

## Done

| | |
|---|---|
| `harmonize_names.py --score-only` | patch17. Runs to the join graph, emits one metrics row, writes no name maps. Reproduces the ad-hoc numbers to 4 dp. |
| `py_scripts/select_reference.py` | patch17. Applies the rule; writes `reference_id.txt` and a ranked table. Tested against real numbers, an over-fusion trap, a near-tie, and an all-passengers cohort. |
| `harmonize_names.py --list-candidates N` | patch18. Ranks candidates from **fais alone** — no PAF — so it runs before alignment while using the same drop-off caller and voter criteria as the main pass. |

All three modes verified independent: normal writes 12 files, `--score-only` 1,
`--list-candidates` 1.

## Remaining: the Nextflow

The 8h / 8 cpu `harmonize_scaffolds` label rules out doing this serially inside one task —
56 alignments at 8 threads is roughly 23 hours. It has to be parallel processes.

```
HARMONIZE_CANDIDATES   fais only, one task, fast      -> reference_candidates.tsv
HARMONIZE_SCORE        one task per candidate, parallel -> <id>.score.tsv
HARMONIZE_SELECT       collects scores                  -> reference_id.txt + table
HARMONIZE_SPECIES      existing, now takes a pinned reference
```

Sketch, all inside `harmonize_scaffolds.nf`:

```groovy
HARMONIZE_CANDIDATES(ch_species, ch_resolver)

ch_score_in = ch_species
    .join(HARMONIZE_CANDIDATES.out.list)              // taxid, ids, fas, candidates_file
    .flatMap { taxid, ids, fas, f ->
        f.readLines().findAll { it && !it.startsWith('#') && !it.startsWith('candidate') }
         .collect { line -> tuple(taxid, line.tokenize('\t')[0], ids, fas) }
    }

HARMONIZE_SCORE(ch_score_in, ch_resolver)
HARMONIZE_SELECT(HARMONIZE_SCORE.out.score.groupTuple(), ch_selector)
HARMONIZE_SPECIES(ch_species.join(HARMONIZE_SELECT.out.ref), ch_resolver)
```

`main.nf` is untouched — it stays under the Groovy 65,535-char limit, which currently has
about 1.1 KB of headroom.

### The risky part

That `flatMap` reads `readLines()` off a staged path to fan the channel out. It either works
or silently yields an empty channel, and an empty channel in Nextflow is not an error — the
downstream processes simply never run and the workflow "succeeds".

**Mitigation:** `nextflow run -stub-run` validates the whole channel topology in ~30 seconds
without running an alignment. Every process in the repo already has a `stub:` block. Write
it, stub-run it, confirm the DAG resolves and every process is instantiated the expected
number of times, and only then spend compute.

Also needed: a guard that fails loudly if the candidate list is empty, rather than letting
an empty channel look like success.

### Changes to `HARMONIZE_SPECIES`

Its input tuple gains a fourth element, the reference id file. The script reads `REF` from
it instead of running the bash selection. When `harmonize_two_pass_selection = false`, a
`NO_SELECTION` sentinel is passed and the existing bash path runs unchanged — so the old
behaviour stays reachable and the new code cannot silently become mandatory.

### New params

```groovy
harmonize_two_pass_selection   = false   // default OFF, see cost below
harmonize_reference_candidates = 8       // cap, applied by descending genome_fraction
```

`harmonize_reference_ids` continues to override everything.

---

## Cost, and why the default is off

With 8 candidates and passengers excluded as queries: **8 × 7 = 56 scoring alignments** plus
9 for the final pass. At 8 cpus each takes ~20–30 min, so ~65 tasks of ~25 min. Heavily
parallel, so under an hour of wall time given cluster capacity — but roughly **25 CPU-hours
against the ~4 spent today**.

On *S. delicatulus* that buys `CMat_203_hap2` (0 graph edges) over `CBau_104_hap2` (5). Real
but modest, given the frames are identical anyway. Six times the CPU for a better-behaved
frame is a call the user should make explicitly, not inherit from a default.

Scaling matters too: 30 assemblies at 8 candidates is 232 alignments, and lifting the cap to
all voters on 30 assemblies is 870.

---

## The escape hatch, which is the right answer for Sde

**None of this is needed for this species.** The two-pass answer is already known:

```groovy
harmonize_reference_ids = 'Sde-CMat_203_hap2'
```

The experiment established it, and pinning costs nothing. Two-pass is for the general case —
a new species where nobody has run the experiment by hand.

That also suggests a cheaper middle path worth considering before building the full thing:
keep the current bash selection, but have the harmonizer **report** the winner's
multiplicity and edge count in the report header. A user who sees `edges=5, mult=1.105` can
run the experiment by hand and pin the result. Most of the value, none of the wiring.

---

## Order of work, if resumed

1. Write the four processes and the rewiring in `harmonize_scaffolds.nf`.
2. `nextflow run -stub-run` — confirm the DAG resolves and process instantiation counts are
   right. Do not proceed on a green run alone; check the counts.
3. Run on Sde with `harmonize_two_pass_selection = true` and confirm it selects
   `CMat_203_hap2`. If it does not, the metric computed on post-teloclip assemblies
   disagrees with the one computed on finalized assemblies, which is itself worth knowing —
   all validation to date used the finalized ones.
4. Confirm `harmonize_two_pass_selection = false` reproduces today's behaviour byte for
   byte.
5. Only then consider flipping the default.
