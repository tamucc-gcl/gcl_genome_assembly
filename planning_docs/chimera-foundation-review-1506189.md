# Chimera foundation review: run 1506189

Reviewed 2026-10-04 from `chimera-review-20261004-061953-1506189.tar.gz`.
Comparison: retained review 1505984, including its previously reviewed HiFi
evidence. This is an output review; no local pipeline tests or biological
classification were performed, and no cuts are authorized by this document.

## Execution and comparison

- SUCCESS; 5h 28m 28s; 63 executed and 237 cached tasks.
- All five hifiasm tasks and all twenty YaHS scaffolding tasks were cached.
- Eight reference-scoring tasks completed, with reported peak RSS 32.2–37.8 GB.
  This is comparable to the earlier successful run, not evidence of a large
  per-task memory increase. Isolation succeeded; the previous kills' exact
  cause remains unproven without historical node/kernel evidence.
- All ten recorded pre-finishing FASTA SHA256 values match review 1505984.
  FASTAs were not uploaded or independently rehashed locally.
- CMat hap1 remains the reference. The composite-candidate table is byte-identical
  to 1505984. Combined name maps have identical sorted lines, with row order changed.
- All pre-existing transition-table fields match for the seven CTlk intervals;
  additional fields now separate location, structural status and action.
- All pre-existing metrics in the seven evidence TSVs match the earlier review;
  local Hi-C search metrics have been added. All seven figures were generated.
- No new HiFi mapping was requested in this run. Earlier evidence can be reused
  against the matching assessment frame; it does not itself adjudicate these events.

## Foundation checks

All twenty called/review tables have the same 42-column header, including empty
tables; no row-width discrepancies were found. All ten coordinate audits report
the AGP/FASTA check passed. The reference audit explicitly reports alignment
transition assessment unavailable. Its empty table must not be interpreted as
having passed an independent structural assessment.

There are four chromosome-scope composites, all in CTlk, producing seven
diagnostic intervals. The 124 CPla composite rows remain out of the cutting scope
because chromosome inference is unresolved; they are not cleared as correct.
This dataset does not exercise a newly admitted recurrent chromosome composite,
so unchanged candidate counts do not test that branch's sensitivity.

Older-join recovery reports 14,007 of 14,094 joins recovered with exact flanks and
gap, with 87 unresolved. These genome-wide recoveries are not 14,007 justified
breakpoints. All seven target intervals remain REVIEW, with location ineligible
and `no_gap_between_alignment_anchors`.

## Target intervals

Positions below are pre-finishing diagnostic midpoints, not cut coordinates.
The Hi-C ratio is the local minimum divided by the local flank median. These
descriptive ratios are not calibrated misassembly probabilities or cut thresholds.

| Assembly | Scaffold | Transition | Midpoint (Mb) | Local Hi-C minimum (Mb) | Minimum/flanks | Minimum inside transition interval? |
|---|---|---|---:|---:|---:|---|
| CTlk hap1 | scaffold_1 | chr9 to chr7 | 62.780 | 62.900 | 0.159 | yes |
| CTlk hap1 | scaffold_1 | chr7 to chr12 | 95.963 | 95.900 | 0.439 | yes |
| CTlk hap1 | scaffold_1 | chr12 to chr4 | 116.029 | 115.200 | 0.639 | no |
| CTlk hap1 | scaffold_1 | chr4 to chr12 | 116.196 | 115.200 | 0.639 | no |
| CTlk hap1 | scaffold_5 | chr4 to chr14 | 35.513 | 35.000 | 0.242 | no |
| CTlk hap2 | scaffold_5 | chr12 to chr7 | 43.144 | 43.200 | 0.486 | yes |
| CTlk hap2 | scaffold_7 | chr14 to chr4 | 39.619 | 39.700 | 0.243 | no |

The chr9/chr7 transition has the strongest relative local depletion and is a useful
first investigation target, not a basis for a sample-specific rule. The short chr4
excursion near 116 Mb should be assessed separately from chromosome-sized arms.
All telomere profiles use tidk; motif window counts do not establish fused telomeres.

Peer-context reports find a qualifying same-individual spanning alignment for
hap1's chr7/chr12 interval, but not the reciprocal hap2 interval under the current
window/record criterion. No interval has a qualifying spanning record from another
individual. These results are not independent fusion confirmation, and lack of a
qualifying record is not demonstrated discontinuity.

## Retained artifacts and origin information

The submitted inventory reports:

- Ten standardized haplotype contig GFAs and five hifiasm logs.
- Twenty final YaHS AGPs, covering both rounds for all ten assemblies.
- Twenty YaHS BIN files and forty scaffolding log artifacts.
- No published unitig GFA, hifiasm overlap cache, graph BED annotation or intermediate
  YaHS AGP. This establishes absence from the inventoried output, not from every
  possible work directory or backup. Cached scratch-only files are not restored
  by adding optional output declarations.

Alternative-partner Hi-C output is present in the cluster inventory:
27,019,695 pairs for CTlk hap1 and 22,615,960 for hap2. The two gzip files are about
508 MB and 426 MB. Their contents were not uploaded. Projection audits report zero
ambiguous/unmapped ends and zero distance mismatches; this is coordinate validation,
not validation of biological placement or independence from scaffolding evidence.

Each of the FOUR suspect scaffolds is a single full-length forward W component
in its round-2 AGP. Therefore round 2 did not create these adjacencies by joining
separate components. That narrows the search to its input and earlier stages;
it does not distinguish earlier YaHS joins, corrections, or raw contig origin.
Do not carry coordinates backwards across correction steps without validation.

## Next bounded implementation

1. Build a generic per-interval origin ledger using retained round-2/round-1 AGPs,
   source FASTAs and hifiasm contig GFAs. Verify source correspondence across
   corrections using sequence anchors. Record gaps, component ends and correction
   boundaries throughout the uncertainty interval, not just at its midpoint.
2. Inventory the specific cached source tasks for missing raw graph/intermediate
   files before contemplating regeneration. Begin with retained final AGPs and
   contig graphs; missing unitig graphs need not block all origin tracing. Only
   regenerate an upstream stage if a particular unresolved question requires it.
3. For each resulting junction hypothesis, measure distinctive oriented anchors,
   competing placements and current versus alternative adjacency support. Reuse
   existing HiFi mappings where the coordinate frame matches, and quantify the
   opportunity for reads to bridge the anchors. Evaluate library-specific Hi-C
   support and alternative partners; aggregate pair counts are not sufficient.
4. Produce one evidence/decision row per junction, with explicit unresolved states,
   breakpoint uncertainty, phase status and action eligibility. Calibrate against
   controls before permitting automatic edits. Include an independent assessment
   route for the chosen reference.

No additional full assembly/harmonization rerun or large-file upload is needed
merely to repeat this foundation check. Preserve the current output and shared
work/cache for the targeted implementation and cluster tests.
