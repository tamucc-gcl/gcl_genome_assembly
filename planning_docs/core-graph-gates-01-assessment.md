# Gate diagnostic acceptance and interpretation

Packet: `core-graph-gates-01.review.tar.gz`.

## Execution

- Exit 0; 876 tasks: 870 cached, 6 completed; no warnings/errors in the Nextflow log.
- New native compiler job: 19.7 seconds. Three policy jobs ran in parallel and took
  2m59s, 3m04s and 3m04s. Summary: 277 ms.
- Six original/corrected-grouping observations each compared 128 native output files;
  every observed/unobserved pair was byte-identical (768 file comparisons).
- Five native-CLI known-answer gate regressions passed.
- Six sets of 1,602 source-coordinate audits passed: 9,612 audits, zero overlapping
  or duplicate assignments. Intermediate checks do not certify final graph accounting.
- Within each filtering policy, original and corrected grouping have identical final
  assignment tables and region coordinates after removing only the case label.

## Exact cause of the missing chromosome arms

Under chromosome-isolated filtering and corrected grouping:

| Source | chr7 candidate range, 0-based half-open | Span | Aligned bases | Density | Longest interval | Required interval |
|---|---|---:|---:|---:|---:|---:|
| CTlk hap1 chr7_1+chr12_1 | [0, 32,828,442) | 32.828 Mb | 31,347,168 | 95.4878% | 49,750 bp | strictly >1,000,000 bp |
| CTlk hap2 chr7_1+chr12_1 | [43,248,070, 83,428,940) | 40.181 Mb | 38,201,061 | 95.0728% | 102,554 bp | strictly >1,000,000 bp |

Both density gates pass the configured 25% minimum. Both whole-query adjusted
coverage gates pass (0.929183 and 0.915991). None of these chr7 intervals is removed
by overlap resolution: counts before/after are identical (57,236 and 60,822).
**The longest-individual-interval gate is the sole rejecting local gate in these
corrected runs.** Disabled deletion filtering produces the same chr7 candidates.

The dominant chr12 assignments do not undergo the secondary-run longest-interval
test. Their own longest intervals are only 117,479 and 103,176 bp, yet they remain
assigned. This explains the asymmetry without invoking poor chr7 alignment density.

Two separate issues remain relevant:

1. Native deletion filtering already removes all chr7 input from the hap1 composite
   and leaves only 2,736,233 aligned chr7 bases at the end of the hap2 composite.
   Regional assignment cannot recover mappings already deleted upstream.
2. After chromosome-isolated filtering restores the input, original grouping creates
   56,608 / 60,152 secondary candidate fragments. Corrected grouping produces one
   coherent chr7 run per composite, but that run still fails the explicit individual
   interval-length rule. The comparison correction therefore fixes grouping behavior
   while the independent policy prevents any change to the final assignments.

The result supports the grouping correction as a field-selection defect; it does
not establish that the separate longest-interval requirement is an implementation
bug. Retain that distinction. Assembly cut decisions are unchanged by this diagnostic.

## Next bounded comparison

Use the same cached mappings to compare existing `-u` settings of 1 Mb, 100 kb and
25 kb, across the full selected cohort, retaining original and corrected grouping.
Keep chromosome-isolated filtering as the principal case and the native-filter
control to separate losses at the two stages. No new threshold default or production
graph policy is accepted yet. In particular, changing the upstream condition to use
total run length would be a policy redesign, not another established typo correction.

Acceptance must examine arm recovery together with ordinary chromosome controls,
source-coordinate ownership and new small/discordant chromosome assignments in all
samples, including fragmented assemblies. A lower interval threshold should not be
accepted merely because these two source scaffolds recover their expected arms.

The diagnostic report's aggregate reason can mention both length and density when
different fragments fail different gates. Per-candidate events distinguish these;
the corrected chr7 runs specifically fail length only. Missing input chromosomes
have no candidate row, so the native hap1 absence must be interpreted using the
filtering evidence rather than described as a regional-gate rejection.
