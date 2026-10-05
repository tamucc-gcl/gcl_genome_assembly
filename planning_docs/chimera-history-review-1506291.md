# Historical assembly comparison: job 1506291

Run SUCCESS; all four current/historical haplotype combinations completed. Current
assessment hashes agree with the origin diagnostic. Interpretation below is from
the archived PAF alignments, not an independent comparison of FASTA letters.

## Five raw-contig transitions predate the current scaffolding

| Current candidate | Historical match across transition | Exact aligned extent |
|---|---|---|
| hap1 scaffold_1, ~62.780 Mb | old hap2 chr7_2, reverse | full 523,749 bp window |
| hap1 scaffold_1, ~95.963 Mb | old hap2 chr12_1, reverse | full 369,712 bp window |
| hap1 scaffold_5, ~35.513 Mb | old hap1 chr15_2, reverse | 303,862 bp of 337,672 bp window; includes full uncertainty interval and both anchors |
| hap2 scaffold_5, ~43.144 Mb | old hap1 chr12_1, forward | full 378,198 bp window |
| hap2 scaffold_7, ~39.619 Mb | old hap2 chr15_1, forward | full 338,862 bp window |

These alignments have only `=` CIGAR operations within the stated extents.
Historical recurrence is not independent biological validation: assemblies may
share reads, contigs and systematic errors. Haplotype correspondence varies by
region, so a global hap1/hap2 relabelling would be inappropriate.

## Current scaffold composition differs from historical scaffold composition

Dominant historical matches change elsewhere on these scaffolds. Examples of
adjacent long alignment endpoints (0-based, half-open current coordinates):

| Current scaffold | Historical comparison | Endpoint interval | Historical targets |
|---|---|---|---|
| hap1 scaffold_1 | old hap2 | 63,051,203–63,051,334 | chr7_2 -> chr6_3 |
| hap2 scaffold_5 | old hap1 | 54,123,484–54,123,676 | chr12_1 -> chr6_2 |
| hap2 scaffold_7 | old hap2 | 46,160,462–46,160,583 | chr15_1 -> chr4_2 |

These are historical target-switch leads, NOT inferred gap coordinates or
authorized cuts. The table is illustrative, not an exhaustive switch detector.
Historical fragmentation can generate this pattern even for correct current
joins. Match switches must be intersected with verified AGP boundaries and
checked against all substantial placements, local sequence and read evidence.
For hap1 scaffold_1, another major chr6_3/chr12_1 transition lies between long
blocks ending at 94,151,483 and starting at 95,221,878; smaller alignments must
be reviewed before narrowing it. Names in different assembly runs are not
assumed to denote the same homologous chromosome.

## Two joins around the 113,941 bp component

Current hap1 scaffold_1 component h1tg002235l_1 lies at
116,057,547–116,171,488, between two 100-N gaps established by the origin run.
The history comparison does not show a continuous historical alignment through
the entire inserted component and its two neighbours. Parts match unplaced_122
in old hap2; upstream sequence has a 164,969 bp exact match to old hap2 chr12_1,
and downstream sequence includes a 100,694 bp exact alignment plus a 1 bp deletion
against old hap1 chr12_1. Different local haplotype affinity is a lead, not proof
of a phase switch. Small repetitive matches and missing exact anchors cannot
establish loss or incorrect placement of the component.

## Next investigation

Expand the YaHS-stage review to all historical target switches on the four
suspect scaffolds, alongside the two already verified component joins. The
current chromosome-label transition list alone does not locate every possible
scaffolding error. Establish actual component/gap coordinates first. Evaluate
current versus alternative adjacency with library-specific Hi-C and local HiFi
opportunity, and investigate possible phase inconsistency. Then revisit the five
raw-contig assignment transitions with the historical recurrence recorded.
No cuts or biological fusion declarations are supported by this comparison alone.
