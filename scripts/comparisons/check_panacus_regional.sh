#!/usr/bin/env bash
# Tiny tool-behaviour probe, not a biological analysis. Run in Panacus 0.5.2 environment.
set -euo pipefail
out=${1:?Supply an unused output directory under comparisons/}
mkdir "$out"
cd "$out"
cat > fixture.gfa <<'GFA'
H	VN:Z:1.0
S	1	AAAA
S	2	CC
S	3	GGG
L	1	+	2	+	0M
L	1	+	3	+	0M
P	A#1#chr1	1+,2+	*
P	A#1#chr2	3+	*
P	A#2#chr1	1+,3+	*
P	B#1#chr1	1+	*
GFA
printf 'A#1#chr1\tA1\nA#1#chr2\tA1\nA#2#chr1\tA2\nB#1#chr1\tB1\n' > groups.tsv
printf 'A#1#chr1\n' > region.txt
panacus --version > version.txt
panacus -t 1 hist -c bp -g groups.tsv fixture.gfa > global.hist.tsv
panacus -t 1 hist -c bp -g groups.tsv -s region.txt fixture.gfa > subset.hist.tsv
panacus -t 1 table -c bp -g groups.tsv fixture.gfa > membership.tsv
panacus -t 1 similarity -c bp -g groups.tsv fixture.gfa > similarity.tsv
cat > expected.md <<'EXPECTED'
# Hand-calculated fixture expectations

Node 1: 4 bp, present in A1/A2/B1, coverage 3.
Node 2: 2 bp, present only in A1, coverage 1.
Node 3: 3 bp, present in A1 (chr2) and A2 (chr1), coverage 2.
Global histogram: coverage 0 = 0 bp, 1 = 2 bp, 2 = 3 bp, 3 = 4 bp.
A1 whole haplotype: private 2 bp, coverage-two 3 bp, core 4 bp; distinct total 9 bp.
A1 chr1 with GLOBAL coverage preserved: private 2 bp, core 4 bp; total 6 bp.
If subset.hist instead reports 6 bp at coverage 1, subsetting changed the counting population and cannot directly supply the required attribution.

The cross-chromosome node 3 deliberately tests why chromosome totals need not be additive for distinct graph sequence.
Expected bp-weighted Jaccard similarities: A1/A2 = 7/9, A1/B1 = 4/9, A2/B1 = 4/7; diagonal 1.
Inspect table units and exported labels before implementing the adapter or interpreting similarity as a distance.
EXPECTED
printf 'Probe files written to %s\n' "$PWD"
