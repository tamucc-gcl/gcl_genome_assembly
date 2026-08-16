#!/usr/bin/env python3
"""
select_reference.py -- choose a harmonization reference from scored candidates.

Consumes the {species}.{id}.score.tsv rows produced by
`harmonize_names.py --score-only`, one per candidate, and writes the winner.

THE RULE, in order
    1. eligible     role == voter. A candidate that failed the voter criteria cannot
                    define a frame for assemblies it is less reliable than.
    2. multiplicity minimise |multiplicity - 1|.

                        multiplicity = placed_scaffolds(other voters)
                                       / (n_consensus x n_other_voters)

                    1.0 is a clean 1:1 chromosome mapping. Above 1, assemblies contribute
                    several scaffolds per chromosome -- a fragmented frame, or redundancy
                    it failed to absorb. Below 1, chromosomes are missing.
    3. edges        fewest join-graph edges required. Zero edges means the graph has
                    nothing to correct, so nothing to get wrong.
    4. gfrac        highest chromosome-set genome fraction.
    5. n50, id      deterministic final tie-breaks.

WHY MULTIPLICITY LEADS AND EDGES ONLY BREAKS TIES
    Edge count alone has a blind spot. A reference that wrongly FUSES two chromosomes needs
    no edges to correct it, so it would score best -- but queries that keep them apart then
    pile onto one piece, and multiplicity rises. Edges detect over-splitting, multiplicity
    detects over-fusion. Neither is sufficient alone.

WHY NOT THE OBVIOUS ONES
    n_chromosomes nearest the batch median: assumes the batch is centred on truth. Two
    failed assemblies shifted the median from 17.5 to 19.5 on this dataset and chose a
    reference split at five junctions.

    Maximise placement: ranks the frames backwards once the reference's own trivially
    placed scaffolds are excluded -- the reference that places the most OTHER scaffolds is
    the one admitting redundant haplotigs into chromosome groups.

    cut_ratio: reported, not used. It measures how cleanly the chromosome set separates
    from the tail, which is a property of the assembly rather than of the frame it induces.

VALIDATED against four references spanning 15 to 26 pieces before being written:
        CMat_203_hap2   mult 0.981  |dev| 0.019  edges  0  gfrac 0.994   selected
        CBau_104_hap2   mult 1.105  |dev| 0.105  edges  5  gfrac 0.989
        CBau_104_hap1   mult 1.124  |dev| 0.124  edges  2  gfrac 0.991
        CTlk_104_hap2   mult 1.163  |dev| 0.163  edges 18  gfrac 0.966
    Frames from the three best were identical (ARI 1.000); only the worst diverged.

Usage:
    select_reference.py --scores *.score.tsv --out-id reference_id.txt \
        --out-table reference_selection.tsv
"""
import argparse
import os
import sys


def read_score(path):
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        line = fh.readline()
        if not line.strip():
            return None
        return dict(zip(hdr, line.rstrip("\n").split("\t")))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", nargs="+", required=True)
    ap.add_argument("--out-id", required=True)
    ap.add_argument("--out-table", required=True)
    ap.add_argument("--require-voter", dest="require_voter", action="store_true",
                    default=True)
    ap.add_argument("--no-require-voter", dest="require_voter", action="store_false",
                    help="allow a non-voter to be selected. Only sensible when every "
                         "candidate failed the voter criteria")
    a = ap.parse_args()

    rows = [r for r in (read_score(p) for p in a.scores) if r]
    if not rows:
        sys.exit("no usable score rows in: %s" % ", ".join(a.scores))

    def num(r, k, d=float("inf")):
        try:
            v = float(r.get(k, ""))
            return d if v != v else v          # NaN -> worst
        except ValueError:
            return d

    pool = [r for r in rows if r.get("role") == "voter"] if a.require_voter else list(rows)
    note = ""
    if not pool:
        pool = list(rows)
        note = ("no candidate passed the voter criteria; selecting among all and flagging "
                "it -- the cohort may be too fragmented for a consensus frame")

    ranked = sorted(pool, key=lambda r: (num(r, "abs_multiplicity_dev"),
                                         num(r, "n_edges"),
                                         -num(r, "genome_fraction", 0.0),
                                         -num(r, "scaffold_n50", 0.0),
                                         r.get("reference_id", "")))
    win = ranked[0]

    # is the choice comfortable, or a coin flip?
    amb = ""
    if len(ranked) > 1:
        d0 = num(ranked[0], "abs_multiplicity_dev")
        d1 = num(ranked[1], "abs_multiplicity_dev")
        if d1 - d0 < 0.02:
            amb = ("multiplicity separates the top two by only %.3f; the choice rests on "
                   "the edge-count tie-break" % (d1 - d0))

    with open(a.out_id, "w") as fh:
        fh.write(win["reference_id"] + "\n")

    cols = ["reference_id", "role", "selected", "rank", "abs_multiplicity_dev",
            "multiplicity", "n_edges", "n_consensus_chromosomes", "n_reference_pieces",
            "genome_fraction", "cut_ratio", "scaffold_n50", "chrom_set_method",
            "chrom_set_flags", "graph_rule", "n_other_voters", "placed_other_voters"]
    order = {r["reference_id"]: i + 1 for i, r in enumerate(ranked)}
    with open(a.out_table, "w") as fh:
        fh.write("# rule: voter -> min |multiplicity-1| -> fewest edges -> "
                 "max genome_fraction -> max n50\n")
        if note:
            fh.write("# WARNING\t%s\n" % note)
        if amb:
            fh.write("# NOTE\t%s\n" % amb)
        fh.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: order.get(r["reference_id"], 999)):
            r = dict(r)
            r["selected"] = "1" if r["reference_id"] == win["reference_id"] else "0"
            r["rank"] = str(order.get(r["reference_id"], ""))
            fh.write("\t".join(r.get(c, "NA") for c in cols) + "\n")

    if note:
        sys.stderr.write("[select-reference] WARNING: %s\n" % note)
    if amb:
        sys.stderr.write("[select-reference] NOTE: %s\n" % amb)
    sys.stderr.write(
        "[select-reference] %s of %d candidates: |mult-1|=%s edges=%s gfrac=%s "
        "consensus=%s\n" % (win["reference_id"], len(rows),
                            win.get("abs_multiplicity_dev"), win.get("n_edges"),
                            win.get("genome_fraction"),
                            win.get("n_consensus_chromosomes")))
    for i, r in enumerate(ranked[1:5], start=2):
        sys.stderr.write("[select-reference]   %d. %-22s |mult-1|=%s edges=%s gfrac=%s\n"
                         % (i, r["reference_id"], r.get("abs_multiplicity_dev"),
                            r.get("n_edges"), r.get("genome_fraction")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
