#!/usr/bin/env python3
"""
select_reference.py -- choose a harmonization reference from scored candidates.

Consumes the {species}.{id}.score.tsv rows produced by
`harmonize_names.py --score-only`, one per candidate, and writes the winner.

THE RULE, in order
    1. eligible     role == voter. A candidate that failed the voter criteria cannot
                    define a frame for assemblies it is less reliable than.
    2. frame        the frame a candidate induces: (n_consensus_chromosomes, n_edges).
                    Candidates are grouped by it, and the group with the FEWEST EDGES at the
                    MODAL chromosome count wins. Zero edges means the graph had nothing to
                    correct, so nothing to get wrong.
    2b. gate       reject |multiplicity - 1| > --max-multiplicity-dev outright. See below.
    3. within group Candidates that induce the SAME frame are interchangeable as far as
                    chromosome identity goes -- the frame experiment measured ARI 1.000
                    between references of 15, 17 and 20 pieces. So the discriminator is
                    placement power: highest chromosome-set genome_fraction.
    4. multiplicity |multiplicity - 1|, only among candidates that also tie on
                    genome_fraction.
    5. n50, id      deterministic final tie-breaks.

WHY MULTIPLICITY IS NO LONGER THE PRIMARY -- A REAL MISTAKE, CAUGHT ON REAL DATA
    It was, and it chose badly. Scoring eight candidates gave four tied at 15 pieces /
    15 chromosomes / 0 edges:

        CLim_110_hap1   |mult-1| 0.0381   placed 109   gfrac 0.9194   <- chosen by mult
        CMat_203_hap2   |mult-1| 0.0667   placed 112   gfrac 0.9944
        CLim_110_hap2   |mult-1| 0.0762   placed 113   gfrac 0.9207
        CMat_203_hap1   |mult-1| 0.0857   placed 114   gfrac 0.9905

    Multiplicity picked the candidate placing the FEWEST query scaffolds, and the one whose
    own chromosome set holds only 91.9% of its sequence against CMat's 99.4%. A reference
    with 8% of itself outside the chromosome set gives queries aligning to that 8% nowhere
    to land, so the low multiplicity is UNDER-PLACEMENT dressed up as a clean 1:1 mapping.

    Multiplicity is still the right instrument for telling a good frame from a fragmented
    one -- CTlk_104_hap2 at 1.16 with 18 edges is correctly last, and CTlk_104_hap1 at 1.24
    with 13 chromosomes worse still. It is the wrong instrument for choosing among
    candidates whose frames are already identical, which is what the earlier four-candidate
    validation could not reveal because there multiplicity, edges and genome fraction all
    agreed.

WHY MULTIPLICITY IS A GATE AND NOT ONLY A TIE-BREAK
    The modal chromosome count catches a fusing reference only when its peers disagree with
    it. Testing found the case where they do not: two candidates both reporting 13
    chromosomes, one of them fusing. The mode cannot discriminate, edges favours the fusing
    one (nothing to merge), and genome fraction favours it too -- so it won outright with a
    multiplicity of 1.65, meaning every other assembly piled 65% extra scaffolds onto its
    chromosomes.

    So multiplicity is applied first as a rejection threshold, and only afterwards as the
    last tie-break. It is the one criterion that sees over-fusion at all.

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
    ap.add_argument("--max-multiplicity-dev", type=float, default=0.30,
                    help="reject a candidate whose |multiplicity-1| exceeds this, before "
                         "any other criterion. Multiplicity has to be a GATE and not just "
                         "a tie-break: a reference that FUSES two real chromosomes needs no "
                         "graph edges to correct it and can carry the best genome fraction, "
                         "so ordering by edges then gfrac would select it outright whenever "
                         "its (wrong) chromosome count happens to be modal. Queries keeping "
                         "those chromosomes apart pile onto one piece, which is what a high "
                         "multiplicity detects")
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

    # ---- multiplicity gate, applied BEFORE anything else.
    # A fusing reference is invisible to every other criterion: it needs no edges (there is
    # nothing to merge), it can have the highest genome fraction, and if its wrong
    # chromosome count is the modal one the mode guard cannot fire either. The only signal
    # is that queries which keep those chromosomes apart pile onto one of its pieces.
    gated = [r for r in pool if num(r, "abs_multiplicity_dev") <= a.max_multiplicity_dev]
    dropped = [r for r in pool if r not in gated]
    if gated:
        pool = gated
    elif pool:
        note = ("every candidate exceeded the multiplicity gate (|mult-1| > %.2f); "
                "selecting among all and flagging it -- no candidate maps 1:1 onto the "
                "cohort, so the frame is suspect whichever is used" % a.max_multiplicity_dev)

    if not pool:
        pool = list(rows)
        note = ("no candidate passed the voter criteria; selecting among all and flagging "
                "it -- the cohort may be too fragmented for a consensus frame")

    # ---- the modal chromosome count. A candidate reporting fewer chromosomes than its
    # peers is either fusing real chromosomes or over-merging in the graph; either way its
    # frame disagrees with the consensus and it should not define it.
    counts = {}
    for r in pool:
        c = r.get("n_consensus_chromosomes")
        if c not in (None, "", "NA"):
            counts[c] = counts.get(c, 0) + 1
    modal = max(counts, key=lambda c: (counts[c], -abs(num({"x": c}, "x", 0)))) if counts \
        else None

    def key(r):
        return (
            0 if r.get("n_consensus_chromosomes") == modal else 1,   # modal frame first
            num(r, "n_edges"),                                       # then fewest edges
            -num(r, "genome_fraction", 0.0),                         # then placement power
            num(r, "abs_multiplicity_dev"),                          # then mapping cleanliness
            -num(r, "scaffold_n50", 0.0),
            r.get("reference_id", ""),
        )

    ranked = sorted(pool, key=key)
    win = ranked[0]

    # how many candidates induce the SAME frame as the winner? If several do, the choice
    # among them was made on placement power, and any of them would give the same
    # chromosome identities -- worth saying so rather than implying the winner was unique.
    same_frame = [r for r in pool
                  if r.get("n_consensus_chromosomes") == win.get("n_consensus_chromosomes")
                  and r.get("n_edges") == win.get("n_edges")]
    amb = ""
    if len(same_frame) > 1:
        amb = ("%d candidates induce the identical frame (%s chromosomes, %s edges): %s. "
               "Chromosome identity is the same whichever is used; the winner was chosen "
               "on genome_fraction."
               % (len(same_frame), win.get("n_consensus_chromosomes"), win.get("n_edges"),
                  ", ".join(sorted(r["reference_id"] for r in same_frame))))

    with open(a.out_id, "w") as fh:
        fh.write(win["reference_id"] + "\n")

    cols = ["reference_id", "role", "selected", "rank", "modal_frame",
            "n_consensus_chromosomes", "n_edges", "genome_fraction",
            "abs_multiplicity_dev", "multiplicity", "n_reference_pieces",
            "cut_ratio", "scaffold_n50", "chrom_set_method",
            "chrom_set_flags", "graph_rule", "n_other_voters", "placed_other_voters"]
    order = {r["reference_id"]: i + 1 for i, r in enumerate(ranked)}
    with open(a.out_table, "w") as fh:
        fh.write("# rule: voter -> modal n_consensus_chromosomes -> fewest edges -> "
                 "max genome_fraction -> min |multiplicity-1| -> max n50\n")
        fh.write("# modal n_consensus_chromosomes across candidates: %s\n" % modal)
        fh.write("# multiplicity gate: |mult-1| <= %.2f\n" % a.max_multiplicity_dev)
        for r in dropped:
            fh.write("# GATED\t%s\t|mult-1|=%s\tlikely fusing real chromosomes\n"
                     % (r["reference_id"], r.get("abs_multiplicity_dev")))
        if note:
            fh.write("# WARNING\t%s\n" % note)
        if amb:
            fh.write("# NOTE\t%s\n" % amb)
        fh.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: order.get(r["reference_id"], 999)):
            r = dict(r)
            r["selected"] = "1" if r["reference_id"] == win["reference_id"] else "0"
            r["rank"] = str(order.get(r["reference_id"], ""))
            r["modal_frame"] = "1" if r.get("n_consensus_chromosomes") == modal else "0"
            fh.write("\t".join(r.get(c, "NA") for c in cols) + "\n")

    for r in dropped:
        sys.stderr.write("[select-reference] GATED %s: |mult-1|=%s exceeds %.2f -- likely "
                         "fusing real chromosomes\n"
                         % (r["reference_id"], r.get("abs_multiplicity_dev"),
                            a.max_multiplicity_dev))
    if note:
        sys.stderr.write("[select-reference] WARNING: %s\n" % note)
    if amb:
        sys.stderr.write("[select-reference] NOTE: %s\n" % amb)
    sys.stderr.write(
        "[select-reference] %s of %d candidates: consensus=%s edges=%s gfrac=%s "
        "|mult-1|=%s\n" % (win["reference_id"], len(rows),
                           win.get("n_consensus_chromosomes"), win.get("n_edges"),
                           win.get("genome_fraction"), win.get("abs_multiplicity_dev")))
    for i, r in enumerate(ranked[1:5], start=2):
        sys.stderr.write("[select-reference]   %d. %-22s consensus=%s edges=%s gfrac=%s "
                         "|mult-1|=%s\n"
                         % (i, r["reference_id"], r.get("n_consensus_chromosomes"),
                            r.get("n_edges"), r.get("genome_fraction"),
                            r.get("abs_multiplicity_dev")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
