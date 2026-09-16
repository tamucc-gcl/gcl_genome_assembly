#!/usr/bin/env python3
# ======================================================================================
# break_chimeras.py
#
# Splits chimeric scaffolds at the junction and rewrites the harmonization name map so the
# two halves get the chromosome names they deserve.
#
# WHERE THIS SITS AND WHY
# -----------------------
# Between HARMONIZE_SCAFFOLDS and FINALIZE_ASSEMBLY. Harmonization emits a NAME MAP, not a
# renamed FASTA -- FINALIZE_ASSEMBLY applies the map. So this is the last point where the
# FASTA is still in its original coordinates and the composite is still one record, and the
# first point where the concordance vote exists to justify cutting it.
#
# It rewrites BOTH:
#   the FASTA -- one record becomes two, named <scaffold>_sub_<start>_<end>
#   the name map -- one row becomes two, mapping those to chr5_1 and chr9_1
#
# so FINALIZE_ASSEMBLY works unchanged.
#
# THE _sub_X_Y CONVENTION IS CACTUS'S, AND IT IS USED HERE FOR A DIFFERENT REASON
# ------------------------------------------------------------------------------
# The cactus manual uses `<name>_sub_<start>_<end>` (0-based, open-ended like BED) for
# manually broken scaffolds, and preserves the offsets into the GFA W-lines. Our scaffolds get
# renamed by harmonization before cactus ever sees them, so that particular benefit does not
# apply. The convention is kept anyway because it makes the name-map rewrite mechanical and
# self-documenting: `scaffold_1_sub_0_37000000 -> chr5_1` states the provenance of the piece
# in the file that already records every rename.
#
# THE CASE THIS IS FOR
# --------------------
#   Sde-CTlk_104_hap1  scaffold_1  chr5_1+chr9_1  111,642,300 bp  vote 1f/7s
#     junction ~37 Mb; ref5 footprint 36.0 Mb, ref9 footprint 73.6 Mb
#
# cactus-graphmap-split assigns each contig to a SINGLE chromosome, so unsplit this scaffold
# put ~73 Mb of chromosome 9 into chr5's subgraph where nothing aligns to it -- absent from
# every graph built so far, and the reason chr9 has no core sequence.
#
# WHICH HALF GETS WHICH NAME comes from the junction finder's left_member/right_member, NOT
# from the order in the composite name. `chr5_1+chr9_1` is sorted by chromosome number; the
# scaffold could run chr9-then-chr5.
#
# USAGE
#   break_chimeras.py --fasta in.fa --name-map in.tsv --candidates cand.tsv \
#       --assembly <id> --out-fasta out.fa --out-name-map out.tsv [--audit a.tsv]
#       [--mode auto|file]
#
#   mode=auto  break rows whose verdict is BREAK_CANDIDATE
#   mode=file  break every row present (the file IS the instruction; edit it to choose)
# ======================================================================================

import argparse
import os
import sys


def read_rows(path):
    """'#'-commented TSV -> list of dicts. Comments stripped by hand and comment.char kept
    OFF: '#' is legitimate inside PanSN names elsewhere in this pipeline, and using it as a
    comment marker has truncated keys twice."""
    rows, hdr = [], None
    if not path or not os.path.isfile(path) or os.path.getsize(path) == 0:
        return rows
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.lstrip().startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = f
                continue
            if len(f) == len(hdr):
                rows.append(dict(zip(hdr, f)))
    return rows


def read_fasta(path):
    """name -> sequence. Whole-file, because a split needs random access to one record and
    these are ~1 Gb assemblies read once."""
    seqs, name, buf = {}, None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def wrap(s, w=60):
    return "\n".join(s[i:i + w] for i in range(0, len(s), w))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fasta", required=True)
    p.add_argument("--name-map", required=True)
    p.add_argument("--candidates", required=True)
    p.add_argument("--assembly", required=True)
    p.add_argument("--out-fasta", required=True)
    p.add_argument("--out-name-map", required=True)
    p.add_argument("--audit", default="")
    p.add_argument("--mode", choices=("auto", "file"), default="auto")
    p.add_argument("--min-piece-bp", type=int, default=1000000,
                   help="refuse a cut that would leave a piece smaller than this (default "
                        "1 Mb). A junction estimated from 1 Mb bins can land near an end; "
                        "cutting there produces a fragment rather than a chromosome.")
    a = p.parse_args()

    nm = read_rows(a.name_map)
    if not nm:
        sys.exit("ERROR: no usable rows in %s" % a.name_map)
    cands = [r for r in read_rows(a.candidates)
             if r.get("assembly") == a.assembly]
    if a.mode == "auto":
        cands = [r for r in cands if r.get("verdict") == "BREAK_CANDIDATE"]

    # index the name map by old_name so a split row can be replaced in place
    by_old = {r["old_name"]: r for r in nm}

    actions = []
    for r in cands:
        sc = r.get("scaffold", "")
        if sc not in by_old:
            actions.append((sc, "SKIP", "not in the name map"))
            continue
        try:
            j = int(float(r.get("junction_bp", "") or 0))
            span = int(float(r.get("span_bp", "") or 0))
        except ValueError:
            actions.append((sc, "SKIP", "junction_bp or span_bp unparseable"))
            continue
        if j <= 0 or j >= span:
            actions.append((sc, "SKIP", "junction %d outside 0..%d" % (j, span)))
            continue
        if j < a.min_piece_bp or (span - j) < a.min_piece_bp:
            actions.append((sc, "SKIP", "a piece would be < %d bp" % a.min_piece_bp))
            continue
        lm, rm = r.get("left_member", ""), r.get("right_member", "")
        if not (lm.startswith("ref") and rm.startswith("ref")):
            actions.append((sc, "SKIP", "left/right member missing"))
            continue
        # the composite name carries the part indices harmonization assigned; take each
        # half's name from the token whose chromosome matches that side
        toks = (by_old[sc].get("new_name") or "").split("+")
        want_l, want_r = "chr%s_" % lm[3:], "chr%s_" % rm[3:]
        nl = next((t for t in toks if t.startswith(want_l)), want_l + "1")
        nr = next((t for t in toks if t.startswith(want_r)), want_r + "1")
        actions.append((sc, "BREAK", "%d|%s,%s" % (j, nl, nr)))

    breaks = {s: v.split("|") for s, k, v in actions if k == "BREAK"}
    if not breaks:
        sys.stderr.write("[break_chimeras] %s: nothing to break (%d candidate rows)\n"
                         % (a.assembly, len(cands)))

    seqs = read_fasta(a.fasta)
    order = [r["old_name"] for r in nm]

    n_break = 0
    with open(a.out_fasta, "w") as out:
        for name in order:
            if name not in seqs:
                continue
            if name in breaks:
                j = int(breaks[name][0])
                s = seqs[name]
                # 0-based, open-ended, matching the _sub_X_Y convention
                out.write(">%s_sub_0_%d\n%s\n" % (name, j, wrap(s[:j])))
                out.write(">%s_sub_%d_%d\n%s\n" % (name, j, len(s), wrap(s[j:])))
                n_break += 1
            else:
                out.write(">%s\n%s\n" % (name, wrap(seqs[name])))
        # anything in the FASTA but not the map (should not happen) is carried through
        for name in seqs:
            if name not in by_old:
                out.write(">%s\n%s\n" % (name, wrap(seqs[name])))
                sys.stderr.write("[break_chimeras] WARNING: %s is in the FASTA but not the "
                                 "name map; carried through unrenamed\n" % name)

    # ---- rewrite the name map: one row becomes two ----------------------------------
    with open(a.out_name_map, "w") as out:
        hdr = list(nm[0].keys())
        out.write("\t".join(hdr) + "\n")
        for r in nm:
            nmold = r["old_name"]
            if nmold in breaks:
                j, names = breaks[nmold]
                nl, nr = names.split(",")
                span = int(float(r.get("length") or 0))
                for newname, piece, lo, hi in ((nl, "_sub_0_%s" % j, 0, int(j)),
                                               (nr, "_sub_%s_%d" % (j, span), int(j), span)):
                    d = dict(r)
                    d["old_name"] = nmold + piece
                    d["new_name"] = newname
                    d["length"] = str(hi - lo)
                    d["class"] = "chromosome"
                    fl = [x for x in (r.get("flags") or "").split(";")
                          if x and x not in ("-",)]
                    fl.append("chimera_broken(from=%s,at=%s)" % (nmold, j))
                    d["flags"] = ";".join(fl)
                    out.write("\t".join(str(d.get(k, "")) for k in hdr) + "\n")
            else:
                out.write("\t".join(str(r.get(k, "")) for k in hdr) + "\n")

    if a.audit:
        with open(a.audit, "w") as out:
            out.write("# What break_chimeras.py did, and to what. A break that is not\n")
            out.write("#   written down is a break nobody can audit.\n")
            out.write("# mode=file means the candidates file IS the instruction -- every row\n")
            out.write("#   in it is broken. mode=auto breaks only verdict=BREAK_CANDIDATE.\n")
            out.write("metric\tvalue\n")
            out.write("assembly\t%s\n" % a.assembly)
            out.write("mode\t%s\n" % a.mode)
            out.write("candidate_rows\t%d\n" % len(cands))
            out.write("scaffolds_broken\t%d\n" % n_break)
            out.write("min_piece_bp\t%d\n" % a.min_piece_bp)
            for sc, kind, why in actions:
                out.write("action.%s\t%s:%s\n" % (sc, kind, why))

    sys.stderr.write("[break_chimeras] %s: broke %d scaffold(s) of %d candidate row(s)\n"
                     % (a.assembly, n_break, len(cands)))
    for sc, kind, why in actions:
        if kind == "SKIP":
            sys.stderr.write("[break_chimeras]   SKIP %s -- %s\n" % (sc, why))


if __name__ == "__main__":
    main()
