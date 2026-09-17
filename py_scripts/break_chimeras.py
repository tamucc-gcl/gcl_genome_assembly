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
    p.add_argument("--candidates", required=True,
                   help="chimera_joins.py output: one row per CHIMERIC JOIN, with cut_bp "
                        "taken from the AGP. A scaffold may have several -- "
                        "Sde-CTlk_104_hap2 scaffold_3 has three components and one chimeric "
                        "join of its two, and Sde-CPla_115_hap1 scaffold_1 has four -- so "
                        "N joins produce N+1 pieces, not two.")
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
        # candidate_verdict is harmonization's, carried through chimera_joins.py. REVIEW
        # rows are excluded deliberately: that verdict means the vote could not separate an
        # artifact from real biology -- in particular a polymorphic fusion whose sister
        # haplotype also carries the junction -- and those must never be cut unattended.
        cands = [r for r in cands
                 if r.get("candidate_verdict", r.get("verdict")) == "BREAK_CANDIDATE"]

    # index the name map by old_name so a split row can be replaced in place
    by_old = {r["old_name"]: r for r in nm}

    # N joins per scaffold, not one. Grouped, sorted, and validated as a SET: a cut is only
    # legal relative to its neighbours, so each piece has to clear min_piece_bp against the
    # adjacent cuts rather than against the scaffold ends.
    actions, per_scaf = [], {}
    for r in cands:
        sc = r.get("scaffold", "")
        if sc not in by_old:
            actions.append((sc, "SKIP", "not in the name map"))
            continue
        if r.get("callable", "yes") == "no":
            actions.append((sc, "SKIP", "chimera_joins marked it not callable: %s"
                            % r.get("reason", "?")))
            continue
        try:
            cut = int(float(r.get("cut_bp", "") or 0))
        except ValueError:
            actions.append((sc, "SKIP", "cut_bp unparseable"))
            continue
        per_scaf.setdefault(sc, []).append((cut, r.get("left_chrom", ""),
                                            r.get("right_chrom", "")))

    breaks = {}
    for sc, rows in per_scaf.items():
        span = 0
        try:
            span = int(float(by_old[sc].get("length") or 0))
        except ValueError:
            pass
        rows.sort()
        cuts = [c for c, _, _ in rows]
        if any(c <= 0 or (span and c >= span) for c in cuts):
            actions.append((sc, "SKIP", "a cut falls outside 0..%d" % span))
            continue
        # piece lengths between consecutive cuts, including both ends
        bounds = [0] + cuts + ([span] if span else [])
        lens = [b - a2_ for a2_, b in zip(bounds, bounds[1:])]
        if span and any(L < a.min_piece_bp for L in lens):
            actions.append((sc, "SKIP", "cutting at %s would leave a piece < %d bp (%s)"
                            % (",".join(map(str, cuts)), a.min_piece_bp,
                               ",".join(map(str, lens)))))
            continue

        # a name per piece, from the chromosome on that side of each cut. The composite name
        # carries harmonization's part indices; where a chromosome appears more than once
        # among the pieces -- scaffold_3 is chr6 / chr12 / chr6 -- the repeats get distinct
        # suffixes so two records cannot claim the same name.
        toks = (by_old[sc].get("new_name") or "").split("+")
        chroms = [rows[0][1]] + [r[2] for r in rows]
        names, seen = [], {}
        for ch in chroms:
            want = "chr%s_" % ch.replace("chr", "")
            base = next((t for t in toks if t.startswith(want)), want + "1")
            seen[base] = seen.get(base, 0) + 1
            names.append(base if seen[base] == 1 else "%s_%d" % (base, seen[base]))
        breaks[sc] = (cuts, names)
        actions.append((sc, "BREAK", "%s -> %s" % (",".join(map(str, cuts)),
                                                   ",".join(names))))
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
                cuts, _ = breaks[name]
                s = seqs[name]
                # 0-based, open-ended, matching cactus's _sub_X_Y convention
                bounds = [0] + list(cuts) + [len(s)]
                for lo, hi in zip(bounds, bounds[1:]):
                    out.write(">%s_sub_%d_%d\n%s\n" % (name, lo, hi, wrap(s[lo:hi])))
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
                cuts, names = breaks[nmold]
                span = int(float(r.get("length") or 0))
                bounds = [0] + list(cuts) + [span]
                for newname, lo, hi in zip(names, bounds, bounds[1:]):
                    d = dict(r)
                    d["old_name"] = "%s_sub_%d_%d" % (nmold, lo, hi)
                    d["new_name"] = newname
                    d["length"] = str(hi - lo)
                    d["class"] = "chromosome"
                    fl = [x for x in (r.get("flags") or "").split(";")
                          if x and x not in ("-",)]
                    fl.append("chimera_broken(from=%s,at=%s,piece=%d_%d)"
                              % (nmold, ",".join(map(str, cuts)), lo, hi))
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
            out.write("pieces_produced\t%d\n"
                      % sum(len(c) + 1 for c, _ in breaks.values()))
            out.write("cuts_total\t%d\n" % sum(len(c) for c, _ in breaks.values()))
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
