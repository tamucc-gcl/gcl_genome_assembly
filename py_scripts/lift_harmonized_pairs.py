#!/usr/bin/env python3
"""
lift_harmonized_pairs.py -- translate Hi-C pairs from pre-harmonization scaffold
coordinates into finalized (harmonized) coordinates.

WHY THIS IS EXACT
-----------------
FINALIZE_ASSEMBLY on the harmonized path performs exactly three operations: rename
old_name -> new_name, reverse-complement where orient == 'rev', and reorder. It does not
filter (the min_scaffold_bp threshold applies only on the NO_HARMONIZE size-rank path) and
does not alter sequence. Scaffold set and lengths are invariant.

A pairs file on pre-harmonization coordinates therefore converts by arithmetic alone:

    chrom  : old_name -> new_name
    pos    : if rev,  pos' = L - pos + 1        (1-based, so no off-by-one)
    strand : if rev,  '+' <-> '-'
    order  : restore the upper-triangular convention after renaming, because the new
             chromosome order is not the old one

No remapping, no AGP, no gap/offset bookkeeping. FINALIZE writes the applied-lift map on
BOTH naming paths -- harmonized (orient may be 'rev') and size-rank (orient always 'fwd')
-- so this is one code path with no options.

WHAT IT DOES NOT HANDLE
-----------------------
Any transformation that changes sequence content or coordinates -- gap filling, telomere
extension, Inspector breaks, decontamination removal. Those all happen UPSTREAM of the
assembly this lift starts from. If you point this at a pairs file mapped to contigs or to a
pre-gap-fill assembly it would produce silently wrong coordinates, so the pairs header
chromsizes are cross-checked against the map lengths and a mismatch is fatal. There is no
flag to bypass that; there is no setting where continuing is correct.

INPUTS
------
--pairs      pairs.gz (or plain) in pre-harmonization coordinates
--name-map   <id>.applied_lift.tsv from FINALIZE_ASSEMBLY
             cols: old_name new_name orient order length

OUTPUTS
-------
--out        translated pairs, UNSORTED, uncompressed (pipe to pairtools sort)
--chrom-sizes  new-name chrom.sizes in finalized order
--stats        TSV of counts, for the validator and for QC

Usage:
    lift_harmonized_pairs.py --pairs in.pairs.gz --name-map nm.tsv \
        --out lifted.pairs --chrom-sizes final.chrom.sizes --stats lift_stats.tsv
"""
import argparse
import gzip
import sys


def openf(path, mode="rt"):
    if path == "-":
        return sys.stdin if "r" in mode else sys.stdout
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def read_name_map(path, invert=False):
    """Read FINALIZE_ASSEMBLY's applied-lift map.

    Columns: old_name  new_name  orient  order  length

    FINALIZE writes this on both naming paths -- harmonized (orient may be 'rev') and
    size-rank (orient always 'fwd') -- so there is exactly one schema and no branch here.
    `order` is the actual output order of the finalized FASTA, dense 1..N.

    invert=True reads the map backwards: new -> old, same orient (reverse-complement is
    its own inverse), order taken as file row order. Used only by the round-trip check.

    Returns (lift, ordered) with lift[src] = (dst, is_rev, length).
    """
    lift, rows = {}, []
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        need = ("old_name", "new_name", "orient", "order", "length")
        missing = [k for k in need if k not in hdr]
        if missing:
            sys.exit(f"applied-lift map {path} missing column(s) {missing}; got {hdr}.\n"
                     f"Expected FINALIZE_ASSEMBLY's <id>.applied_lift.tsv, not the "
                     f"4-column <id>.name_map.tsv (which has no orient).")
        I = {k: hdr.index(k) for k in need}
        for rn, line in enumerate(fh, start=1):
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            a, b = f[I["old_name"]], f[I["new_name"]]
            src, dst = (b, a) if invert else (a, b)
            rev = f[I["orient"]].strip().lower() == "rev"
            order = rn if invert else int(f[I["order"]])
            length = int(f[I["length"]])
            if src in lift:
                sys.exit(f"duplicate source name {src!r} in {path}")
            lift[src] = (dst, rev, length)
            rows.append((order, dst, length))
    rows.sort()
    return lift, [(n, L) for _o, n, L in rows]


def parse_header(lines):
    """Split a 4DN pairs header into the pieces we keep, rewrite, or drop."""
    fmt, columns, keep, chromsize = None, None, [], {}
    for ln in lines:
        s = ln.rstrip("\n")
        low = s.lower()
        if low.startswith("## pairs format"):
            fmt = s
        elif low.startswith("#columns:"):
            columns = s
        elif low.startswith("#chromsize:"):
            parts = s.split(":", 1)[1].split()
            if len(parts) >= 2:
                chromsize[parts[0]] = int(parts[1])
        elif low.startswith("#samheader:") and "\t@sq\t" in s.replace(" ", "\t").lower():
            pass                      # stale @SQ records -- drop, they name old scaffolds
        elif low.startswith("#sorted:") or low.startswith("#shape:"):
            pass                      # re-established below / by pairtools sort
        else:
            keep.append(s)
    return fmt, columns, keep, chromsize


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--name-map", required=True)
    ap.add_argument("--out", required=True, help="translated pairs, unsorted")
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--stats", required=True)
    ap.add_argument("--sample-id", default="")
    ap.add_argument("--invert", action="store_true", default=False,
                    help="apply the map backwards (finalized -> pre-finalize). Only for the "
                         "round-trip check: invert then forward must recover the input "
                         "exactly, which validates the arithmetic on real data")
    a = ap.parse_args()

    lift, ordered = read_name_map(a.name_map, invert=a.invert)
    new_idx = {n: i for i, (n, _L) in enumerate(ordered)}

    with open(a.chrom_sizes, "w") as fh:
        for n, L in ordered:
            fh.write(f"{n}\t{L}\n")

    fin = openf(a.pairs)
    header_lines = []
    first_body = None
    for line in fin:
        if line.startswith("#"):
            header_lines.append(line)
        else:
            first_body = line
            break
    fmt, columns, keep, chromsize = parse_header(header_lines)

    # ---- consistency gate: the pairs must be on the assembly the name map describes
    problems = []
    if chromsize:
        missing = sorted(set(chromsize) - set(lift))
        if missing:
            problems.append("%d scaffold(s) in the pairs header are absent from the name "
                            "map, e.g. %s" % (len(missing), ", ".join(missing[:5])))
        badlen = [(c, chromsize[c], lift[c][2]) for c in chromsize
                  if c in lift and chromsize[c] != lift[c][2]]
        if badlen:
            problems.append("%d length mismatch(es), e.g. %s" % (
                len(badlen), "; ".join("%s pairs=%d map=%d" % b for b in badlen[:3])))
    if problems:
        sys.exit("[lift] REFUSING: " + " | ".join(problems)
                 + "\n[lift] These pairs are not on the assembly this name map describes. "
                   "Lifting them would produce silently wrong coordinates. Check that the "
                   "Hi-C was mapped to the same assembly FINALIZE renamed.")

    # ---- header out (#columns must come last)
    out = openf(a.out, "wt")
    out.write((fmt or "## pairs format v1.0") + "\n")
    out.write("#shape: upper triangle\n")
    if a.sample_id:
        out.write(f"#genome_assembly: {a.sample_id}\n")
    out.write("#lift: harmonized rename+revcomp from %s via %s\n"
              % (a.pairs, a.name_map))
    for n, L in ordered:
        out.write(f"#chromsize: {n} {L}\n")
    for s in keep:
        out.write(s + "\n")
    out.write((columns or
               "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type") + "\n")

    # ---- body: tight loop
    n_in = n_out = n_drop = n_rev1 = n_rev2 = n_swap = 0
    FLIP = {"+": "-", "-": "+"}
    write = out.write
    lget = lift.get
    body_iter = fin if first_body is None else _chain(first_body, fin)
    for line in body_iter:
        n_in += 1
        f = line.rstrip("\n").split("\t")
        e1 = lget(f[1])
        e2 = lget(f[3])
        if e1 is None or e2 is None:
            sys.exit("[lift] scaffold %r or %r absent from the applied-lift map (body line "
                     "%d). FINALIZE keeps every scaffold, so this means the pairs and the "
                     "map are from different assemblies." % (f[1], f[3], n_in))
        n1, r1, L1 = e1
        n2, r2, L2 = e2
        p1 = int(f[2])
        p2 = int(f[4])
        s1 = f[5]
        s2 = f[6]
        if r1:
            p1 = L1 - p1 + 1
            s1 = FLIP.get(s1, s1)
            n_rev1 += 1
        if r2:
            p2 = L2 - p2 + 1
            s2 = FLIP.get(s2, s2)
            n_rev2 += 1
        # restore upper-triangular order in the NEW chromosome ordering
        i1 = new_idx[n1]
        i2 = new_idx[n2]
        if i1 > i2 or (i1 == i2 and p1 > p2):
            n1, n2 = n2, n1
            p1, p2 = p2, p1
            s1, s2 = s2, s1
            n_swap += 1
        f[1] = n1
        f[2] = str(p1)
        f[3] = n2
        f[4] = str(p2)
        f[5] = s1
        f[6] = s2
        write("\t".join(f) + "\n")
        n_out += 1

    if out is not sys.stdout:
        out.close()
    fin.close()

    with open(a.stats, "w") as fh:
        fh.write("metric\tvalue\n")
        for k, v in (("pairs_in", n_in), ("pairs_out", n_out), ("pairs_dropped", n_drop),
                     ("ends_flipped_side1", n_rev1), ("ends_flipped_side2", n_rev2),
                     ("pairs_reordered", n_swap),
                     ("scaffolds", len(ordered)),
                     ("scaffolds_reversed", sum(1 for o in lift if lift[o][1])),
                     ("inverted", int(a.invert))):
            fh.write(f"{k}\t{v}\n")
    sys.stderr.write(
        "[lift] %s: in=%d out=%d dropped=%d flipped=%d/%d reordered=%d "
        "scaffolds=%d (reversed=%d)\n" % (
            a.sample_id or a.pairs, n_in, n_out, n_drop, n_rev1, n_rev2, n_swap,
            len(ordered), sum(1 for o in lift if lift[o][1])))


def _chain(first, rest):
    yield first
    for x in rest:
        yield x


if __name__ == "__main__":
    main()
