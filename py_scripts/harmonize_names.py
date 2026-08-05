#!/usr/bin/env python3
"""
harmonize_names.py -- give homologous chromosomes consistent names across a set
of same-species assemblies, using one in-batch assembly as the naming reference.

RENAME-AND-REORIENT ONLY. This never re-scaffolds: scaffold count, boundaries and
internal gaps are invariant. The reference's chromosome-scale scaffolds define the
chromosome slots (chrN, N = descending-size rank); every other assembly inherits
those names by minimap2 homology.

Name grammar (always suffixed, so the output is uniform for later graph partitioning):

    chr{K}_{part}                  one scaffold covering chromosome K, piece <part>
    chr{A}_{i}+chr{B}_{j}[+...]    one scaffold spanning >1 chromosome (fusion / chimera)
    unplaced_{n}                   no confident chromosome placement

  * Split case (2 scaffolds cover chr K): both get primary K, ordered by their
    reference start -> chrK_1, chrK_2. Each is reverse-complemented to match K.
  * Merge / fusion (1 scaffold covers chr A and chr B): composite name carrying a
    per-chromosome part token for each chromosome it touches; NOT reoriented (its two
    segments can disagree on strand), so left native. Flagged, with batch concordance.

Placement is decided by alignment, not by size: a small scaffold that aligns
confidently to chr1 is chr1; a large scaffold that aligns nowhere is unplaced.

Inputs via a manifest TSV (header: id  role  fai  paf):
    - one 'ref' row  (paf column = '-')
    - one 'query' row per other assembly (paf = minimap2 <reference> <query> PAF; .gz ok)
Each 'fai' is a samtools .fai (name<TAB>length<TAB>...); only the first two cols are read.

Outputs (into --outdir):
    {id}.harmonized_name_map.tsv   per assembly (incl. the reference)
    {species}.harmonization_report.tsv   one row per (assembly, scaffold) + flags
"""

import argparse
import csv
import gzip
import os
import sys
from collections import defaultdict
from itertools import combinations


def openf(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fai(path):
    """Return [(name, length), ...] in file order."""
    out = []
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            out.append((f[0], int(f[1])))
    return out


def parse_paf(path, chrom_of):
    """Accumulate per (query_scaffold, chromosome) alignment evidence.

    chrom_of maps a *reference scaffold name* -> chromosome number, and only
    contains chromosome-scale reference scaffolds. Alignments to sub-threshold
    reference sequences are ignored (they confer no chromosome name).

    Returns:
        ev:   qname -> chr_num -> {bp, tstart, tend, plus, minus}
    """
    ev = defaultdict(lambda: defaultdict(
        lambda: {"bp": 0, "tstart": None, "tend": 0, "plus": 0, "minus": 0}))
    if not path or path == "-" or not os.path.exists(path):
        return ev
    with openf(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 12:
                continue
            q = c[0]
            strand = c[4]
            t = c[5]
            ts = int(c[7])
            te = int(c[8])
            blk = int(c[10])          # alignment block length -> aligned bp
            if t not in chrom_of:
                continue
            k = chrom_of[t]
            e = ev[q][k]
            e["bp"] += blk
            e["tstart"] = ts if e["tstart"] is None else min(e["tstart"], ts)
            e["tend"] = max(e["tend"], te)
            if strand == "-":
                e["minus"] += blk
            else:
                e["plus"] += blk
    return ev


def classify(chr_ev, scaf_len, min_frac, sec_frac):
    """Decide placement of one query scaffold from its per-chromosome evidence."""
    total_bp = sum(e["bp"] for e in chr_ev.values())
    if scaf_len == 0 or total_bp == 0 or (total_bp / scaf_len) < min_frac:
        return {"class": "unplaced",
                "aligned_frac": (total_bp / scaf_len) if scaf_len else 0.0}
    by_bp = sorted(chr_ev.items(), key=lambda kv: kv[1]["bp"], reverse=True)
    primary = by_bp[0][0]
    partners = [k for k, e in by_bp[1:] if e["bp"] >= sec_frac * total_bp]
    members = sorted([primary] + partners)
    orient = {k: ("-" if chr_ev[k]["minus"] > chr_ev[k]["plus"] else "+")
              for k in members}
    return {
        "class": "composite" if partners else "chromosome",
        "primary": primary,
        "members": members,
        "tstart": {k: (chr_ev[k]["tstart"] or 0) for k in members},
        "tend": {k: chr_ev[k]["tend"] for k in members},
        "orient": orient,
        "aligned_frac": total_bp / scaf_len,
    }


def select_chromosome_set(ref_fai, min_scaffold_bp, method, dropoff_ratio, dropoff_min_frac):
    """Pick the reference chromosome set (descending by length).

    method='threshold' -> every reference scaffold >= min_scaffold_bp.
    method='dropoff'   -> cut at the sharpest adjacent length ratio among boundaries where
                          the last chromosome is >= min_scaffold_bp (floor) and the
                          cumulative genome fraction has reached dropoff_min_frac (so we
                          don't cut after a single dominant chromosome). The floor also
                          excludes deep-tail cuts. If no boundary reaches dropoff_ratio,
                          fall back to the plain threshold and flag it.

    Returns (chrom_list, meta) where meta records method / n / cut ratio / genome fraction
    / flags so the choice is auditable in the report.
    """
    ordered = sorted(ref_fai, key=lambda x: (-x[1], x[0]))
    n = len(ordered)
    total = sum(L for _, L in ordered) or 1
    meta = {"method": method, "n_chrom": 0, "cut_ratio": None,
            "genome_fraction": 0.0, "flags": []}

    def threshold_set():
        return [(nm, L) for nm, L in ordered if L >= min_scaffold_bp]

    if method == "threshold" or n <= 1:
        cs = threshold_set()
        meta["method"] = "threshold"
    else:
        cum, cf = 0.0, []
        for _, L in ordered:
            cum += L
            cf.append(cum / total)
        best_k, best_ratio = None, -1.0
        for k in range(1, n):                       # ordered[k] must exist for the ratio
            last_len, nxt_len = ordered[k-1][1], ordered[k][1]
            if last_len < min_scaffold_bp or cf[k-1] < dropoff_min_frac:
                continue
            ratio = last_len / nxt_len if nxt_len > 0 else float("inf")
            if ratio > best_ratio:
                best_ratio, best_k = ratio, k
        if best_k is not None and best_ratio >= dropoff_ratio:
            cs = ordered[:best_k]
            meta.update(method="dropoff", cut_ratio=best_ratio)
        else:
            cs = threshold_set()
            meta.update(method="threshold_fallback",
                        cut_ratio=(best_ratio if best_k is not None else None))
            meta["flags"].append("no_sharp_dropoff")

    meta["n_chrom"] = len(cs)
    meta["genome_fraction"] = sum(L for _, L in cs) / total
    if cs and meta["genome_fraction"] < 0.7:
        meta["flags"].append("low_genome_fraction")
    return cs, meta


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True, help="TSV: id  role  fai  paf")
    ap.add_argument("--reference-id", required=True)
    ap.add_argument("--species", default="species", help="label for the report file")
    ap.add_argument("--min-scaffold-bp", type=int, default=1_000_000,
                    help="floor for a chromosome-scale scaffold; also the threshold used "
                         "when the drop-off method finds no clear boundary")
    ap.add_argument("--chromosome-set-method", choices=("dropoff", "threshold"),
                    default="dropoff",
                    help="pick the reference chromosome set by data-driven length "
                         "drop-off (default) or a plain --min-scaffold-bp threshold")
    ap.add_argument("--dropoff-ratio", type=float, default=2.0,
                    help="min adjacent length ratio to accept a drop-off boundary")
    ap.add_argument("--dropoff-min-frac", type=float, default=0.5,
                    help="cumulative genome fraction reached before a drop-off cut is "
                         "allowed (guards against cutting after a dominant chromosome)")
    ap.add_argument("--min-aligned-frac", type=float, default=0.5,
                    help="min (aligned-to-chromosomes bp / scaffold length) to place")
    ap.add_argument("--secondary-frac", type=float, default=0.2,
                    help="a 2nd chromosome contributing >= this share -> composite")
    ap.add_argument("--overlap-tol", type=int, default=0,
                    help="bp of reference overlap between two pieces before flagging")
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()

    with open(a.manifest) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    by_id = {r["id"]: r for r in rows}
    if a.reference_id not in by_id:
        sys.exit(f"reference id {a.reference_id!r} not in manifest")

    # ---- reference chromosome set (chrN = descending-size rank) ----
    ref_fai = read_fai(by_id[a.reference_id]["fai"])
    chrom, chrom_meta = select_chromosome_set(
        ref_fai, a.min_scaffold_bp, a.chromosome_set_method,
        a.dropoff_ratio, a.dropoff_min_frac)
    chrom_of = {n: i for i, (n, _) in enumerate(chrom, start=1)}  # ref scaffold -> chr num
    n_chrom = len(chrom)
    if n_chrom == 0:
        sys.stderr.write(
            f"[harmonize {a.species}] WARNING: reference {a.reference_id} yields no "
            f"chromosome-scale scaffolds; nothing will be placed. Lower --min-scaffold-bp "
            f"or use a chromosome-scale reference.\n")

    # ---- resolve every assembly into placements ----
    place = {}          # id -> {scaffold -> info}
    present_chr = {}    # id -> set(chr nums present as a single primary)
    fused_pairs = {}    # id -> set(frozenset({a,b}))

    for rid, row in by_id.items():
        fai = read_fai(row["fai"])
        if rid == a.reference_id:
            p = {}
            for n, L in fai:
                if n in chrom_of:
                    k = chrom_of[n]
                    p[n] = {"class": "chromosome", "primary": k, "members": [k],
                            "tstart": {k: 0}, "tend": {k: L}, "orient": {k: "+"},
                            "aligned_frac": 1.0, "length": L}
                else:
                    p[n] = {"class": "unplaced", "length": L, "aligned_frac": 0.0}
        else:
            ev = parse_paf(row["paf"], chrom_of)
            p = {}
            for n, L in fai:
                info = classify(ev.get(n, {}), L, a.min_aligned_frac, a.secondary_frac)
                info["length"] = L
                p[n] = info
        place[rid] = p
        pres, fused = set(), set()
        for info in p.values():
            if info["class"] == "chromosome":
                pres.add(info["primary"])
            elif info["class"] == "composite":
                for x, y in combinations(info["members"], 2):
                    fused.add(frozenset((x, y)))
        present_chr[rid] = pres
        fused_pairs[rid] = fused

    def concordance(x, y):
        pair = frozenset((x, y))
        n_f = sum(1 for s in fused_pairs.values() if pair in s)
        n_s = sum(1 for r in by_id if x in present_chr[r] and y in present_chr[r])
        return n_f, n_s

    # ---- assign per-chromosome part indices, names, order ----
    results = {}
    for rid in by_id:
        p = place[rid]
        claimants = defaultdict(list)          # chr -> [(tstart, scaffold)]
        for n, info in p.items():
            if info["class"] in ("chromosome", "composite"):
                for k in info["members"]:
                    claimants[k].append((info["tstart"][k], n))
        part = defaultdict(dict)               # scaffold -> {chr: index}
        overlap = defaultdict(set)             # scaffold -> set(chr)
        for k, lst in claimants.items():
            lst.sort(key=lambda x: (x[0], x[1]))
            prev_end, prev = None, None
            for idx, (ts, n) in enumerate(lst, start=1):
                part[n][k] = idx
                if prev_end is not None and ts < prev_end - a.overlap_tol:
                    overlap[n].add(k)
                    overlap[prev].add(k)
                prev_end, prev = p[n]["tend"][k], n

        unplaced = sorted([(n, info["length"]) for n, info in p.items()
                           if info["class"] == "unplaced"],
                          key=lambda x: (-x[1], x[0]))
        unplaced_idx = {n: i for i, (n, _) in enumerate(unplaced, start=1)}

        out = []
        for n, info in p.items():
            flags = []
            if info["class"] == "chromosome":
                k = info["primary"]
                new = f"chr{k}_{part[n][k]}"
                orient = "rev" if info["orient"][k] == "-" else "fwd"
                span = f"chr{k}:{info['tstart'][k]}-{info['tend'][k]}"
                key = (0, k, part[n][k])
                if k in overlap[n]:
                    flags.append("overlap")
            elif info["class"] == "composite":
                toks = [f"chr{k}_{part[n][k]}" for k in info["members"]]
                new = "+".join(toks)
                orient = "fwd"                 # never reorient a composite
                span = ";".join(f"chr{k}:{info['tstart'][k]}-{info['tend'][k]}"
                                for k in info["members"])
                prim = info["primary"]
                key = (0, prim, part[n][prim])
                flags.append("fusion")
                for x, y in combinations(info["members"], 2):
                    n_f, n_s = concordance(x, y)
                    tag = "concordant" if (n_f > 1 or n_s == 0) else "chimera_suspect"
                    flags.append(f"{tag}(chr{min(x, y)}+chr{max(x, y)}:{n_f}f/{n_s}s)")
                if any(k in overlap[n] for k in info["members"]):
                    flags.append("overlap")
            else:
                u = unplaced_idx[n]
                new = f"unplaced_{u}"
                orient, span = "fwd", "-"
                key = (1, u, 0)
                if info.get("aligned_frac", 0) > 0:
                    flags.append(f"low_cov({info['aligned_frac']:.2f})")
            out.append({"old": n, "new": new, "orient": orient, "length": info["length"],
                        "class": info["class"], "span": span,
                        "flags": ";".join(flags) if flags else "-", "_key": key})
        out.sort(key=lambda d: d["_key"])
        for i, d in enumerate(out, start=1):
            d["order"] = i
        results[rid] = out

    # ---- write outputs ----
    os.makedirs(a.outdir, exist_ok=True)
    for rid, out in results.items():
        with open(os.path.join(a.outdir, f"{rid}.harmonized_name_map.tsv"), "w") as fh:
            fh.write("old_name\tnew_name\torient\torder\tlength\tclass\tref_span\tflags\n")
            for d in out:
                fh.write("\t".join(str(x) for x in (
                    d["old"], d["new"], d["orient"], d["order"],
                    d["length"], d["class"], d["span"], d["flags"])) + "\n")

    rep = os.path.join(a.outdir, f"{a.species}.harmonization_report.tsv")
    with open(rep, "w") as fh:
        fh.write(f"# reference_id\t{a.reference_id}\n")
        fh.write(f"# n_reference_chromosomes\t{n_chrom}\n")
        fh.write(f"# chromosome_set_method\t{chrom_meta['method']}\n")
        fh.write(f"# chromosome_set_cut_ratio\t"
                 f"{('%.2f' % chrom_meta['cut_ratio']) if chrom_meta['cut_ratio'] is not None else 'NA'}\n")
        fh.write(f"# chromosome_set_genome_fraction\t{chrom_meta['genome_fraction']:.3f}\n")
        fh.write(f"# chromosome_set_flags\t{';'.join(chrom_meta['flags']) if chrom_meta['flags'] else '-'}\n")
        fh.write(f"# params\tmin_scaffold_bp={a.min_scaffold_bp}\t"
                 f"dropoff_ratio={a.dropoff_ratio}\tdropoff_min_frac={a.dropoff_min_frac}\t"
                 f"min_aligned_frac={a.min_aligned_frac}\t"
                 f"secondary_frac={a.secondary_frac}\n")
        fh.write("assembly\told_name\tnew_name\tclass\tlength\torient\tref_span\tflags\n")
        for rid, out in results.items():
            for d in out:
                fh.write("\t".join(str(x) for x in (
                    rid, d["old"], d["new"], d["class"], d["length"],
                    d["orient"], d["span"], d["flags"])) + "\n")

    n_flagged = sum(1 for out in results.values() for d in out if d["flags"] != "-")
    sys.stderr.write(
        f"[harmonize {a.species}] reference={a.reference_id} "
        f"chromosomes={n_chrom} ({chrom_meta['method']}, "
        f"{chrom_meta['genome_fraction']:.2f} of genome) "
        f"assemblies={len(results)} flagged_scaffolds={n_flagged}\n")


if __name__ == "__main__":
    main()
