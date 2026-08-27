#!/usr/bin/env python3
# ======================================================================================
# rearrange_from_untangle.py
#
# Turns per-chromosome `odgi untangle` output into rearrangement calls, and -- more
# importantly -- into a CANDIDATE TABLE that surfaces a segregating inversion without
# anyone having to write awk.
#
# WHY THE CANDIDATE TABLE IS THE POINT
# ------------------------------------
# A candidate segregating inversion on chr10 was found during development, and finding it
# took several rounds of ad-hoc awk over the raw untangle output: summing inverted bp per
# query, spotting that only 2 of 10 haplotypes carried anything, noticing both carriers'
# sister haplotypes were clean (so both individuals are heterozygous, AF 0.2), then
# separately checking that the carriers were not simply mis-oriented scaffolds. None of
# that was visible in any published output. This script makes it a row.
#
# THE ARTIFACT DISCRIMINATOR IS HARMONIZATION'S FLAGS, NOT A PROJECTION STATISTIC
# ------------------------------------------------------------------------------
# A scaffold whose two halves belong on different chromosomes has NO correct single
# orientation, so one half projects as 100% inverted. Two such cases were found:
#
#   Sde-CPla_115#1#chr4_7+chr10_10    10 bp fwd, 715,908 bp inv  (100.0%)
#   Sde-CPla_115#1#chr10_17+chr11_5    0 bp fwd, 1,415,923 bp inv (100.0%)
#
# Both were already `fwd` in harmonization and did not flip under the merged-footprint
# patch, so this is NOT an orientation-call error -- it is a chimeric join. Harmonization
# already says so: `unsupported(ref4+ref10:0f/8s)` means that junction was seen fused in 0
# voters and split in 8. Zero fused voters is a far sharper tell than "100% inverted", so
# the join on harmonization's flags is the primary discriminator and per-query %inverted is
# secondary. Both flag spellings matter: `unsupported(...)` for no voter support and
# `chimera_suspect(...)` for a fused minority (e.g. `1f/7s`).
#
# WHAT THE %INVERTED COLUMN IS STILL FOR
# --------------------------------------
# It separates whole-scaffold orientation problems from internal inversions. The chr10
# carriers ran 3.8%, 15.3% and 45.9% inverted -- a failed reverse-complement would read
# ~100% -- which is how the candidate survived scrutiny.
#
# TWO FOOTPRINT MEASURES, DELIBERATELY
# ------------------------------------
# union    union of inverted reference intervals. Tightens as -e improves (13.95 -> 12.85 Mb
#          on chr10 across a 100x resolution change) because coarse segments swallow
#          non-inverted interstitial sequence.
# span     same intervals merged across gaps below --run-merge-gap. Held at 14.45-14.49 Mb
#          across that same sweep.
# Their ratio (fill) is a quality signal: 96.5% at -e 1 Mb, 88.7% at 10 kb. A stable span
# with a tightening union is what a real polymorphism looks like under increasing
# resolution; a segmentation artifact does not behave that way.
#
# Dependencies: python3 stdlib.
# ======================================================================================

import argparse
import glob
import gzip
import os
import re
import sys
from collections import defaultdict

FLAG_RE = re.compile(r"(unsupported|chimera_suspect)\(([^)]*)\)")


def vopen(p):
    with open(p, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return (gzip.open if gz else open)(p, "rt", encoding="utf-8", errors="replace")


def merged(iv, gap=0):
    """Union length and block count of [start,end) intervals, merging gaps <= gap."""
    if not iv:
        return 0, 0
    iv = sorted(iv)
    total, blocks = 0, 1
    s, e = iv[0]
    for a, b in iv[1:]:
        if a > e + gap:
            total += e - s
            blocks += 1
            s, e = a, b
        elif b > e:
            e = b
    return total + (e - s), blocks


def pansn_to_report(q):
    """Sde-CBau_104#1#chr10_1#0[...] -> ('Sde-CBau_104_hap1', 'chr10_1').

    Harmonization reports assembly ids as <sample>_hap<N>; the graph uses PanSN
    <sample>#<hap>#<contig>. The reference carries #0 and no _hapN suffix.
    """
    p = q.split("#")
    if len(p) < 3:
        return q, q
    sample, hap, contig = p[0], p[1], p[2]
    contig = contig.split("[")[0]
    asm = sample if hap == "0" else "%s_hap%s" % (sample, hap)
    return asm, contig


def load_harmonization(path):
    """assembly|new_name -> (class, orient, flags). Tolerant of column drift: locates the
    header row rather than assuming fixed positions."""
    if not path or not os.path.isfile(path):
        return {}
    out = {}
    hdr = None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if hdr is None:
                if f and f[0] == "assembly":
                    hdr = {name: i for i, name in enumerate(f)}
                continue
            def g(name, dflt=""):
                i = hdr.get(name, -1)
                return f[i] if 0 <= i < len(f) else dflt
            key = "%s|%s" % (g("assembly"), g("new_name"))
            out[key] = (g("class"), g("orient"), g("flags"))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--untangle", nargs="+", required=True,
                   help="per-chromosome untangle TSV(.gz) files")
    p.add_argument("--outdir", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--flavor", default="full")
    p.add_argument("--harmonization", default=None,
                   help="<taxid>.harmonization_report.tsv, for the artifact join")
    p.add_argument("--min-jaccard", type=float, default=0.1)
    p.add_argument("--min-seg-bp", type=int, default=1000)
    p.add_argument("--run-merge-gap", type=int, default=200000)
    p.add_argument("--artifact-frac", type=float, default=0.95,
                   help="per-query inverted fraction at or above which the query is flagged "
                        "as a probable whole-scaffold orientation artifact (default 0.95)")
    a = p.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    stem = "%s.%s" % (a.label, a.flavor)

    def op(s):
        return os.path.join(a.outdir, stem + s)

    harm = load_harmonization(a.harmonization)

    files = []
    for pat in a.untangle:
        files.extend(sorted(glob.glob(pat)) if any(c in pat for c in "*?[") else [pat])

    q_fwd = defaultdict(int)
    q_inv = defaultdict(int)
    # unfiltered counterparts: pct_inv computed on FILTERED rows can be inflated when -j
    # removes a query's forward segments, which would make an ordinary scaffold look like an
    # orientation artifact. Both are reported so the filter cannot silently shape the call.
    q_fwd_all = defaultdict(int)
    q_inv_all = defaultdict(int)
    q_dup = defaultdict(int)
    inv_iv = defaultdict(list)          # (refname, query) -> reference intervals
    dup_rows = []
    n_rows = n_kept = n_lowscore = n_small = 0
    refnames = defaultdict(set)

    fbed = open(op(".inversions.bed"), "w")
    fbed.write('track name="%s_inversions" description="untangle minus-strand segments"\n'
               % stem)

    for path in files:
        with vopen(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 10:
                    continue
                n_rows += 1
                q, qs, qe = f[0], int(f[1]), int(f[2])
                rn, rs, re_ = f[3], int(f[4]), int(f[5])
                try:
                    score = float(f[6])
                except ValueError:
                    score = 0.0
                strand, selfcov = f[7], f[8]
                qbp = qe - qs
                if strand == "-":
                    q_inv_all[q] += qbp
                else:
                    q_fwd_all[q] += qbp
                if score < a.min_jaccard:
                    n_lowscore += 1
                    continue
                if qbp < a.min_seg_bp:
                    n_small += 1
                    continue
                n_kept += 1
                refnames[q].add(rn)
                if strand == "-":
                    q_inv[q] += qbp
                    inv_iv[(rn, q)].append((rs, re_))
                    fbed.write("%s\t%d\t%d\t%s|inv\t%d\t-\n"
                               % (rn, rs, re_, q, min(1000, int(score * 1000))))
                else:
                    q_fwd[q] += qbp
                try:
                    if float(selfcov) > 1.0:
                        q_dup[q] += qbp
                        dup_rows.append((rn, rs, re_, q, float(selfcov), qbp))
                except ValueError:
                    pass
    fbed.close()

    # ---- per-query orientation summary: the artifact discriminator -------------------
    with open(op(".query_orientation.tsv"), "w") as out:
        out.write("# pct_inv near 100 is a whole-scaffold ORIENTATION problem, not an\n")
        out.write("#   internal inversion. A scaffold spanning two reference chromosomes has\n")
        out.write("#   no correct single orientation, so one half projects fully inverted.\n")
        out.write("# harm_flags is the PRIMARY discriminator: unsupported(...:0f/Ns) means no\n")
        out.write("#   voter supports that junction, chimera_suspect(...) a fused minority.\n")
        out.write("# pct_inv uses FILTERED rows; pct_inv_unfiltered uses every row. A large\n")
        out.write("#   gap between them means the -j/-size filters, not orientation, are\n")
        out.write("#   driving the number. ORIENTATION_SUSPECT requires BOTH to be high.\n")
        out.write("query\tassembly\tcontig\tfwd_bp\tinv_bp\tpct_inv\tpct_inv_unfiltered\t"
                  "dup_bp\tn_ref_targets\tharm_class\tharm_orient\tharm_flags\t"
                  "artifact_flag\n")
        for q in sorted(set(q_fwd) | set(q_inv) | set(q_fwd_all) | set(q_inv_all)):
            fwd, inv = q_fwd.get(q, 0), q_inv.get(q, 0)
            tot = fwd + inv
            pct = 100.0 * inv / tot if tot else 0.0
            tot_all = q_fwd_all.get(q, 0) + q_inv_all.get(q, 0)
            pct_all = 100.0 * q_inv_all.get(q, 0) / tot_all if tot_all else 0.0
            asm, contig = pansn_to_report(q)
            hc, ho, hf = harm.get("%s|%s" % (asm, contig), ("", "", ""))
            flags = []
            # both measures must agree, so a filtering effect alone cannot raise the flag
            if pct >= 100.0 * a.artifact_frac and pct_all >= 100.0 * a.artifact_frac:
                flags.append("ORIENTATION_SUSPECT")
            if FLAG_RE.search(hf or ""):
                flags.append("CHIMERIC_JUNCTION")
            if "composite" in (hc or ""):
                flags.append("COMPOSITE")
            out.write("%s\t%s\t%s\t%d\t%d\t%.4f\t%.4f\t%d\t%d\t%s\t%s\t%s\t%s\n"
                      % (q, asm, contig, fwd, inv, pct, pct_all, q_dup.get(q, 0),
                         len(refnames.get(q, ())), hc, ho, hf, ";".join(flags) or "."))

    # ---- candidate loci: the table that should have surfaced chr10 -------------------
    # Group inverted intervals into loci per reference contig, then report which
    # haplotypes carry each locus. Two carriers out of ten with clean sister haplotypes is
    # a segregating polymorphism; ten carriers is reference divergence; one is an assembly
    # artifact until proven otherwise.
    loci = defaultdict(list)
    for (rn, q), ivs in inv_iv.items():
        for s, e in ivs:
            loci[rn].append((s, e, q))
    with open(op(".rearrangement_candidates.tsv"), "w") as out:
        out.write("# inverted segments merged into loci per reference contig (gap <= %d).\n"
                  % a.run_merge_gap)
        out.write("# SORTED SO A LOW-CARRIER-COUNT LOCUS SURFACES FIRST -- a locus carried by\n")
        out.write("#   2 of N haplotypes whose sister haplotypes are clean is a candidate\n")
        out.write("#   segregating inversion; one carried by all is reference divergence.\n")
        out.write("# union_bp tightens with finer -e; span_bp is stable. fill = union/span is\n")
        out.write("#   a quality signal (0.887 at -e 10kb vs 0.965 at 1Mb on chr10).\n")
        out.write("# any_artifact_flag set means at least one carrier is chimeric, composite,\n")
        out.write("#   or ~fully inverted -- treat the locus as suspect, not as biology.\n")
        out.write("ref_contig\tlocus_start\tlocus_end\tspan_bp\tunion_bp\tfill\t"
                  "n_carriers\tcarriers\tany_artifact_flag\n")
        rows = []
        for rn, items in loci.items():
            items.sort()
            cs, ce = items[0][0], items[0][1]
            members = [items[0]]
            groups = []
            for s, e, q in items[1:]:
                if s > ce + a.run_merge_gap:
                    groups.append((cs, ce, members))
                    cs, ce, members = s, e, [(s, e, q)]
                else:
                    ce = max(ce, e)
                    members.append((s, e, q))
            groups.append((cs, ce, members))
            for gs, ge, mem in groups:
                carriers = sorted({q for _, _, q in mem})
                union, _ = merged([(s, e) for s, e, _ in mem])
                span, _ = merged([(s, e) for s, e, _ in mem], a.run_merge_gap)
                art = []
                for q in carriers:
                    asm, contig = pansn_to_report(q)
                    hc, ho, hf = harm.get("%s|%s" % (asm, contig), ("", "", ""))
                    fwd, inv = q_fwd.get(q, 0), q_inv.get(q, 0)
                    pct = 100.0 * inv / (fwd + inv) if (fwd + inv) else 0.0
                    if FLAG_RE.search(hf or "") or "composite" in (hc or "") \
                            or pct >= 100.0 * a.artifact_frac:
                        art.append(q)
                # sort key: unflagged before flagged, then fewest carriers, then largest.
                # A clean 2-of-N locus is the interesting case; a flagged 1-carrier locus is
                # an artifact and must not occupy the top of the table.
                rows.append((1 if art else 0, len(carriers), -span,
                             rn, gs, ge, span, union, carriers, art))
        rows.sort()
        for _flagged, nc, negspan, rn, gs, ge, span, union, carriers, art in rows:
            out.write("%s\t%d\t%d\t%d\t%d\t%s\t%d\t%s\t%s\n"
                      % (rn, gs, ge, span, union,
                         ("%.4f" % (union / span)) if span else "NA",
                         nc, ",".join(carriers), ",".join(art) or "."))

    with open(op(".duplications.tsv"), "w") as out:
        out.write("# self.cov > 1. THE ONLY duplication signal available: AT traversals show\n")
        out.write("# just 27 node re-visits across 3,268,312 SV alleles, so tandem\n")
        out.write("# duplications are structurally unrepresentable in this graph's bubbles.\n")
        out.write("ref_contig\tref_start\tref_end\tquery\tself_cov\tquery_bp\n")
        for r in sorted(dup_rows, key=lambda x: -x[5]):
            out.write("%s\t%d\t%d\t%s\t%.4f\t%d\n" % r)

    with open(op(".untangle_audit.tsv"), "w") as out:
        out.write("metric\tvalue\n")
        out.write("flavor\t%s\nfiles\t%d\n" % (a.flavor, len(files)))
        out.write("rows_read\t%d\nrows_kept\t%d\n" % (n_rows, n_kept))
        out.write("dropped_score_below_%.3f\t%d\n" % (a.min_jaccard, n_lowscore))
        out.write("dropped_shorter_than_%d\t%d\n" % (a.min_seg_bp, n_small))
        out.write("queries\t%d\n" % len(set(q_fwd) | set(q_inv)))
        out.write("inverted_bp_total\t%d\n" % sum(q_inv.values()))
        out.write("forward_bp_total\t%d\n" % sum(q_fwd.values()))
        out.write("inverted_bp_total_unfiltered\t%d\n" % sum(q_inv_all.values()))
        out.write("forward_bp_total_unfiltered\t%d\n" % sum(q_fwd_all.values()))
        out.write("dup_bp_total\t%d\n" % sum(q_dup.values()))
        out.write("harmonization_rows_loaded\t%d\n" % len(harm))
        out.write("# filters are REPORTED, not silent: a large dropped count means the\n")
        out.write("# thresholds, not the biology, are shaping the output.\n")

    sys.stderr.write("[rearrange] %s: %d rows -> %d kept (%d low score, %d short)\n"
                     % (stem, n_rows, n_kept, n_lowscore, n_small))
    sys.stderr.write("[rearrange] inverted bp %d, forward bp %d, dup bp %d\n"
                     % (sum(q_inv.values()), sum(q_fwd.values()), sum(q_dup.values())))
    sys.stderr.write("[rearrange] harmonization rows loaded: %d%s\n"
                     % (len(harm), "" if harm else "  (NO ARTIFACT JOIN -- flags unavailable)"))


if __name__ == "__main__":
    main()
