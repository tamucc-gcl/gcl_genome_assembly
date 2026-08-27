#!/usr/bin/env python3
# ======================================================================================
# classify_variants.py
#
# Replaces the length-based awk classifier in pangenome_variants.nf with a topological one
# derived from the graph's own allele traversals, plus an allele-frequency layer.
#
# Promoted from audit_variant_layers.py v2 after validation on Spratelloides:
#   - the length classifier reproduces the published variant_summary.tsv exactly once two
#     awk branches added since that artifact are accounted for
#   - AT is present on 41,920,948/41,920,948 raw records, 0 unparseable
#
# THREE THINGS THIS ENFORCES THAT THE PROTOTYPE LEARNED THE HARD WAY
#
#   1. AT IS ONLY VALID PRE-DECOMPOSITION. vcfwave and `bcftools norm -m` rewrite REF/ALT
#      while AT is inherited from the parent record, so alt index i stops corresponding to
#      traversal i+1. The guard below refuses topology on any VCF whose header shows
#      decomposition rather than trusting the wiring to feed the right file. Length class
#      and AF are still computed, so the fine view runs through the same code path.
#
#   2. DUP MUST BE RESTRICTED TO NODES THE REFERENCE ALSO VISITS. Testing
#      `alt_count[n] > ref_count.get(n, 0)` fires for any novel node, making DUP a synonym
#      for NOVEL_INS -- on Spratelloides both returned an identical 493,623. DUP_NOVEL
#      (a novel node visited more than once) is reported separately.
#
#   3. SUMMING PER-ALLELE bp IS NOT A GENOME-WIDE TOTAL. 41.9M records carry 55.3M alt
#      alleles, so a 7-alt record's reference span is counted seven times; SUBST alleles
#      >=100 kb alone sum to 1.10 Gb of longer-allele bp in a ~1 Gb genome. Per-allele bp is
#      correct for a SIZE SPECTRUM (x axis) and meaningless as a total. This emits the
#      MERGED, NON-OVERLAPPING REFERENCE FOOTPRINT per class as the defensible total, and
#      labels the per-allele column so it cannot be mistaken for one.
#
# LABELS -- non-exclusive, because a locus can be an inverted duplication
#   INV_PATH_EXPLICIT  contiguous run of >=2 shared nodes, flipped orientation AND reversed
#                      order. Named for what it measures: per Romain et al. (bioRxiv
#                      2025.03.14.643331) only the path-explicit motif is visible to
#                      topology. Alignment-rescued inversions are disjoint paths and land in
#                      SUBST; rescue_inversions.py finds those. On Spratelloides this label
#                      is 30 alleles while odgi untangle finds ~21.5 Mb inverted on chr10
#                      alone, so THIS NUMBER IS A FLOOR, NOT AN INVERSION COUNT.
#   DUP / DUP_NOVEL    copy-number gain of reference content / tandem gain of novel content
#   NOVEL_INS / DEL    node ids present only in alt / only in ref
#   REORDER            shared nodes non-monotone in ref order, orientations unchanged
#
# primary_class is the EXCLUSIVE partition. Anything that sums uses it.
#
# Dependencies: python3 stdlib.
#
# USAGE
#   classify_variants.py --vcf parents.vcf.gz --outdir . --label <taxid> \
#       --flavor clip --tier parent [--chrom chr1_1]
# ======================================================================================

import argparse
import gzip
import os
import sys
from collections import Counter, defaultdict

try:
    import numpy as _np
except ImportError:
    _np = None

# --------------------------------------------------------------------------------------
_DECOMP_MARKERS = ("vcfwave", "ID=ORIGIN", "bcftools_normCommand",
                   "bcftools_normVersion", "vcfallelicprimitives")

SV_CLASSES = ("SV_INS", "SV_DEL", "SV_COMPLEX", "SV_BLOCKSUB")
TOPO_LABELS = ("INV_PATH_EXPLICIT", "INV_WEAK", "DUP", "DUP_NOVEL",
               "NOVEL_INS", "DEL", "REORDER")


def vopen(path):
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return (gzip.open if gz else open)(path, "rt", encoding="utf-8", errors="replace")


def parse_info(s):
    d = {}
    for kv in s.split(";"):
        if not kv:
            continue
        k, sep, v = kv.partition("=")
        d[k] = v if sep else True
    return d


def parse_traversal(t):
    out = []
    i, n = 0, len(t)
    while i < n:
        c = t[i]
        if c == ">":
            o = 1
        elif c == "<":
            o = -1
        else:
            return None
        j = i + 1
        while j < n and t[j] not in "><":
            j += 1
        tok = t[i + 1:j]
        if not tok.isdigit():
            return None
        out.append((int(tok), o))
        i = j
    return out or None


# --------------------------------------------------------------------------------------
def classify_traversal(ref_t, alt_t):
    labels = set()
    ref_ids = [n for n, _ in ref_t]
    alt_ids = [n for n, _ in alt_t]
    rc, ac = Counter(ref_ids), Counter(alt_ids)
    shared = set(rc) & set(ac)

    if any(n not in rc for n in alt_ids):
        labels.add("NOVEL_INS")
    if any(n not in ac for n in ref_ids):
        labels.add("DEL")

    # restricted to shared nodes -- novel sequence is NOT a duplication (prototype bug 2)
    if any(ac[n] > rc[n] for n in shared):
        labels.add("DUP")
    if any(ac[n] > 1 for n in ac if n not in rc):
        labels.add("DUP_NOVEL")

    ro, ramb = {}, set()
    for n, o in ref_t:
        if n in ro and ro[n] != o:
            ramb.add(n)
        ro.setdefault(n, o)
    ridx = {}
    for k, n in enumerate(ref_ids):
        ridx.setdefault(n, k)

    seq = [(ridx[n], o != ro[n]) for n, o in alt_t if n in shared and n not in ramb]
    nflip = sum(1 for _, f in seq if f)
    run = best = 0
    prev = None
    for idx, flip in seq:
        if flip and prev is not None and idx < prev:
            run = run + 1 if run else 2
        else:
            run = 1 if flip else 0
        best = max(best, run)
        prev = idx
    if best >= 2:
        labels.add("INV_PATH_EXPLICIT")
    elif nflip == 1:
        labels.add("INV_WEAK")

    fwd = [idx for idx, f in seq if not f]
    if len(fwd) >= 3 and any(b < a for a, b in zip(fwd, fwd[1:])) \
            and "INV_PATH_EXPLICIT" not in labels:
        labels.add("REORDER")
    return labels


def length_class(rl, al, minsv):
    d, mn = abs(rl - al), min(rl, al)
    if rl == 1 and al == 1:
        return "SNP"
    if d < minsv and mn < minsv:
        return "INDEL"
    if d < minsv:
        return "SV_COMPLEX"
    if mn >= minsv:
        return "SV_BLOCKSUB"
    return "SV_INS" if al > rl else "SV_DEL"


def primary_class(labels, lc):
    """Exclusive partition. SUBST is deliberately NOT subdivided here.

    rescue_inversions.py measures frac_fwd/frac_rev per SUBST allele and the join step
    splits SUBST into SUBST_HOMOLOGOUS / SUBST_PARTIAL / SUBST_UNRELATED from that
    measurement. On Spratelloides 4,015 of 4,329 SUBST alleles >=100 kb align to their
    reference allele at >=0.5 with no inversion, so SUBST is largely a statement about how
    the GRAPH chose to represent homologous sequence, not about the sequences themselves.
    """
    if lc in ("SNP", "INDEL"):
        return lc
    dup = ("DUP" in labels) or ("DUP_NOVEL" in labels)
    if "INV_PATH_EXPLICIT" in labels:
        return "INV_DUP" if dup else "INV_PATH_EXPLICIT"
    if dup:
        return "DUP"
    if "NOVEL_INS" in labels and "DEL" in labels:
        return "SUBST"
    if "NOVEL_INS" in labels:
        return "INS"
    if "DEL" in labels:
        return "DEL"
    if "REORDER" in labels:
        return "REORDER"
    return "UNCLASSIFIED"


# --------------------------------------------------------------------------------------
PRIMARY_CLASSES = ("SNP", "INDEL", "SUBST", "INS", "DEL",
                   "INV_PATH_EXPLICIT", "INV_DUP", "DUP", "REORDER", "UNCLASSIFIED")
_CLS_IDX = {c: i for i, c in enumerate(PRIMARY_CLASSES)}


class NodeFootprint:
    """Per-class union of graph node ids, as bitmasks over the global node id space.

    WHY NODES AND NOT REFERENCE INTERVALS
    A reference footprint answers "how much of THE REFERENCE is under variants of this
    class", which is reference-biased by construction: an insertion barely touches the
    reference, so INS scores 5.2 Mb across 469,008 alleles while DEL scores 57.8 Mb across
    249,021. The pangenome-native question is how much of the GRAPH is involved, which is a
    union over node ids and is bounded by graph length (1,926,892,905 for the clip graph)
    rather than by the length of one chosen haplotype.

    Two measures per class:
      involved  union of nodes in the ref OR alt traversal -- the graph territory the class
                occupies
      novel     union of nodes in the alt traversal and ABSENT from that record's ref
                traversal -- non-reference sequence. Deduplicated, because summing per
                allele double-counts any node novel at more than one allele.

    Bitmask rather than sets: node ids run to ~153M, so ten Python sets would be
    prohibitive. Two bytearrays hold up to 16 class bits each at 1 byte per node.
    """

    def __init__(self, max_id):
        n = max_id + 1
        self.n = n
        self.inv_lo = bytearray(n)
        self.inv_hi = bytearray(n)
        self.nov_lo = bytearray(n)
        self.nov_hi = bytearray(n)

    def add(self, cls_idx, involved_ids, novel_ids):
        bit = 1 << (cls_idx & 7)
        ilo, ihi = (self.inv_lo, self.inv_hi) if cls_idx < 8 else (self.inv_hi, self.inv_lo)
        nlo = self.nov_lo if cls_idx < 8 else self.nov_hi
        a = ilo if cls_idx < 8 else ihi
        n = self.n
        for nid in involved_ids:
            if nid < n:
                a[nid] |= bit
        for nid in novel_ids:
            if nid < n:
                nlo[nid] |= bit

    def totals(self, lengths):
        """{class: (involved_bp, novel_bp)} plus ('__ALL__', (any_involved, any_novel))."""
        if _np is None:
            return {}
        out = {}
        L = lengths
        for c, i in _CLS_IDX.items():
            bit = 1 << (i & 7)
            inv = _np.frombuffer(self.inv_lo if i < 8 else self.inv_hi, dtype=_np.uint8)
            nov = _np.frombuffer(self.nov_lo if i < 8 else self.nov_hi, dtype=_np.uint8)
            out[c] = (int(L[(inv[:len(L)] & bit) != 0].sum()),
                      int(L[(nov[:len(L)] & bit) != 0].sum()))
        ai = _np.frombuffer(self.inv_lo, dtype=_np.uint8)[:len(L)] \
            | _np.frombuffer(self.inv_hi, dtype=_np.uint8)[:len(L)]
        an = _np.frombuffer(self.nov_lo, dtype=_np.uint8)[:len(L)] \
            | _np.frombuffer(self.nov_hi, dtype=_np.uint8)[:len(L)]
        out["__ALL__"] = (int(L[ai != 0].sum()), int(L[an != 0].sum()))
        return out


def load_node_lengths(gfa):
    """Single S-line pass. Returns (numpy int32 lengths indexed by node id, max_id)."""
    if _np is None:
        sys.stderr.write("[classify] numpy unavailable; node footprints disabled\n")
        return None, 0
    with open(gfa, "rb") as probe:
        gz = probe.read(2) == b"\x1f\x8b"
    opener = gzip.open if gz else open
    # Single pass, filling a numpy array that doubles on demand. Accumulating Python lists
    # of 151.7M ids would cost several GB in int objects alone before conversion.
    cap = 1 << 24
    arr = _np.zeros(cap, dtype=_np.int32)
    max_id = 0
    n_seg = 0
    with opener(gfa, "rb") as fh:
        for line in fh:
            if line[:2] != b"S\t":
                continue
            f = line.split(b"\t", 3)
            try:
                nid = int(f[1])
            except ValueError:
                continue
            if nid >= cap:
                while nid >= cap:
                    cap <<= 1
                grown = _np.zeros(cap, dtype=_np.int32)
                grown[:len(arr)] = arr
                arr = grown
            arr[nid] = len(f[2].rstrip())
            n_seg += 1
            if nid > max_id:
                max_id = nid
    if not n_seg:
        return None, 0
    arr = arr[:max_id + 1]
    ids, lens = n_seg, None
    sys.stderr.write("[classify] node lengths: %d segments, max id %d, graph bp %d\n"
                     % (n_seg, max_id, int(arr.sum())))
    return arr, max_id


def merged_length_by_chrom(iv_by_key, want_class=None):
    """Sum of per-chromosome merged lengths. Intervals from different chromosomes must
    NEVER be merged together; they are different coordinate spaces."""
    total = 0
    for (cls, chrom), ivs in iv_by_key.items():
        if want_class is None or cls == want_class:
            total += merged_length(ivs)
    return total


def merged_length_all_chroms(iv_by_key):
    """Merged footprint pooling all classes, still per chromosome."""
    per = defaultdict(list)
    for (cls, chrom), ivs in iv_by_key.items():
        per[chrom].extend(ivs)
    return sum(merged_length(v) for v in per.values()), len(per)


def merged_length(intervals):
    if not intervals:
        return 0
    iv = sorted(intervals)
    total, s, e = 0, iv[0][0], iv[0][1]
    for a, b in iv[1:]:
        if a > e:
            total += e - s
            s, e = a, b
        elif b > e:
            e = b
    return total + (e - s)


def bins_from(spec):
    edges = [0] + [int(x) for x in spec.split(",") if x.strip()]
    edges.append(float("inf"))
    return edges


def bin_label(x, edges):
    for lo, hi in zip(edges, edges[1:]):
        if lo <= x < hi:
            return (">=%d" % lo) if hi == float("inf") else ("%d-%d" % (lo, hi))
    return "unbinned"


# --------------------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--flavor", default="clip", choices=["clip", "full"])
    p.add_argument("--tier", default="parent", choices=["parent", "fine"])
    p.add_argument("--chrom", default="all", help="for output naming only")
    p.add_argument("--minsv", type=int, default=50)
    p.add_argument("--sv-bins",
                   default="50,100,500,1000,5000,10000,50000,100000,500000,1000000")
    p.add_argument("--gfa", default=None,
                   help="GFA(.gz) for the same graph. Enables the PANGENOME-NATIVE node "
                        "footprint and non-reference node bp. Without it only the "
                        "reference-anchored footprint is emitted, which is reference-biased.")
    p.add_argument("--max-records", type=int, default=0)
    p.add_argument("--progress", type=int, default=5_000_000)
    p.add_argument("--force-topology", action="store_true",
                   help="compute topology despite decomposition markers. Produces "
                        "meaningless labels; debugging only.")
    a = p.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    stem = "%s.%s.%s.%s" % (a.label, a.flavor, a.tier, a.chrom)

    def op(suf):
        return os.path.join(a.outdir, stem + suf)

    edges = bins_from(a.sv_bins)

    node_lengths, node_max = (None, 0)
    nodefp = None
    if a.gfa:
        node_lengths, node_max = load_node_lengths(a.gfa)
        if node_lengths is not None:
            nodefp = NodeFootprint(node_max)

    topology_on = True
    decomp = []
    samples = []
    header = []

    n_rec = n_alt = 0
    lv_hist = Counter()
    at_present = at_bad = 0

    cls_n, cls_bp = Counter(), Counter()
    prim_n, prim_bp = Counter(), Counter()
    prim_bins = Counter()
    prim_bins_bp = Counter()
    label_n, label_bp = Counter(), Counter()
    cooc = Counter()
    xtab = Counter()
    # KEYED BY (class, chrom). Keying by class alone overlays every chromosome onto one
    # coordinate axis, so chr1:5,000,000 and chr7:5,000,000 merge into a single interval and
    # every footprint is silently bounded by the longest chromosome. That produced SUBST
    # 89,538,307 and all-classes 90,307,914 against chr1_1's 94.7 Mb -- and a subset of SUBST
    # (>=500 kb alleles, merged per chromosome) came to 146,413,587, exceeding the supposed
    # whole. Third instance of the same error class in this project: a footprint computed
    # without respecting the coordinate frame.
    ref_iv = defaultdict(list)        # (primary_class, chrom) -> reference intervals
    af_spec = Counter()
    ac_present = 0
    priv_sample = Counter()
    gt_called, gt_missing = Counter(), Counter()

    fsv = open(op(".sv_sizes.tsv"), "w")
    fsv.write("# longer_allele_bp is PER ALLELE. Correct for a size spectrum; NOT a total.\n"
              "# Summing it double-counts multiallelic reference spans. Use the merged\n"
              "# reference footprint in .variant_summary.tsv for totals.\n"
              "# topology_labels are NON-EXCLUSIVE; primary_class is the exclusive partition.\n")
    fsv.write("chrom\tpos\tallele_idx\tref_len\talt_len\tlonger_allele_bp\tlength_class\t"
              "primary_class\ttopology_labels\tAC\tAN\tAF\tflavor\ttier\n")

    fbed = open(op(".variants.bed"), "w")
    fbed.write('track name="%s" description="variant reference footprints"\n' % stem)

    with vopen(a.vcf) as fh:
        for line in fh:
            if line.startswith("##"):
                header.append(line)
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                samples = cols[9:] if len(cols) > 9 else []
                decomp = sorted({m for m in _DECOMP_MARKERS
                                 if any(m in h for h in header)})
                if decomp and not a.force_topology:
                    topology_on = False
                    sys.stderr.write(
                        "[classify] TOPOLOGY DISABLED -- header shows decomposition (%s).\n"
                        "           AT is inherited from the pre-decomposition parent, so\n"
                        "           alt index i does not correspond to traversal i+1.\n"
                        "           Length class and AF still computed.\n" % ", ".join(decomp))
                if a.tier == "parent" and not topology_on:
                    sys.exit("[classify] ERROR: --tier parent requires a pre-decomposition "
                             "VCF, but this one is decomposed. Wiring error.")
                continue

            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            n_rec += 1
            if a.progress and n_rec % a.progress == 0:
                sys.stderr.write("[classify] %d records\n" % n_rec)
            if a.max_records and n_rec > a.max_records:
                n_rec -= 1
                break

            chrom, pos, ref = f[0], int(f[1]), f[3]
            info = parse_info(f[7])
            lv = info.get("LV")
            lv_hist["LV=%s" % (lv if lv is not None else "missing")] += 1

            rl = len(ref)
            r0, r1 = pos - 1, pos - 1 + rl

            trav = None
            at = info.get("AT")
            if isinstance(at, str):
                at_present += 1
                if topology_on:
                    parsed = [parse_traversal(t) for t in at.split(",")]
                    if parsed and parsed[0] and all(x is not None for x in parsed):
                        trav = parsed
                    else:
                        at_bad += 1

            an = info.get("AN")
            acs = info.get("AC")
            try:
                an_i = int(an)
            except (TypeError, ValueError):
                an_i = 0
            ac_list = []
            if isinstance(acs, str):
                for x in acs.split(","):
                    try:
                        ac_list.append(int(x))
                    except ValueError:
                        ac_list.append(-1)
            if an_i and ac_list:
                ac_present += 1

            alts = f[4].split(",")
            n_alt += len(alts)
            for i, alt in enumerate(alts):
                al = len(alt)
                longer = max(rl, al)
                lc = length_class(rl, al, a.minsv)
                cls_n[lc] += 1
                cls_bp[lc] += longer

                labels = set()
                if trav is not None and i + 1 < len(trav) and lc in SV_CLASSES:
                    labels = classify_traversal(trav[0], trav[i + 1])
                    for lab in labels:
                        label_n[lab] += 1
                        label_bp[lab] += longer
                        xtab[(lc, lab)] += 1
                    cooc[",".join(sorted(labels)) or "NONE"] += 1

                pc = primary_class(labels, lc)
                if nodefp is not None and trav is not None and i + 1 < len(trav):
                    ref_ids = [n for n, _ in trav[0]]
                    alt_ids = [n for n, _ in trav[i + 1]]
                    rset = set(ref_ids)
                    nodefp.add(_CLS_IDX.get(pc, len(PRIMARY_CLASSES) - 1),
                               ref_ids + alt_ids,
                               [n for n in alt_ids if n not in rset])
                prim_n[pc] += 1
                prim_bp[pc] += longer
                _b = bin_label(longer, edges)
                prim_bins[(pc, _b)] += 1
                prim_bins_bp[(pc, _b)] += longer
                ref_iv[(pc, chrom)].append((r0, r1))

                ac = ac_list[i] if i < len(ac_list) else -1
                af = (ac / an_i) if (an_i and ac >= 0) else -1.0
                if ac >= 0 and an_i:
                    af_spec[(lc, ac)] += 1

                fsv.write("%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%.6f\t%s\t%s\n"
                          % (chrom, pos, i, rl, al, longer, lc, pc,
                             ";".join(sorted(labels)) or ".", ac, an_i, af,
                             a.flavor, a.tier))
                if lc in SV_CLASSES:
                    fbed.write("%s\t%d\t%d\t%s|AF=%.4f\t%d\t.\n"
                               % (chrom, r0, r1, pc, af if af >= 0 else -1,
                                  min(1000, longer)))

            if samples and len(f) > 9:
                fmt = f[8].split(":")
                gi = fmt.index("GT") if "GT" in fmt else None
                if gi is not None:
                    for sname, scol in zip(samples, f[9:]):
                        g = scol.split(":")
                        gt = g[gi] if gi < len(g) else "."
                        toks = gt.replace("|", "/").split("/")
                        if all(t in (".", "") for t in toks):
                            gt_missing[sname] += 1
                        else:
                            gt_called[sname] += 1
                            if any(t.isdigit() and int(t) > 0 for t in toks):
                                # attribute only when this is the sole carrier (AC==1)
                                for k, ac in enumerate(ac_list):
                                    if ac == 1 and any(t == str(k + 1) for t in toks):
                                        priv_sample[sname] += 1

    fsv.close()
    fbed.close()

    # ---- summary, with BOTH bp measures clearly distinguished -------------------------
    with open(op(".variant_summary.tsv"), "w") as out:
        out.write("# per_allele_bp: sum of max(|REF|,|ALT|) over alleles. A SIZE SPECTRUM\n")
        out.write("#   input only. NOT a total: multiallelic sites reuse one reference span\n")
        out.write("#   and max() counts ALT length when ALT is longer. On Spratelloides the\n")
        out.write("#   SV sum is 2.95 Gb against a ~1 Gb reference, 28x inflated for SUBST.\n")
        out.write("# pangenome_node_bp: union of GRAPH NODES in the ref or alt traversal.\n")
        out.write("#   THE PANGENOME-NATIVE TOTAL. Bounded by graph length, not by any one\n")
        out.write("#   haplotype's length. Requires --gfa.\n")
        out.write("# novel_node_bp: union of alt-traversal nodes ABSENT from the ref\n")
        out.write("#   traversal. Non-reference sequence, deduplicated. Requires --gfa.\n")
        out.write("# merged_ref_footprint_bp: union of REFERENCE intervals. Reference-biased\n")
        out.write("#   by construction -- an insertion barely touches the reference, so INS\n")
        out.write("#   scores far lower here than its node footprint. Kept for readers who\n")
        out.write("#   want reference-anchored numbers; prefer pangenome_node_bp.\n")
        out.write("primary_class\tn_alleles\tper_allele_bp\tpangenome_node_bp\t"
                  "novel_node_bp\tmerged_ref_footprint_bp\n")
        nt = nodefp.totals(node_lengths) if nodefp is not None else {}
        for k in sorted(prim_n, key=lambda x: -prim_n[x]):
            inv, nov = nt.get(k, (-1, -1))
            out.write("%s\t%d\t%d\t%d\t%d\t%d\n"
                      % (k, prim_n[k], prim_bp[k], inv, nov,
                         merged_length_by_chrom(ref_iv, k)))

    # separate file: two tables in one file cannot be read by read.delim, and the
    # length-based classes are legacy context rather than the primary partition
    with open(op(".length_class_summary.tsv"), "w") as out:
        out.write("# LEGACY length-based classes, kept for comparison with published runs.\n")
        out.write("# SV_COMPLEX and SV_BLOCKSUB are not variant classes -- they are the\n")
        out.write("# buckets the old length heuristic used when it could not decide. Use\n")
        out.write("# primary_class in .variant_summary.tsv instead.\n")
        out.write("length_class\tn_alleles\tper_allele_bp\n")
        for k in ("SNP", "INDEL") + SV_CLASSES:
            out.write("%s\t%d\t%d\n" % (k, cls_n[k], cls_bp[k]))
        out.write("SV\t%d\t%d\n" % (sum(cls_n[k] for k in SV_CLASSES),
                                      sum(cls_bp[k] for k in SV_CLASSES)))

    with open(op(".ref_footprint_by_chrom.tsv"), "w") as out:
        out.write("# Merged reference footprint per CHROMOSOME per class. Emitted so a wrong\n")
        out.write("# total is decomposable: if one row's footprint exceeds its chromosome's\n")
        out.write("# length, or the table has one row where it should have fifteen, the\n")
        out.write("# intervals were merged across coordinate spaces.\n")
        out.write("chrom\tprimary_class\tmerged_ref_footprint_bp\n")
        for (cls, chrom), ivs in sorted(ref_iv.items(), key=lambda x: (x[0][1], x[0][0])):
            out.write("%s\t%s\t%d\n" % (chrom, cls, merged_length(ivs)))

    with open(op(".size_spectrum.tsv"), "w") as out:
        out.write("# PRE-AGGREGATED plot input. sv_sizes.tsv has one row per allele -- 31.3M\n")
        out.write("#   on this cohort -- which is too large to read in R, so binning happens\n")
        out.write("#   here. Bin edges come from --sv-bins so they match the private-segment\n")
        out.write("#   spectrum and the two can be compared on one axis.\n")
        out.write("# per_allele_bp is the SIZE-SPECTRUM measure: correct as a y axis, NOT a\n")
        out.write("#   genome-wide total. See .variant_summary.tsv for the defensible totals.\n")
        out.write("primary_class\tsize_bin\tn_alleles\tper_allele_bp\n")
        for (k, b), c in sorted(prim_bins.items()):
            out.write("%s\t%s\t%d\t%d\n" % (k, b, c, prim_bins_bp[(k, b)]))

    with open(op(".topology_xtab.tsv"), "w") as out:
        out.write("# labels NON-EXCLUSIVE: row sums exceed the class total by design.\n")
        out.write("# INV_PATH_EXPLICIT IS A FLOOR, NOT AN INVERSION COUNT -- topology sees\n")
        out.write("# only the path-explicit motif (Romain et al. 2025). Alignment-rescued\n")
        out.write("# inversions are in SUBST; see rescue_inversions.py.\n")
        out.write("length_class\ttopology_label\tn_alleles\tper_allele_bp\n")
        for (lc, lab), c in sorted(xtab.items()):
            out.write("%s\t%s\t%d\t%d\n" % (lc, lab, c, label_bp.get(lab, 0)))

    with open(op(".label_cooccurrence.tsv"), "w") as out:
        out.write("# upset-plot input. label_set is the full non-exclusive combination.\n")
        out.write("label_set\tn_alleles\n")
        for k, c in sorted(cooc.items(), key=lambda x: -x[1]):
            out.write("%s\t%d\n" % (k, c))

    with open(op(".af_spectrum.tsv"), "w") as out:
        out.write("# AC from the VCF's own AC/AN, not recounted from sample columns.\n")
        out.write("# Sample columns have unequal ploidy (4 diploid + 1 haploid on this\n")
        out.write("# cohort, reference absent), so counting carriers is not an allele count.\n")
        out.write("length_class\tAC\tn_alleles\n")
        for (lc, ac), c in sorted(af_spec.items()):
            out.write("%s\t%d\t%d\n" % (lc, ac, c))

    with open(op(".private_variants.tsv"), "w") as out:
        out.write("# AC==1 alleles attributed to their sole carrier SAMPLE.\n")
        out.write("# NOT comparable across samples without normalising by path count:\n")
        out.write("# a diploid column has two chances to carry a singleton, a haploid one.\n")
        out.write("sample\tcalled\tmissing\tpct_missing\tprivate_alleles\n")
        for s in samples:
            c, m = gt_called[s], gt_missing[s]
            t = c + m
            out.write("%s\t%d\t%d\t%s\t%d\n"
                      % (s, c, m, ("%.4f" % (100.0 * m / t)) if t else "NA", priv_sample[s]))

    with open(op(".representation_audit.tsv"), "w") as out:
        out.write("# Emitted EVERY run. The prototype had four bugs and all four were caught\n")
        out.write("# because a printed total was impossible on its face. These are the\n")
        out.write("# quantities that made them visible; they are not debug output.\n")
        out.write("metric\tvalue\n")
        out.write("vcf\t%s\n" % os.path.basename(a.vcf))
        out.write("flavor\t%s\ntier\t%s\nchrom\t%s\n" % (a.flavor, a.tier, a.chrom))
        out.write("records\t%d\nalt_alleles\t%d\n" % (n_rec, n_alt))
        out.write("alt_alleles_per_record\t%.4f\n" % (n_alt / n_rec if n_rec else 0))
        out.write("topology_enabled\t%s\n" % topology_on)
        out.write("decomposition_markers\t%s\n" % (",".join(decomp) or "-"))
        out.write("records_with_AT\t%d\nAT_unparseable\t%d\n" % (at_present, at_bad))
        out.write("records_with_AC_AN\t%d\n" % ac_present)
        out.write("samples\t%d\n" % len(samples))
        sv_alleles = sum(cls_n[k] for k in SV_CLASSES)
        out.write("sv_alleles\t%d\n" % sv_alleles)
        out.write("sv_per_allele_bp_SUM_DO_NOT_USE_AS_TOTAL\t%d\n"
                  % sum(cls_bp[k] for k in SV_CLASSES))
        allref, n_chroms = merged_length_all_chroms(ref_iv)
        out.write("merged_ref_footprint_all_classes\t%d\n" % allref)
        out.write("chromosomes_seen\t%d\n" % n_chroms)
        out.write("# A footprint total near the length of ONE chromosome when many are\n")
        out.write("# present means intervals were merged across coordinate spaces.\n")
        if nodefp is not None:
            ai, an = nodefp.totals(node_lengths).get("__ALL__", (-1, -1))
            out.write("graph_total_bp\t%d\n" % int(node_lengths.sum()))
            out.write("pangenome_node_bp_all_classes\t%d\n" % ai)
            out.write("novel_node_bp_all_classes\t%d\n" % an)
            out.write("# CROSS-CHECK: novel_node_bp_all_classes should approximate\n")
            out.write("# graph_total_bp minus the reference path length, and should\n")
            out.write("# reconcile with the private-sequence total from gfa_hap_coverage.py.\n")
        for k, v in sorted(lv_hist.items()):
            out.write("%s\t%d\n" % (k, v))

    sys.stderr.write("\n[classify] %s\n" % stem)
    sys.stderr.write("[classify] records=%d alt_alleles=%d (%.2f per record)\n"
                     % (n_rec, n_alt, n_alt / n_rec if n_rec else 0))
    sys.stderr.write("[classify] topology %s%s\n"
                     % ("ENABLED" if topology_on else "DISABLED",
                        "" if topology_on else " (" + ",".join(decomp) + ")"))
    if topology_on:
        sys.stderr.write("[classify] primary class: %s\n"
                         % ", ".join("%s=%d" % (k, v) for k, v in
                                     sorted(prim_n.items(), key=lambda x: -x[1])))
        sys.stderr.write("[classify] INV_PATH_EXPLICIT=%d -- a FLOOR, not an inversion count\n"
                         % label_n["INV_PATH_EXPLICIT"])
    _allref, _nchrom = merged_length_all_chroms(ref_iv)
    sys.stderr.write("[classify] merged REFERENCE footprint, all classes: %d bp "
                     "across %d chromosomes\n" % (_allref, _nchrom))
    if nodefp is not None:
        ai, an = nodefp.totals(node_lengths).get("__ALL__", (-1, -1))
        sys.stderr.write("[classify] PANGENOME node footprint, all classes: %d bp\n" % ai)
        sys.stderr.write("[classify] novel (non-reference) node bp        : %d bp\n" % an)
        sys.stderr.write("[classify] graph total bp                       : %d bp\n"
                         % int(node_lengths.sum()))
    sys.stderr.write("[classify] outputs in %s\n" % a.outdir)


if __name__ == "__main__":
    main()
