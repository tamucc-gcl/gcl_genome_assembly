#!/usr/bin/env python3
# ======================================================================================
# rescue_inversions.py
#
# READ-ONLY. Finds ALIGNMENT-RESCUED inversions hiding in the SUBST bucket of a cactus
# pangenome raw VCF, and profiles whatever is left over so the unresolved-tier thresholds
# can be set from data instead of guessed.
#
# WHY THIS EXISTS
# ---------------
# Romain et al. (bioRxiv 2025.03.14.643331) describe two topological motifs for inversion
# bubbles in pangenome graphs:
#   path-explicit    ancestral and inverted alleles traverse SHARED nodes in opposite
#                    directions. Detectable from AT alone. On Spratelloides: 30 alleles.
#   alignment-rescued  the alleles are DISJOINT node sets, because the aligner never
#                    recognised the inverted copy as homologous. Topology cannot see these
#                    at all -- they look like allele replacement, i.e. SUBST. Recovering
#                    them requires aligning one allele against the REVERSE COMPLEMENT of
#                    the other. For Minigraph-Cactus the path-explicit set skews large
#                    (boundary ~500 kb), so the alignment-rescued population is expected to
#                    be both the majority and the smaller events.
#
# SUBST on Spratelloides is 1,956,621 alleles. odgi untangle independently found ~21.5 Mb
# of inverted sequence on chr10 alone against 11.4 Mb genome-wide from AT, confirming a
# large hidden population -- but untangle works at segment scale and `-e` is a resolution
# ceiling, not a discovery knob, so it does not reach the sub-500 kb events. This does.
#
# METHOD -- three stages, cheapest first
#   1. TRIAGE      stream the VCF, keep alt alleles that are SUBST by AT topology (both
#                  novel and lost nodes, no INV/DUP) and >= --min-bp.
#   2. PREFILTER   sampled k-mer strand bias. Plain (NOT canonical) k-mers, so forward and
#                  reverse-complement are distinguishable. Containment of ALT k-mers in the
#                  REF k-mer set, computed forward and revcomp. Unrelated alleles score low
#                  both ways; a true inversion scores high revcomp and low forward.
#                  Sampling stride adapts to length so cost per allele is bounded.
#   3. CONFIRM     minimap2 of revcomp(ALT) against REF, in batches. Aligned fraction is
#                  computed from the MERGED, non-overlapping query intervals -- not summed
#                  block lengths, which double-count overlapping alignments (the same error
#                  as an alignment envelope). Call at --min-frac against the SHORTER allele,
#                  matching how SUBST itself is defined and vg's own SV evaluation
#                  convention (>80% for alleles over 10 bp).
#
# WHAT IT ANSWERS
#   - how many SUBST alleles, and how much bp, are actually inversions
#   - the size spectrum of those, i.e. whether they are the sub-500 kb population predicted
#   - for everything NOT rescued: a span x ref-to-alt-ratio profile, which is the evidence
#     for pangenome_unresolved_span and pangenome_unresolved_ratio. Those are currently
#     500000 and 100 in the config and were invented, not measured.
#
# REQUIREMENTS: python3 stdlib, minimap2 on PATH. No new dependencies.
#
# USAGE
#   python3 rescue_inversions.py --vcf <label>.raw.vcf.gz --outdir rescue --label rescue
#
#   # bounded first pass to check throughput before committing
#   python3 rescue_inversions.py --vcf <label>.raw.vcf.gz --outdir /tmp/rs --label smoke \
#       --max-candidates 5000
#
# Streams the VCF once; peak memory is one record. Submit it.
# ======================================================================================

import argparse
import gzip
import os
import subprocess
import sys
import tempfile
import zlib
from collections import Counter, defaultdict

# --------------------------------------------------------------------------------------
# io
# --------------------------------------------------------------------------------------
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def vopen(path):
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    if gz:
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def parse_info(s):
    d = {}
    for kv in s.split(";"):
        if not kv:
            continue
        k, sep, v = kv.partition("=")
        d[k] = v if sep else True
    return d


_DECOMP_MARKERS = ("vcfwave", "ID=ORIGIN", "bcftools_normCommand",
                   "bcftools_normVersion", "vcfallelicprimitives")


# --------------------------------------------------------------------------------------
# AT traversal -- only enough to identify SUBST, reusing audit_variant_layers.py's logic
# --------------------------------------------------------------------------------------
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


def topology_of(ref_t, alt_t):
    """SUBST | INV_PATH_EXPLICIT | OTHER -- only the distinctions this script needs.

    INV_PATH_EXPLICIT alleles are inversions by graph topology alone. They are the natural
    POSITIVE CONTROL for the alignment test, and is_subst() excluded them by construction,
    so --control exists to route them through as well. If the alignment test does not
    recover them, its thresholds are wrong.
    """
    if _has_path_explicit_inv(ref_t, alt_t):
        return "INV_PATH_EXPLICIT"
    return "SUBST" if is_subst(ref_t, alt_t) else "OTHER"


def _has_path_explicit_inv(ref_t, alt_t):
    ref_ids = [n for n, _ in ref_t]
    alt_ids = [n for n, _ in alt_t]
    shared = set(ref_ids) & set(alt_ids)
    ro, ramb = {}, set()
    for n, o in ref_t:
        if n in ro and ro[n] != o:
            ramb.add(n)
        ro.setdefault(n, o)
    ridx = {}
    for k, n in enumerate(ref_ids):
        ridx.setdefault(n, k)
    seq = [(ridx[n], o != ro[n]) for n, o in alt_t if n in shared and n not in ramb]
    run, best, prev = 0, 0, None
    for idx, flip in seq:
        if flip and prev is not None and idx < prev:
            run = run + 1 if run else 2
        else:
            run = 1 if flip else 0
        best = max(best, run)
        prev = idx
    return best >= 2


def is_subst(ref_t, alt_t):
    """True when the alt both gains and loses node content and shows no INV or DUP.

    Deliberately the same predicate as classify_variants' primary_class == SUBST, so this
    script operates on exactly the population that classifier assigns there.
    """
    ref_ids = [n for n, _ in ref_t]
    alt_ids = [n for n, _ in alt_t]
    rc, ac = Counter(ref_ids), Counter(alt_ids)
    shared = set(rc) & set(ac)

    novel = [n for n in alt_ids if n not in rc]
    lost = [n for n in ref_ids if n not in ac]
    if not (novel and lost):
        return False
    # DUP: a node the reference also visits, traversed more times in alt
    if any(ac[n] > rc[n] for n in shared):
        return False
    # DUP_NOVEL: a novel node visited more than once
    if any(ac[n] > 1 for n in ac if n not in rc):
        return False

    # INV (path-explicit): contiguous run of >=2 shared nodes, flipped AND reversed
    ro, ramb = {}, set()
    for n, o in ref_t:
        if n in ro and ro[n] != o:
            ramb.add(n)
        ro.setdefault(n, o)
    ridx = {}
    for k, n in enumerate(ref_ids):
        ridx.setdefault(n, k)
    seq = [(ridx[n], o != ro[n]) for n, o in alt_t if n in shared and n not in ramb]
    run, best, prev = 0, 0, None
    for idx, flip in seq:
        if flip and prev is not None and idx < prev:
            run = run + 1 if run else 2
        else:
            run = 1 if flip else 0
        best = max(best, run)
        prev = idx
    return best < 2


# --------------------------------------------------------------------------------------
# stage 2 -- sampled k-mer strand bias
# --------------------------------------------------------------------------------------
def frac_kmers(seq, k, scale):
    """FracMinHash: keep k-mers whose hash is 0 mod scale.

    Both sequences of a pair MUST use the same scale, or the retained sets sample different
    positions and never intersect. That was the v1 bug: stride was derived per sequence from
    its own length, so a 20000 bp REF sampled at stride 19 while a 19000 bp ALT sampled at
    stride 18, and the sets diverged immediately. SUBST alleles differ in length by
    definition, so the prefilter failed on essentially the whole population and passed only
    the >5 Mb alleles that bypass it. Hash-based selection is position independent, so it is
    immune to that.
    """
    n = len(seq) - k + 1
    if n <= 0:
        return set()
    if scale <= 1:
        return {seq[i:i + k] for i in range(n)}
    out = set()
    for i in range(n):
        km = seq[i:i + k]
        if zlib.crc32(km.encode()) % scale == 0:
            out.add(km)
    return out


def strand_bias(ref, alt, k, target_n):
    """(forward containment, revcomp containment) of ALT k-mers in the REF k-mer set.

    Plain (non-canonical) k-mers: canonical k-mers collapse a sequence and its reverse
    complement to the same key, which would erase the signal being measured.
    """
    ref, alt = ref.upper(), alt.upper()
    scale = max(1, max(len(ref), len(alt)) // max(1, target_n))   # SHARED, from the pair
    rk = frac_kmers(ref, k, scale)
    if not rk:
        return 0.0, 0.0
    fk = frac_kmers(alt, k, scale)
    vk = frac_kmers(revcomp(alt), k, scale)
    f = (len(fk & rk) / len(fk)) if fk else 0.0
    v = (len(vk & rk) / len(vk)) if vk else 0.0
    return f, v


# --------------------------------------------------------------------------------------
# stage 3 -- minimap2 confirmation, merged intervals
# --------------------------------------------------------------------------------------
def merged_len(intervals):
    """Union length of [start, end) intervals. Summed block lengths would double-count
    overlapping alignments -- the same error as using an alignment envelope."""
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


def run_minimap2(pairs, threads, preset, workdir, keep):
    """pairs: list of (id, ref_seq, alt_seq). Returns {id: (fwd_merged_bp, rev_merged_bp)}.

    Returns four merged coverages per id: query-forward, query-reverse, target-forward,
    target-reverse. Coverage must be read in the SAME frame as the denominator it is divided
    by, or the fraction can exceed 1 -- v3 merged query intervals but divided by
    min(ref_len, alt_len), producing frac_rev up to 2.8 whenever ALT was the longer allele.

    ALT is aligned to REF AS-IS -- no reverse complement. minimap2 searches both strands by
    default, so pre-reverse-complementing the query is a no-op: it finds the same homology
    and merely reports it on the opposite strand. v2 did exactly that and then counted any
    alignment as evidence of inversion, which called 66% of SUBST and produced 1.74 Gb of
    "inverted" sequence in a 1 Gb genome.

    PAF column 5 IS the signal. Coverage is merged separately for '+' and '-' rows; an
    inversion is homology predominantly on the minus strand.

    Merging matters too: summed block lengths would double-count overlapping alignments,
    the same error as using an alignment envelope.
    """
    if not pairs:
        return {}
    tgt = os.path.join(workdir, "t.fa")
    qry = os.path.join(workdir, "q.fa")
    with open(tgt, "w") as ft, open(qry, "w") as fq:
        for pid, rseq, aseq in pairs:
            ft.write(">%s\n%s\n" % (pid, rseq))
            fq.write(">%s\n%s\n" % (pid, aseq))
    cmd = ["minimap2", "-x", preset, "-t", str(threads), "--secondary=yes",
           "-N", "50", tgt, qry]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    except FileNotFoundError:
        sys.exit("ERROR: minimap2 not on PATH")
    if proc.returncode != 0:
        sys.stderr.write("[rescue] minimap2 rc=%d: %s\n"
                         % (proc.returncode, proc.stderr[-400:]))
        return {}
    qf, qr, tf, tr = (defaultdict(list) for _ in range(4))
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 9:
            continue
        if f[0] != f[5]:          # cross-hit between different candidates in the batch
            continue
        minus = (f[4] == "-")
        # PAF query and target coordinates are both on the forward strand of their own
        # sequence regardless of alignment strand, so both can be merged directly.
        (qr if minus else qf)[f[0]].append((int(f[2]), int(f[3])))
        (tr if minus else tf)[f[0]].append((int(f[7]), int(f[8])))
    if not keep:
        for pth in (tgt, qry):
            try:
                os.unlink(pth)
            except OSError:
                pass
    ids = set(qf) | set(qr) | set(tf) | set(tr)
    return {pid: (merged_len(qf.get(pid, [])), merged_len(qr.get(pid, [])),
                  merged_len(tf.get(pid, [])), merged_len(tr.get(pid, []))) for pid in ids}


# --------------------------------------------------------------------------------------
BIN_EDGES = [0, 1_000, 5_000, 10_000, 50_000, 100_000, 500_000,
             1_000_000, 10_000_000, float("inf")]


def bin_label(x):
    for lo, hi in zip(BIN_EDGES, BIN_EDGES[1:]):
        if lo <= x < hi:
            return (">=%d" % lo) if hi == float("inf") else ("%d-%d" % (lo, hi))
    return "unbinned"


def ratio_bin(r):
    for t in (1.5, 3, 10, 100, 1000):
        if r < t:
            return "<%g" % t
    return ">=1000"


# --------------------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf", required=True, help="raw (PRE-decomposition) VCF")
    p.add_argument("--outdir", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--min-bp", type=int, default=1000,
                   help="size floor on the longer allele (default 1000)")
    p.add_argument("--minsv", type=int, default=50)
    p.add_argument("--kmer", type=int, default=21)
    p.add_argument("--kmer-sample", type=int, default=1000,
                   help="k-mers sampled per sequence; bounds cost per allele (default 1000)")
    p.add_argument("--kmer-max-bp", type=int, default=5_000_000,
                   help="alleles longer than this skip the prefilter and go straight to "
                        "minimap2 (default 5Mb)")
    p.add_argument("--prefilter-rev", type=float, default=0.20,
                   help="minimum revcomp containment to reach minimap2 (default 0.20)")
    p.add_argument("--prefilter-ratio", type=float, default=2.0,
                   help="revcomp containment must exceed forward by this factor (default 2)")
    p.add_argument("--min-strand-margin", type=float, default=2.0,
                   help="minus-strand fraction must exceed plus-strand by this factor "
                        "(default 2.0). Without it, ordinary forward homology is called an "
                        "inversion, because minimap2 reports both strands.")
    p.add_argument("--min-frac", type=float, default=0.8,
                   help="merged aligned fraction of the SHORTER allele to call (default 0.8)")
    p.add_argument("--preset", default="asm20")
    p.add_argument("--batch", type=int, default=200,
                   help="max pairs per minimap2 invocation (default 200)")
    p.add_argument("--batch-bp", type=int, default=20_000_000,
                   help="max total bp per minimap2 invocation; a few megabase alleles "
                        "must not blow up one batch (default 20Mb)")
    p.add_argument("--prefilter", action="store_true",
                   help="enable the k-mer strand-bias prefilter. OFF by default: at this "
                        "scale minimap2 on everything is affordable and cannot introduce a "
                        "sampling bias. Enable only for much larger cohorts.")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--max-candidates", type=int, default=0,
                   help="stop after N SUBST candidates (0 = all); throughput testing")
    p.add_argument("--progress", type=int, default=5_000_000)
    p.add_argument("--subst-homologous", type=float, default=0.5,
                   help="max(frac_fwd, frac_rev) at or above which a non-inverted SUBST is "
                        "SUBST_HOMOLOGOUS: alleles that DO align, so the graph fragmented a "
                        "homologous region into novel-plus-lost node content (default 0.5). "
                        "On Spratelloides 4,015 of 4,329 SUBST alleles >=100 kb are here.")
    p.add_argument("--subst-unrelated", type=float, default=0.01,
                   help="max(frac_fwd, frac_rev) below which a SUBST is SUBST_UNRELATED: "
                        "genuinely unrelated sequence at one locus. Concentrated in SMALL "
                        "alleles -- 22.9% at 1-10 kb, only 3.5% above 1 Mb -- which is why "
                        "the unresolved tier is NOT a span threshold (default 0.01).")
    p.add_argument("--control", action="store_true",
                   help="also run the alignment test on INV_PATH_EXPLICIT alleles, which are "
                        "inversions by topology alone. Positive control: if these are not "
                        "recovered, the thresholds are too strict.")
    p.add_argument("--keep-temp", action="store_true")
    a = p.parse_args()

    os.makedirs(a.outdir, exist_ok=True)

    def op(suf):
        return os.path.join(a.outdir, a.label + suf)

    n_rec = 0
    n_subst = 0
    n_control = 0
    n_cand = 0
    n_prefilter_pass = 0
    n_called = 0
    called_bp = 0

    size_hist = Counter()          # rescued, by size bin
    subst_hist = Counter()         # all SUBST >= min_bp, by size bin
    subst_hist_cls = Counter()     # SUBST subdivided by measured homology
    subst_bp_cls = Counter()
    unresolved_profile = Counter() # (size bin, ratio bin) for NON-rescued
    header_lines = []

    batch = []
    batch_bp = 0
    called_aln_bp = 0
    control_hit = [0]
    tmpdir = tempfile.mkdtemp(prefix="rescue_", dir=a.outdir)

    fcand = open(op(".rescue_candidates.tsv"), "w")
    fcand.write("# every SUBST allele >= --min-bp, with prefilter and confirmation outcome\n")
    fcand.write("# aligned_bp is the MERGED query interval union, not summed block lengths\n")
    fcand.write("# aligned_*_bp are MERGED query interval unions per strand, not summed\n")
    fcand.write("# block lengths. frac_* are against the SHORTER allele.\n")
    fcand.write("# q_* are merged coverage on the ALT allele, t_* on the REF allele.\n")
    fcand.write("# frac_* use whichever frame belongs to the SHORTER allele.\n")
    fcand.write("# topology: SUBST = candidate population; INV_PATH_EXPLICIT = positive\n")
    fcand.write("#   control, already an inversion by AT, included only with --control.\n")
    fcand.write("chrom\tpos\tallele_idx\tref_len\talt_len\tlonger_bp\tratio\t"
                "kmer_fwd\tkmer_rev\tq_fwd_bp\tq_rev_bp\tt_fwd_bp\tt_rev_bp\t"
                "frac_fwd\tfrac_rev\ttopology\tcall\n")

    fsub = open(op(".subst_homology.tsv"), "w")
    fsub.write("# SUBST subdivided by MEASURED homology between its own alleles, replacing\n")
    fsub.write("# the span/ratio thresholds originally planned for the unresolved tier.\n")
    fsub.write("# Those were wrong: 93% of SUBST alleles >=100 kb align to their reference\n")
    fsub.write("# allele at >=0.5 with no inversion, so large SUBST is NOT an alignment hole\n")
    fsub.write("# -- it is homologous sequence the GRAPH chose to represent as replacement.\n")
    fsub.write("#   SUBST_HOMOLOGOUS  max frac >= %.2f : graph representation artifact\n"
               % a.subst_homologous)
    fsub.write("#   SUBST_PARTIAL     between            : mixed\n")
    fsub.write("#   SUBST_UNRELATED   max frac <  %.2f : genuinely unrelated sequence\n"
               % a.subst_unrelated)
    fsub.write("# INV_ALN_RESCUED rows are in .rescue_calls.bed, not here.\n")
    fsub.write("chrom\tpos\tallele_idx\tlonger_bp\tfrac_fwd\tfrac_rev\tsubst_class\n")

    fbed = open(op(".rescue_calls.bed"), "w")
    fbed.write('track name="INV_ALN_RESCUED" description="alignment-rescued inversions"\n')

    def flush(force=False):
        nonlocal batch, batch_bp, n_called, called_bp, called_aln_bp
        if not batch:
            return
        if not force and len(batch) < a.batch and batch_bp < a.batch_bp:
            return
        res = run_minimap2([(b["id"], b["ref"], b["alt"]) for b in batch],
                           a.threads, a.preset, tmpdir, a.keep_temp)
        for b in batch:
            qfw, qrv, tfw, trv = res.get(b["id"], (0, 0, 0, 0))
            # measure the SHORTER allele's coverage in its OWN coordinate frame
            if b["alt_len"] <= b["ref_len"]:
                afwd, arev, denom = qfw, qrv, b["alt_len"] or 1
            else:
                afwd, arev, denom = tfw, trv, b["ref_len"] or 1
            ffwd, frev = min(1.0, afwd / denom), min(1.0, arev / denom)
            # inversion: minus-strand homology clears the floor AND dominates plus-strand.
            # The margin matters because repeats give real inversions some forward homology.
            call = ("INV_ALN_RESCUED"
                    if (frev >= a.min_frac and frev >= a.min_strand_margin * max(ffwd, 1e-9))
                    else "SUBST")
            if call == "INV_ALN_RESCUED" and b["src"] == "INV_PATH_EXPLICIT":
                control_hit[0] += 1
            if call == "INV_ALN_RESCUED" and b["src"] != "SUBST":
                continue          # control alleles are reported, not added to the catalog
            if call == "INV_ALN_RESCUED":
                n_called += 1
                called_bp += b["longer"]
                called_aln_bp += arev
                size_hist[bin_label(b["longer"])] += 1
                fbed.write("%s\t%d\t%d\tINV_ALN_RESCUED_%s\t%d\t-\n"
                           % (b["chrom"], b["pos"] - 1,
                              b["pos"] - 1 + b["ref_len"], b["id"],
                              int(min(1000, frev * 1000))))
            else:
                unresolved_profile[(bin_label(b["longer"]), ratio_bin(b["ratio"]))] += 1
                mx = max(ffwd, frev)
                sc = ("SUBST_HOMOLOGOUS" if mx >= a.subst_homologous
                      else "SUBST_UNRELATED" if mx < a.subst_unrelated
                      else "SUBST_PARTIAL")
                subst_hist_cls[sc] += 1
                subst_bp_cls[sc] += b["longer"]
                if b["src"] == "SUBST":
                    fsub.write("%s\t%d\t%d\t%d\t%.4f\t%.4f\t%s\n"
                               % (b["chrom"], b["pos"], b["ai"], b["longer"],
                                  ffwd, frev, sc))
            fcand.write("%s\t%d\t%d\t%d\t%d\t%d\t%.3f\t%.4f\t%.4f\t"
                        "%d\t%d\t%d\t%d\t%.4f\t%.4f\t%s\t%s\n"
                        % (b["chrom"], b["pos"], b["ai"], b["ref_len"], b["alt_len"],
                           b["longer"], b["ratio"], b["kf"], b["kr"],
                           qfw, qrv, tfw, trv, ffwd, frev, b["src"], call))
        batch = []
        batch_bp = 0

    with vopen(a.vcf) as fh:
        for line in fh:
            if line.startswith("##"):
                header_lines.append(line)
                continue
            if line.startswith("#CHROM"):
                bad = [m for m in _DECOMP_MARKERS if any(m in h for h in header_lines)]
                if bad:
                    sys.exit("ERROR: this VCF shows decomposition (%s). AT is inherited from\n"
                             "       the parent record, so allele i does not correspond to\n"
                             "       traversal i+1 and SUBST cannot be identified. Use the\n"
                             "       raw VCF." % ", ".join(bad))
                continue

            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            n_rec += 1
            if a.progress and n_rec % a.progress == 0:
                sys.stderr.write("[rescue] %d rec | %d SUBST>=%d | %d aligned | %d called\n"
                                 % (n_rec, n_subst, a.min_bp, n_prefilter_pass, n_called))

            info = parse_info(f[7])
            at = info.get("AT")
            if not isinstance(at, str):
                continue
            trav = [parse_traversal(t) for t in at.split(",")]
            if not trav or trav[0] is None or any(t is None for t in trav):
                continue

            chrom, pos, ref = f[0], int(f[1]), f[3]
            rl = len(ref)
            for ai, alt in enumerate(f[4].split(",")):
                al = len(alt)
                longer = max(rl, al)
                if longer < a.min_bp:
                    continue
                d, mn = abs(rl - al), min(rl, al)
                if rl == 1 and al == 1:
                    continue
                if d < a.minsv and mn < a.minsv:
                    continue
                if ai + 1 >= len(trav):
                    continue
                topo = topology_of(trav[0], trav[ai + 1])
                if topo == "OTHER":
                    continue
                if topo == "INV_PATH_EXPLICIT" and not a.control:
                    continue
                if topo == "SUBST":
                    n_subst += 1
                else:
                    n_control += 1
                if topo == "SUBST":
                    subst_hist[bin_label(longer)] += 1
                ratio = (rl / al) if al else float("inf")

                if a.prefilter and longer <= a.kmer_max_bp:
                    kf, kr = strand_bias(ref, alt, a.kmer, a.kmer_sample)
                    if not (kr >= a.prefilter_rev and kr >= a.prefilter_ratio * max(kf, 1e-9)):
                        unresolved_profile[(bin_label(longer), ratio_bin(ratio))] += 1
                        fcand.write("%s\t%d\t%d\t%d\t%d\t%d\t%.3f\t%.4f\t%.4f\t"
                                "0\t0\t0\t0\t0.0000\t0.0000\tSUBST\tSUBST\n"
                                % (chrom, pos, ai, rl, al, longer, ratio, kf, kr))
                        continue
                else:
                    kf, kr = -1.0, -1.0     # prefilter disabled or allele too long

                n_prefilter_pass += 1
                batch.append({"id": "c%d" % n_prefilter_pass, "chrom": chrom, "pos": pos,
                              "ai": ai, "ref_len": rl, "alt_len": al, "longer": longer,
                              "ratio": ratio, "kf": kf, "kr": kr,
                              "ref": ref.upper(), "alt": alt.upper(), "src": topo})
                batch_bp += rl + al
                flush()

                n_cand += 1
                if a.max_candidates and n_cand >= a.max_candidates:
                    flush(force=True)
                    sys.stderr.write("[rescue] --max-candidates reached\n")
                    n_rec = -n_rec
                    break
            if n_rec < 0:
                n_rec = -n_rec
                break

    flush(force=True)
    fcand.close()
    fbed.close()
    fsub.close()

    with open(op(".rescue_summary.tsv"), "w") as out:
        out.write("metric\tvalue\n")
        out.write("records_scanned\t%d\n" % n_rec)
        out.write("subst_alleles_ge_min_bp\t%d\n" % n_subst)
        out.write("control_path_explicit_tested\t%d\n" % n_control)
        out.write("control_path_explicit_recovered\t%d\n" % control_hit[0])
        out.write("reached_minimap2\t%d\n" % n_prefilter_pass)
        out.write("called_INV_ALN_RESCUED\t%d\n" % n_called)
        out.write("called_bp_longer_allele\t%d\n" % called_bp)
        out.write("called_bp_merged_aligned_minus_strand\t%d\n" % called_aln_bp)
        out.write("min_strand_margin\t%.2f\n" % a.min_strand_margin)
        out.write("\nsubst_class\tn_alleles\tper_allele_bp\n")
        for k in ("SUBST_HOMOLOGOUS", "SUBST_PARTIAL", "SUBST_UNRELATED"):
            out.write("%s\t%d\t%d\n" % (k, subst_hist_cls[k], subst_bp_cls[k]))
        out.write("min_bp\t%d\n" % a.min_bp)
        out.write("min_frac\t%.3f\n" % a.min_frac)
        out.write("\nsize_bin\tsubst_total\trescued\n")
        for b in sorted(set(subst_hist) | set(size_hist),
                        key=lambda s: BIN_EDGES.index(int(s.split("-")[0].lstrip(">=")))
                        if s.split("-")[0].lstrip(">=").isdigit() else 99):
            out.write("%s\t%d\t%d\n" % (b, subst_hist[b], size_hist[b]))

    with open(op(".rescue_unresolved_profile.tsv"), "w") as out:
        out.write("# NON-rescued SUBST alleles, by size and ref-to-alt ratio.\n")
        out.write("# THIS SETS pangenome_unresolved_span AND pangenome_unresolved_ratio.\n")
        out.write("# Look for where the distribution separates: a dense low-ratio mass is\n")
        out.write("# ordinary allele replacement; a sparse high-ratio or very-large-span\n")
        out.write("# tail is a path skipping the region and belongs in the unresolved tier.\n")
        out.write("size_bin\tref_to_alt_ratio_bin\tn_alleles\n")
        for (sb, rb), c in sorted(unresolved_profile.items()):
            out.write("%s\t%s\t%d\n" % (sb, rb, c))

    if not a.keep_temp:
        try:
            os.rmdir(tmpdir)
        except OSError:
            pass

    sys.stderr.write("\n[rescue] ===== headline =====\n")
    sys.stderr.write("[rescue] SUBST >= %d bp        : %d\n" % (a.min_bp, n_subst))
    sys.stderr.write("[rescue] passed prefilter      : %d (%.2f%%)\n"
                     % (n_prefilter_pass, 100.0 * n_prefilter_pass / n_subst if n_subst else 0))
    if n_control:
        sys.stderr.write("[rescue] POSITIVE CONTROL: %d/%d path-explicit inversions recovered"
                         " by alignment\n" % (control_hit[0], n_control))
    sys.stderr.write("[rescue] called INV_ALN_RESCUED: %d\n" % n_called)
    sys.stderr.write("[rescue]   longer-allele bp     : %d\n" % called_bp)
    sys.stderr.write("[rescue]   minus-strand aln bp  : %d  <- the inverted homologous seq\n"
                     % called_aln_bp)
    sys.stderr.write("[rescue] SUBST subdivision: %s\n"
                     % ", ".join("%s=%d" % (k, subst_hist_cls[k])
                                 for k in ("SUBST_HOMOLOGOUS", "SUBST_PARTIAL",
                                           "SUBST_UNRELATED")))
    sys.stderr.write("[rescue] reports in %s\n" % a.outdir)


if __name__ == "__main__":
    main()
