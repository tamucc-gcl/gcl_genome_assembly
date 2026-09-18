#!/usr/bin/env python3
# ======================================================================================
# chimera_evidence.py
#
# Independent evidence for ONE chimeric-scaffold candidate: Hi-C cross-contact profile,
# interstitial and terminal telomere signal, nearest N-gap, and a refined cut point.
#
# WHAT THIS IS FOR
# ----------------
# harmonize_names.py detects composites and votes on them across haplotypes, which is the
# only signal independent of how the scaffold was BUILT. That vote is the gate. This adds the
# evidence that locates the cut precisely and lets a human confirm the call:
#
#   Hi-C depletion   locates the junction to ~100 kb, far better than the PAF's 1 Mb bins
#   telomeres        independent of Hi-C; diagnostic of END-TO-END fusion specifically
#   N-gap            the cut point that severs nothing
#
# THE HI-C CAVEAT, STATED ONCE
# ----------------------------
# Hi-C made the join, so it cannot independently justify breaking it. It is used to LOCATE,
# not to justify. And it must be measured correctly: the first attempt on the chr5+chr9
# scaffold returned 1.209 (elevated) from measuring at 35 Mb -- a 5 Mb-bin PAF estimate --
# with edge-truncated windows. Scanned properly with full windows only, the minimum is at
# 37.4 Mb with a ratio of 0.748. Depleted. A surprising number got explained rather than
# re-checked, and a design decision was built on the artifact.
#
# So: FULL WINDOWS ONLY, and SCAN for the minimum rather than assuming its position.
#
# TELOMERE ABSENCE IS NOT EVIDENCE AGAINST
# ----------------------------------------
# A mid-arm fusion leaves no interstitial array. Measured: chr5_1+chr9_1 has arrays in BOTH
# orientations over 37.6-43.7 Mb peaking at 148 while its own termini are near-silent --
# two chromosome ends fused back to back. chr6_3+chr12_1 has nothing at its junction but
# fwd=353 at position 0 -- an intact chromosome arm from its own telomere, joined to a
# truncated partner. Both are real chimeras; only one has telomere corroboration.
#
# Dependencies: cooler, numpy, matplotlib.
#
# USAGE
#   chimera_evidence.py --cool <per-scaffold .cool> --fasta <mini reference .fa> \
#       --candidates <chimera_candidates.tsv> --assembly <id> --scaffold <name> \
#       --telomere-motif CCCTAA --outdir . [--label <prefix>]
# ======================================================================================

import argparse
import os
import sys

import numpy as np


# --------------------------------------------------------------------------------------
# inputs
# --------------------------------------------------------------------------------------
def read_rows(path):
    """'#'-commented TSV -> list of dicts.

    Comments stripped by hand, comment.char never used: '#' is a legitimate character inside
    PanSN haplotype names elsewhere in this pipeline and using it as a comment marker has
    silently truncated keys twice.
    """
    rows, hdr = [], None
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


def read_one_fasta(path):
    """The mini reference holds exactly one record; return (name, sequence)."""
    name, buf = None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    break
                name = line[1:].strip().split()[0]
            else:
                buf.append(line.strip())
    return name, "".join(buf)


# --------------------------------------------------------------------------------------
# Hi-C
# --------------------------------------------------------------------------------------
def hic_profile(cool_path, scaffold, win_bins=20):
    """Cross-contact for every position, FULL WINDOWS ONLY.

    At position a the statistic is the mean of the block [a-w, a) x [a, a+w) -- contact
    BETWEEN the two sides. A real junction depletes it; ordinary sequence does not.

    Positions closer than w to either end are skipped rather than truncated. Truncating them
    (`max(0, a-w)`) averages a smaller block and produced a spurious minimum 2 Mb from the
    start of a scaffold whose first 10 kb carry a genuine telomere.

    cooler's region parser CANNOT be used: it reads the text after '+' in `chr5_1+chr9_1` as
    a coordinate and raises. Bins are located by index instead.
    """
    import cooler
    c = cooler.Cooler(cool_path)
    res = c.binsize
    bins = c.bins()[:]
    names = bins["chrom"].astype(str).values
    idx = np.where(names == scaffold)[0]
    if idx.size == 0:
        return None, res, "scaffold %r not in the cool (have: %s)" % (
            scaffold, ", ".join(sorted(set(names))[:5]))
    i0, i1 = int(idx.min()), int(idx.max()) + 1
    M = c.matrix(balance=False, sparse=False)[i0:i1, i0:i1].astype(float)
    n = M.shape[0]
    if n < 4 * win_bins:
        return None, res, "only %d bins at %d bp; need >= %d for full windows" % (
            n, res, 4 * win_bins)

    prof = np.full(n, np.nan)
    for a in range(win_bins, n - win_bins):
        prof[a] = np.nanmean(M[a - win_bins:a, a:a + win_bins])
    return prof, res, ""


def summarise_profile(prof, res, n_low=5):
    """Minimum position, its ratio to the median, and how many of the n_low lowest windows
    are CONTIGUOUS with the minimum.

    Contiguity is the discriminating statistic. Five consecutive low windows is a boundary;
    one isolated low window is noise. Measured: chr5_1+chr9_1 had all five lowest consecutive
    over 37.2-37.6 Mb; chr6_3+chr12_1 had three of five over 44.2-44.5 Mb.
    """
    ok = np.isfinite(prof)
    if not ok.any():
        return None
    med = float(np.nanmedian(prof[ok]))
    # kind="stable": the default quicksort leaves the order among TIES arbitrary, and a flat
    # low region is common in sparse bins. Without this, which window counts as "the
    # minimum" varies between runs on the same data.
    order = np.argsort(np.where(ok, prof, np.inf), kind="stable")
    lowest = [int(x) for x in order[:n_low]]
    mn = lowest[0]
    # The LONGEST contiguous run among the lowest windows -- not the run through whichever
    # one sorted first. With ties that choice is arbitrary even under a stable sort (the
    # lowest value may appear at either end of a flat region), so a run measured from it
    # varies. The question being asked is "is there a contiguous block of low windows", and
    # the longest run answers it regardless of which tie leads.
    low = sorted(set(lowest))
    contiguous, run = 1, 1
    for i in range(1, len(low)):
        run = run + 1 if low[i] == low[i - 1] + 1 else 1
        contiguous = max(contiguous, run)
    return {
        "min_bin": mn,
        "min_bp": mn * res,
        "min_value": float(prof[mn]),
        "median": med,
        "ratio": (float(prof[mn]) / med) if med else float("nan"),
        "n_low_contiguous": contiguous,
        "low_bp": [int(k) * res for k in sorted(lowest)],
    }


# --------------------------------------------------------------------------------------
# telomeres, computed here rather than taken from tidk
# --------------------------------------------------------------------------------------
def revcomp(s):
    return s.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def read_tidk_windows(path, scaffold):
    """tidk search output -> [(window_start, forward, reverse), ...] for one scaffold.

    THE SOURCE OF RECORD. tidk is already a pipeline dependency and its output is what
    FINAL_VIZ reports, so computing the same quantity a second way means that if the two ever
    disagree there is no way to know which is right. It also normalises the canonical repeat,
    which a plain motif count does not.

    tidk runs in FINAL_VIZ, downstream of harmonization, so its output does not exist when
    this script runs in the pipeline -- CHIMERA_EVIDENCE invokes tidk itself on the
    single-scaffold mini reference and passes the result here.

    Columns: id / window / forward_repeat_number / reverse_repeat_number / telomeric_repeat
    """
    rows, hdr = [], None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.lstrip().startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = f
                continue
            if len(f) != len(hdr):
                continue
            d = dict(zip(hdr, f))
            if d.get("id") != scaffold:
                continue
            try:
                rows.append((int(d["window"]), int(d["forward_repeat_number"]),
                             int(d["reverse_repeat_number"])))
            except (KeyError, ValueError):
                continue
    rows.sort()
    # tidk reports the window END; convert to START so positions match the sequence-derived
    # fallback and the figure's x axis
    if len(rows) >= 2:
        w = rows[1][0] - rows[0][0]
        rows = [(max(0, p - w), f, r) for p, f, r in rows]
    return rows


def telomere_windows(seq, motif, win=10000):
    """FALLBACK ONLY: motif count per window in both orientations, from the sequence.

    Used when no tidk output is supplied. Equivalent for CCCTAA -- that motif has no
    prefix/suffix self-overlap, so non-overlapping str.count is exact -- but it does not
    normalise the canonical repeat, so for an unusual motif it can diverge from tidk. The
    audit records which source was used.
    """
    fwd, rev = motif.upper(), revcomp(motif.upper())
    up = seq.upper()
    out = []
    for s in range(0, len(up), win):
        chunk = up[s:s + win]
        out.append((s, chunk.count(fwd), chunk.count(rev)))
    return out


def summarise_telomere(tw, junction_bp, flank=2000000, term=100000):
    """Interstitial signal near the junction against the scaffold background, plus the
    TERMINAL signal at both ends.

    The terminal reading is what distinguishes the two measured cases: chr6_3+chr12_1 has
    fwd=353 in its first window, so that arm runs from a real chromosome end -- an intact
    chromosome joined to a truncated partner. chr5_1+chr9_1 has near-silent termini and a
    148-peak array in the middle -- two ends fused, with the scaffold's own ends incomplete.
    """
    if not tw:
        return None
    tot = [f + r for _, f, r in tw]
    bg = float(np.median(tot)) if tot else 0.0
    near = [(p, f, r) for p, f, r in tw if abs(p - junction_bp) <= flank]
    span = tw[-1][0] + 10000
    # A scaffold shorter than 3x the terminal window has no distinguishable "ends", and the
    # head/tail slices would overlap the middle -- reporting the interstitial array as
    # terminal signal. Candidates are >= 20 Mb so this never fires in practice, but without
    # it the readings are silently wrong on anything small.
    if span < 3 * term:
        head = tail = []
    else:
        head = [(p, f, r) for p, f, r in tw if p < term]
        tail = [(p, f, r) for p, f, r in tw if p > span - term]
    mx = max(near, key=lambda x: x[1] + x[2]) if near else (0, 0, 0)
    return {
        "background_median": bg,
        "junction_max": mx[1] + mx[2],
        "junction_max_bp": mx[0],
        "junction_max_fwd": mx[1],
        "junction_max_rev": mx[2],
        "junction_over_background": ((mx[1] + mx[2]) / bg) if bg else float("nan"),
        "both_orientations": bool(mx[1] > 0 and mx[2] > 0),
        "terminal_start": max((f + r for _, f, r in head), default=0),
        "terminal_end": max((f + r for _, f, r in tail), default=0),
    }


# --------------------------------------------------------------------------------------
# N-gaps
# --------------------------------------------------------------------------------------
def gaps(seq, min_bp=50):
    out, st = [], None
    up = seq.upper()
    for i, c in enumerate(up):
        if c == "N":
            if st is None:
                st = i
        elif st is not None:
            if i - st >= min_bp:
                out.append((st, i))
            st = None
    if st is not None and len(up) - st >= min_bp:
        out.append((st, len(up)))
    return out


def snap(gs, target, window=500000):
    """Nearest gap to the target, if one is within `window`.

    Cutting inside a gap severs nothing; cutting in sequence severs real bases. Which of
    those happened has to be recorded, because a reviewer cannot tell otherwise -- and after
    gap filling the answer varies per junction. Measured: chr6_3+chr12_1 has a gap 8 kb from
    its Hi-C minimum, chr5_1+chr9_1's nearest is 114 kb away.
    """
    if not gs:
        return None
    best = min(gs, key=lambda g: abs((g[0] + g[1]) // 2 - target))
    mid = (best[0] + best[1]) // 2
    if abs(mid - target) > window:
        return None
    return {"gap_start": best[0], "gap_end": best[1], "cut_bp": mid,
            "offset_from_target": mid - target}


# --------------------------------------------------------------------------------------
# figure
# --------------------------------------------------------------------------------------
def figure(path, cool_path, scaffold, res, prof, paf_bp, hic_bp, cut_bp, tw, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import cooler

    c = cooler.Cooler(cool_path)
    bins = c.bins()[:]
    idx = np.where(bins["chrom"].astype(str).values == scaffold)[0]
    i0, i1 = int(idx.min()), int(idx.max()) + 1
    M = c.matrix(balance=False, sparse=False)[i0:i1, i0:i1].astype(float)
    n = M.shape[0]
    mb = n * res / 1e6

    fig, ax = plt.subplots(3, 1, figsize=(9, 11),
                           gridspec_kw={"height_ratios": [4, 1.1, 1.1]})
    ax[0].imshow(np.log10(M + 1), cmap="YlOrRd", extent=[0, mb, mb, 0])
    for v, col, lab in ((paf_bp, "tab:blue", "PAF"), (hic_bp, "tab:green", "Hi-C min"),
                        (cut_bp, "black", "cut")):
        if v is not None:
            ax[0].axhline(v / 1e6, color=col, lw=0.8, ls="--")
            ax[0].axvline(v / 1e6, color=col, lw=0.8, ls="--", label=lab)
    ax[0].set_title(title, fontsize=10)
    ax[0].set_xlabel("Mb"); ax[0].set_ylabel("Mb")
    ax[0].legend(fontsize=7, loc="upper right")

    x = np.arange(n) * res / 1e6
    ax[1].plot(x, prof, lw=0.7, color="tab:purple")
    ax[1].axhline(np.nanmedian(prof), color="grey", lw=0.6, ls=":")
    for v, col in ((paf_bp, "tab:blue"), (hic_bp, "tab:green"), (cut_bp, "black")):
        if v is not None:
            ax[1].axvline(v / 1e6, color=col, lw=0.8, ls="--")
    ax[1].set_ylabel("cross-contact")
    ax[1].set_xlim(0, mb)

    if tw:
        tp = np.array([p for p, _, _ in tw]) / 1e6
        ax[2].plot(tp, [f for _, f, _ in tw], lw=0.6, color="tab:red", label="forward")
        ax[2].plot(tp, [-r for _, _, r in tw], lw=0.6, color="tab:blue", label="reverse")
        ax[2].axhline(0, color="grey", lw=0.5)
        for v, col in ((cut_bp, "black"),):
            if v is not None:
                ax[2].axvline(v / 1e6, color=col, lw=0.8, ls="--")
        ax[2].set_ylabel("telomere motif")
        ax[2].legend(fontsize=7)
        ax[2].set_xlim(0, mb)
    ax[2].set_xlabel("position (Mb)")

    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()


# --------------------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cool", default="",
                   help="per-scaffold .cool at ~100 kb. OPTIONAL: telomere and N-gap evidence "
                        "are independent of Hi-C and cheap, while a contact map at detection "
                        "time requires re-mapping reads to a mini reference (the scaffold-"
                        "stage BAM is in contig or round-1 space, and round-1/round-2 "
                        "scaffold names COLLIDE). Without it those two still run and the "
                        "Hi-C panel is skipped.")
    p.add_argument("--fasta", required=True, help="mini reference: this scaffold only")
    p.add_argument("--candidates", required=True,
                   help="chimera_joins.py output: the called joins, with cut_bp from the AGP")
    p.add_argument("--cut-bp", type=int, default=0,
                   help="confirm THIS cut rather than re-deriving one. The AGP already gives "
                        "the position exactly; this script confirms the decision, it does not "
                        "locate it.")
    p.add_argument("--assembly", required=True)
    p.add_argument("--scaffold", required=True)
    p.add_argument("--telomere-windows", default="",
                   help="tidk search output for this scaffold. PREFERRED: tidk is the tool of "
                        "record and already a dependency. Without it the motif is counted "
                        "from the sequence instead, which is equivalent for CCCTAA but does "
                        "not normalise the canonical repeat.")
    p.add_argument("--telomere-motif", default="CCCTAA",
                   help="only used for the sequence-derived fallback")
    p.add_argument("--outdir", required=True)
    p.add_argument("--label", default="")
    p.add_argument("--win-bins", type=int, default=20,
                   help="half-window in bins for the cross-contact statistic (default 20, "
                        "so 2 Mb either side at 100 kb)")
    p.add_argument("--max-paf-hic-disagree", type=int, default=2000000,
                   help="if the PAF junction and the Hi-C minimum differ by more than this, "
                        "the two signals point at different places and the candidate goes to "
                        "REVIEW (default 2 Mb). Not zero: the PAF junction comes from 1 Mb "
                        "bins, so +/-1 bin is expected -- the measured disagreements were "
                        "0.4 and 1.4 Mb.")
    p.add_argument("--gap-snap-window", type=int, default=500000)
    a = p.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    lab = a.label or ("%s.%s" % (a.assembly, a.scaffold.replace("+", "_")))
    op = lambda s: os.path.join(a.outdir, lab + s)

    rows = [r for r in read_rows(a.candidates)
            if r.get("assembly") == a.assembly and r.get("scaffold") == a.scaffold]
    if not rows:
        sys.exit("ERROR: no row for %s / %s in %s" % (a.assembly, a.scaffold, a.candidates))
    if a.cut_bp:
        rows = [r for r in rows
                if str(r.get("cut_bp", "")) == str(a.cut_bp)] or rows
    row = rows[0]
    try:
        paf_bp = int(float(row.get("cut_bp") or row.get("junction_bp") or 0)) or None
    except ValueError:
        paf_bp = None

    name, seq = read_one_fasta(a.fasta)
    if name != a.scaffold:
        sys.stderr.write("[chimera_evidence] NOTE: fasta record %r, expected %r; using the "
                         "record present\n" % (name, a.scaffold))
    span = len(seq)

    prof, res, hic, err = None, 100000, None, ""
    if a.cool and os.path.isfile(a.cool):
        prof, res, err = hic_profile(a.cool, name, a.win_bins)
        hic = summarise_profile(prof, res) if prof is not None else None
        if err:
            sys.stderr.write("[chimera_evidence] Hi-C unavailable: %s\n" % err)
    else:
        err = "no cool supplied"
        sys.stderr.write("[chimera_evidence] no contact map: telomere and N-gap evidence "
                         "only. Those are independent of Hi-C, which cannot justify breaking "
                         "a join it made anyway -- it confirms, and confirmation can follow.\n")

    if a.telomere_windows and os.path.isfile(a.telomere_windows) \
            and os.path.getsize(a.telomere_windows) > 0:
        tw = read_tidk_windows(a.telomere_windows, name)
        telo_src = "tidk"
        if not tw:
            sys.stderr.write("[chimera_evidence] WARNING: no tidk rows for %r; falling back "
                             "to a sequence-derived motif count\n" % name)
            tw = telomere_windows(seq, a.telomere_motif)
            telo_src = "sequence_fallback_no_tidk_rows"
    else:
        tw = telomere_windows(seq, a.telomere_motif)
        telo_src = "sequence_fallback"
        sys.stderr.write("[chimera_evidence] no tidk output supplied; counting %s from the "
                         "sequence\n" % a.telomere_motif)
    # the AGP cut is the position of record; the Hi-C minimum is a cross-check on it
    target = paf_bp or (hic["min_bp"] if hic else 0)
    telo = summarise_telomere(tw, target)
    gs = gaps(seq)
    sn = snap(gs, target, a.gap_snap_window)

    # ---- the cut point, and the verdict -------------------------------------------
    # the AGP position already sits in a 100 bp scaffolding gap, so snapping should be a
    # no-op or a few kb -- measured offsets on the known candidates were -9 kb, +0 kb, -0 kb,
    # the residue of gap filling upstream of the join.
    cut_bp = sn["cut_bp"] if sn else target
    notes = []
    verdict = row.get("candidate_verdict", row.get("verdict", "?"))
    if hic is None:
        notes.append("no_hic")
    elif paf_bp and abs(hic["min_bp"] - paf_bp) > a.max_paf_hic_disagree:
        # NOT a downgrade to REVIEW any more: the AGP records the join, and Hi-C contact is
        # depleted approaching a junction as well as at it -- measured, the Hi-C minimum was
        # 256 kb and 705 kb from the AGP join on the two candidates. Disagreement is worth
        # flagging, not worth overriding an exact position with an inferred one.
        notes.append("hic_min_%+d_from_agp_cut" % (hic["min_bp"] - paf_bp))
    if hic and hic["ratio"] >= 1.0:
        notes.append("hic_not_depleted=%.3f" % hic["ratio"])
    if sn:
        notes.append("snapped_to_gap=%d(%+d)" % (sn["cut_bp"], sn["offset_from_target"]))
    else:
        notes.append("no_gap_within_%d" % a.gap_snap_window)
    if telo and telo["both_orientations"] and telo["junction_over_background"] >= 3:
        notes.append("interstitial_telomere=%dx" % round(telo["junction_over_background"]))

    with open(op(".chimera_evidence.tsv"), "w") as out:
        out.write("# Independent evidence for one candidate. The concordance vote in the\n")
        out.write("#   candidates file remains the GATE -- it is the only signal independent\n")
        out.write("#   of how this scaffold was built. Hi-C made the join, so it is used here\n")
        out.write("#   to LOCATE the junction, not to justify breaking it.\n")
        out.write("# hic_ratio < 1 means contact ACROSS the junction is below the scaffold\n")
        out.write("#   median at matched distance -- the data never supported the join.\n")
        out.write("# hic_n_low_contiguous is the discriminating statistic: five consecutive\n")
        out.write("#   low windows is a boundary, one isolated low window is noise.\n")
        out.write("# telomere ABSENCE is not evidence against breaking -- a mid-arm fusion\n")
        out.write("#   leaves none. Presence in BOTH orientations indicates fused ends.\n")
        out.write("metric\tvalue\n")
        for k, v in (("assembly", a.assembly), ("scaffold", a.scaffold),
                     ("name", row.get("name", ".")), ("span_bp", span),
                     ("vote", row.get("vote", ".")),
                     ("verdict_from_candidates", row.get("verdict", ".")),
                     ("paf_junction_bp", paf_bp if paf_bp else "NA")):
            out.write("%s\t%s\n" % (k, v))
        if hic:
            out.write("hic_resolution_bp\t%d\n" % res)
            for k in ("min_bp", "min_value", "median", "ratio", "n_low_contiguous"):
                out.write("hic_%s\t%s\n" % (k, hic[k]))
            out.write("hic_low_windows_bp\t%s\n"
                      % ",".join(str(x) for x in hic["low_bp"]))
        else:
            out.write("hic_available\tno\nhic_reason\t%s\n" % (err or "unknown"))
        out.write("telomere_source\t%s\n" % telo_src)
        if telo:
            for k, v in sorted(telo.items()):
                out.write("telomere_%s\t%s\n" % (k, v))
        out.write("n_gaps_total\t%d\n" % len(gs))
        if sn:
            for k in ("gap_start", "gap_end", "cut_bp", "offset_from_target"):
                out.write("gap_%s\t%s\n" % (k, sn[k]))
        out.write("cut_bp\t%d\n" % cut_bp)
        out.write("verdict\t%s\n" % verdict)
        out.write("notes\t%s\n" % (";".join(notes) or "."))

    if prof is not None:
        figure(op(".chimera_evidence.png"), a.cool, name, res, prof, paf_bp,
               hic["min_bp"] if hic else None, cut_bp, tw,
               "%s  %s\n%s  vote %s" % (a.assembly, row.get("name", a.scaffold),
                                        a.scaffold, row.get("vote", "?")))

    sys.stderr.write("[chimera_evidence] %s %s: cut %d, verdict %s%s\n"
                     % (a.assembly, a.scaffold, cut_bp, verdict,
                        ("  [" + ";".join(notes) + "]") if notes else ""))


if __name__ == "__main__":
    main()
