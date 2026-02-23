#!/usr/bin/env python3
"""
Make a multi-page PDF "book" of TAD-style plots (per contig) from an .mcool.

For each contig >= --min_contig_bp:
  - Compute insulation + boundaries (cooltools)
  - Plot per-contig contact map (log1p) with boundaries overlaid
  - Plot insulation track underneath with boundary ticks

Outputs:
  - <out_prefix>.insulation.tsv
  - <out_prefix>.boundaries.tsv
  - <out_prefix>.tad_book.pdf

Requirements:
  - cooler
  - cooltools
  - pandas
  - numpy
  - matplotlib
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import cooler
import cooltools


def open_mcool_at_resolution(mcool_path: str, res: int) -> cooler.Cooler:
    candidates = [
        f"{mcool_path}::resolutions/{res}",
        f"{mcool_path}::/resolutions/{res}",
    ]
    last_err = None
    for uri in candidates:
        try:
            return cooler.Cooler(uri)
        except Exception as e:
            last_err = e
    raise RuntimeError(
        f"Could not open mcool at resolution {res}. Tried:\n  " +
        "\n  ".join(candidates) +
        f"\nLast error: {repr(last_err)}"
    )


def pick_insulation_col(df: pd.DataFrame) -> str:
    # cooltools versions vary; try common names
    for c in ["log2_insulation_score", "insulation_score"]:
        if c in df.columns:
            return c
    # sometimes window-specific
    for c in df.columns:
        if "log2_insulation" in c:
            return c
    # last resort
    for c in df.columns:
        if "insulation" in c:
            return c
    raise RuntimeError(f"No insulation score-like column found. Columns: {list(df.columns)}")


def boundary_sites(bound_df: pd.DataFrame) -> pd.DataFrame:
    # cooltools versions vary
    for c in ["is_boundary", "boundary"]:
        if c in bound_df.columns:
            return bound_df[bound_df[c] == True]
    if "boundary_strength" in bound_df.columns:
        return bound_df[~bound_df["boundary_strength"].isna()]
    return bound_df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mcool", required=True)
    ap.add_argument("--resolution", type=int, required=True)
    ap.add_argument("--window_bp", type=int, required=True)
    ap.add_argument("--assembly_id", required=True)
    ap.add_argument("--out_prefix", required=True)
    ap.add_argument("--min_contig_bp", type=int, default=5_000_000)
    ap.add_argument("--max_contigs", type=int, default=0,
                    help="0 means no cap. Otherwise plot only the largest N contigs passing min_contig_bp.")
    ap.add_argument("--balance", action="store_true")
    ap.add_argument("--ignore_diags", type=int, default=2)
    ap.add_argument("--dpi", type=int, default=200)
    args = ap.parse_args()

    clr = open_mcool_at_resolution(args.mcool, args.resolution)
    print(f"Opened: {clr.uri}")
    print(f"Binsize: {clr.binsize}")

    # balancing check
    use_balance = args.balance
    if use_balance:
        bins_head = clr.bins()[:5]
        if "weight" not in bins_head.columns:
            print("WARNING: --balance requested but no 'weight' column found; proceeding unbalanced.")
            use_balance = False

    # Compute insulation
    ins = cooltools.insulation(
        clr,
        window_bp=args.window_bp,
        ignore_diags=args.ignore_diags,
        clr_weight_name="weight" if use_balance else None
    )
    ins_path = f"{args.out_prefix}.insulation.tsv"
    ins.to_csv(ins_path, sep="\t", index=False)

    # Call boundaries
    try:
        bnd = cooltools.call_boundaries(ins)
    except AttributeError:
        from cooltools.insulation import call_boundaries
        bnd = call_boundaries(ins)
    bnd_path = f"{args.out_prefix}.boundaries.tsv"
    bnd.to_csv(bnd_path, sep="\t", index=False)

    score_col = pick_insulation_col(ins)

    # Choose contigs >= threshold
    chromsizes = clr.chromsizes.to_dict()
    contigs = [(c, chromsizes[c]) for c in clr.chromnames if chromsizes[c] >= args.min_contig_bp]
    if not contigs:
        raise RuntimeError(f"No contigs >= {args.min_contig_bp} bp found.")

    # Sort by length, optionally cap
    contigs.sort(key=lambda x: x[1], reverse=True)
    if args.max_contigs and args.max_contigs > 0:
        contigs = contigs[:args.max_contigs]

    # Pre-slice insulation/boundaries by contig for speed
    ins_by = {c: ins[ins["chrom"] == c].copy() for c, _L in contigs}
    bnd_by = {c: bnd[bnd["chrom"] == c].copy() for c, _L in contigs}

    pdf_path = f"{args.out_prefix}.tad_book.pdf"
    with PdfPages(pdf_path) as pdf:
        for contig, L in contigs:
            ins_c = ins_by[contig]
            bnd_c = bnd_by[contig]
            bsites = boundary_sites(bnd_c)

            # Fetch contig matrix
            M = clr.matrix(balance=use_balance).fetch(contig)
            M = np.asarray(M)
            M = np.log1p(M)  # stabilize dynamic range

            # Figure: contact map + insulation track
            fig = plt.figure(figsize=(8.5, 10.5))  # letter-ish portrait
            gs = fig.add_gridspec(2, 1, height_ratios=[4.2, 1.2], hspace=0.25)

            ax0 = fig.add_subplot(gs[0])
            ax0.imshow(M, origin="lower", aspect="auto")
            ax0.set_title(
                f"{args.assembly_id}\n{contig}  |  {L:,} bp  |  res={args.resolution:,}  window={args.window_bp:,}"
            )
            ax0.set_xlabel("Bin")
            ax0.set_ylabel("Bin")

            # Overlay boundary lines on heatmap
            # Convert boundary starts to bin indices
            if not bsites.empty and "start" in bsites.columns:
                b_bins = (bsites["start"].values // args.resolution).astype(int)
                b_bins = b_bins[(b_bins >= 0) & (b_bins < M.shape[0])]
                for bb in b_bins:
                    ax0.axvline(bb, linewidth=0.5)
                    ax0.axhline(bb, linewidth=0.5)

            ax1 = fig.add_subplot(gs[1])
            mids = (ins_c["start"].values + ins_c["end"].values) / 2
            ax1.plot(mids, ins_c[score_col].values, linewidth=0.8)
            ax1.axhline(0, linewidth=1)
            ax1.set_xlabel(f"Position on {contig} (bp)")
            ax1.set_ylabel(score_col)

            # Mark boundaries as vertical ticks on insulation track
            if not bsites.empty and "start" in bsites.columns:
                bx = (bsites["start"].values + bsites["end"].values) / 2
                ymin = np.nanmin(ins_c[score_col].values)
                ymax = np.nanmax(ins_c[score_col].values)
                ax1.vlines(bx, ymin=ymin, ymax=ymax, linewidth=0.5)

            fig.tight_layout()
            pdf.savefig(fig, dpi=args.dpi)
            plt.close(fig)

    print("Wrote:")
    print(f"  {ins_path}")
    print(f"  {bnd_path}")
    print(f"  {pdf_path}")


if __name__ == "__main__":
    main()