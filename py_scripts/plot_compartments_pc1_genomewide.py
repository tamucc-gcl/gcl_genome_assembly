#!/usr/bin/env python3
"""
Compute and plot A/B compartments (PC1/E1) from a .mcool.

Outputs:
  1) PC1 track as bedGraph-like TSV (chrom, start, end, E1)
  2) Genome-wide concatenated PC1 plot (one dot/segment per bin) with scaffold boundaries

Requires:
  - cooler
  - cooltools
  - pandas
  - numpy
  - matplotlib

Usage:
  plot_compartments_pc1_genomewide.py \
    --mcool sample.mcool \
    --resolution 250000 \
    --assembly_id Sde_CLim_110_hap1 \
    --out_prefix Sde_CLim_110_hap1.comp_250kb \
    --min_contig_bp 5000000 \
    --max_contigs 30
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import cooler
import cooltools


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mcool", required=True, help="Input .mcool")
    ap.add_argument("--resolution", type=int, required=True, help="Resolution in bp (e.g., 250000)")
    ap.add_argument("--assembly_id", required=True, help="Label for titles")
    ap.add_argument("--out_prefix", required=True, help="Prefix for outputs")
    ap.add_argument("--min_contig_bp", type=int, default=5_000_000,
                    help="Only label/plot contigs >= this length (bp). Smaller contigs still included but not labeled.")
    ap.add_argument("--max_contigs", type=int, default=30,
                    help="Max number of contigs to label on x-axis (largest first).")
    ap.add_argument("--ignore_diags", type=int, default=2, help="Ignore first N diagonals for eigenvector calc")
    return ap.parse_args()


def main():
    args = parse_args()

    clr = cooler.Cooler(f"{args.mcool}::resolutions/{args.resolution}")

    # Compute eigenvectors per chromosome/contig (cis).
    # Returns:
    #   eigvecs: bins with columns chrom,start,end,E1,E2,...
    #   eigvals: per-chrom eigenvalues
    eigvals, eigvecs = cooltools.eigs_cis(
        clr,
        n_eigs=2,
        ignore_diags=args.ignore_diags
    )

    # Keep PC1
    if "E1" not in eigvecs.columns:
        raise RuntimeError("Expected column 'E1' from cooltools.eigs_cis, but did not find it.")

    pc1 = eigvecs[["chrom", "start", "end", "E1"]].copy()
    pc1_path = f"{args.out_prefix}.pc1.bedGraph"
    pc1.to_csv(pc1_path, sep="\t", header=True, index=False)

    # Build genome-wide concatenated coordinate system
    chromsizes = clr.chromsizes.to_dict()  # length in bp
    chrom_order = list(clr.chromnames)

    # Compute offsets
    offsets = {}
    running = 0
    boundaries = []  # (xpos, chrom)
    for chrom in chrom_order:
        offsets[chrom] = running
        running += chromsizes[chrom]
        boundaries.append((running, chrom))

    # Add midpoints for bins to plot
    pc1["mid"] = (pc1["start"].values + pc1["end"].values) / 2
    pc1["gpos"] = pc1.apply(lambda r: offsets[r["chrom"]] + r["mid"], axis=1)

    # Choose contigs to label: largest >= min_contig_bp; cap at max_contigs
    chrom_df = pd.DataFrame({"chrom": chrom_order, "length": [chromsizes[c] for c in chrom_order]})
    chrom_df = chrom_df.sort_values("length", ascending=False)
    label_df = chrom_df[chrom_df["length"] >= args.min_contig_bp].head(args.max_contigs)

    # For labels, compute midpoint in genome coords
    label_positions = []
    for _, row in label_df.iterrows():
        c = row["chrom"]
        start = offsets[c]
        end = offsets[c] + chromsizes[c]
        label_positions.append(( (start + end) / 2, c ))

    # Plot
    fig = plt.figure(figsize=(14, 4.5))
    ax = fig.add_subplot(111)

    ax.plot(pc1["gpos"].values, pc1["E1"].values, linewidth=0.6)
    ax.axhline(0, linewidth=1)

    # Scaffold boundaries (thin verticals)
    for xpos, _chrom in boundaries:
        ax.axvline(xpos, linewidth=0.4)

    # X-axis labeling: only the biggest contigs
    if label_positions:
        xticks = [p for p, c in label_positions]
        xlabels = [c for p, c in label_positions]
        ax.set_xticks(xticks)
        ax.set_xticklabels(xlabels, rotation=90)
    else:
        ax.set_xticks([])

    ax.set_ylabel("Compartment score (PC1 / E1)")
    ax.set_xlabel("Genome coordinate (concatenated scaffolds)")
    ax.set_title(f"{args.assembly_id} A/B compartments (PC1) @ {args.resolution} bp")

    # Make the plot tighter
    ax.margins(x=0)
    fig.tight_layout()

    out_png = f"{args.out_prefix}.pc1.genomewide.png"
    # out_pdf = f"{args.out_prefix}.pc1.genomewide.pdf"
    fig.savefig(out_png, dpi=300)
    # fig.savefig(out_pdf)
    plt.close(fig)

    # Also write the eigenvalues table (helpful QC)
    eigvals_path = f"{args.out_prefix}.eigvals.tsv"
    eigvals.to_csv(eigvals_path, sep="\t", index=False)

    print(f"Wrote:\n  {pc1_path}\n  {out_png}\n  {eigvals_path}")


if __name__ == "__main__":
    main()