#!/usr/bin/env python3
import argparse
import math
import numpy as np
import pyBigWig
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def parse_fai(fai_path):
    contigs = []
    with open(fai_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            contig = parts[0]
            length = int(parts[1])
            contigs.append((contig, length))
    return contigs

def contig_means(bw, chrom, length, bin_size, stat="mean"):
    n_bins = int(math.ceil(length / bin_size))
    starts = np.arange(0, n_bins * bin_size, bin_size, dtype=int)
    ends = np.minimum(starts + bin_size, length)

    # pyBigWig.stats returns a list of floats/None
    vals = bw.stats(chrom, 0, length, nBins=n_bins, type=stat)
    y = np.array([v if v is not None else np.nan for v in vals], dtype=float)
    x = (starts + ends) / 2.0
    return x, y

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bw", required=True, help="Coverage bigWig (e.g. sample.cov.1kb.bw)")
    ap.add_argument("--fai", required=True, help="FASTA index (.fai) for contig order/lengths")
    ap.add_argument("--out_pdf", required=True, help="Output multi-page PDF")
    ap.add_argument("--bin_size", type=int, default=1000, help="Bin size used for plotting (default 1000)")
    ap.add_argument("--max_contigs", type=int, default=0, help="If >0, plot only first N contigs")
    ap.add_argument("--min_len", type=int, default=0, help="Skip contigs shorter than this length")
    ap.add_argument("--cap_y", type=float, default=0.0, help="If >0, cap y-axis at this value")
    ap.add_argument("--also_png_dir", default="", help="If set, also write one PNG per contig into this dir")
    args = ap.parse_args()

    contigs = parse_fai(args.fai)
    if args.max_contigs and args.max_contigs > 0:
        contigs = contigs[:args.max_contigs]
    if args.min_len and args.min_len > 0:
        contigs = [(c,l) for c,l in contigs if l >= args.min_len]

    bw = pyBigWig.open(args.bw)

    with PdfPages(args.out_pdf) as pdf:
        for chrom, length in contigs:
            if chrom not in bw.chroms():
                continue

            x, y = contig_means(bw, chrom, length, args.bin_size, stat="mean")

            fig = plt.figure(figsize=(11, 3))
            ax = fig.add_subplot(111)
            ax.plot(x, y, linewidth=0.7)
            ax.set_title(f"{chrom}  (len={length:,} bp)  bin={args.bin_size} bp")
            ax.set_xlabel("Position (bp)")
            ax.set_ylabel("Coverage")
            ax.set_xlim(0, length)

            if args.cap_y and args.cap_y > 0:
                ax.set_ylim(0, args.cap_y)
            else:
                # Robust autoscale ignoring crazy spikes
                finite = y[np.isfinite(y)]
                if finite.size > 0:
                    cap = np.nanpercentile(finite, 99.5)
                    ax.set_ylim(0, cap if cap > 0 else 1)

            ax.grid(True, linewidth=0.3, alpha=0.4)

            pdf.savefig(fig, bbox_inches="tight")
            if args.also_png_dir:
                import os
                os.makedirs(args.also_png_dir, exist_ok=True)
                fig.savefig(os.path.join(args.also_png_dir, f"{chrom}.coverage.png"),
                            dpi=200, bbox_inches="tight")
            plt.close(fig)

    bw.close()

if __name__ == "__main__":
    main()
