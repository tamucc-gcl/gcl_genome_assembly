#!/usr/bin/env python3
"""
plot_mito_circular.py
=====================
Generates a circular gene map of an annotated mitochondrial genome
from a GenBank file using pyCirclize.

Usage:
    python plot_mito_circular.py \
        --genbank sample_mitogenome.gb \
        --sample_id sample_name \
        --output sample_mito_circular.png \
        [--dpi 300] [--width 10]
"""

import argparse
import sys
from pathlib import Path

from pycirclize import Circos
from pycirclize.parser import Genbank
from Bio import SeqIO
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt


# ── Colour palette ──────────────────────────────────────────────────────────
FEATURE_COLOURS = {
    "CDS":    "#E64B35",   # red-orange
    "tRNA":   "#4DBBD5",   # teal
    "rRNA":   "#F9C74F",   # gold
}


def parse_args():
    p = argparse.ArgumentParser(
        description="Circular mitogenome map from GenBank annotation")
    p.add_argument("--genbank", required=True,
                   help="Input GenBank (.gb) file")
    p.add_argument("--sample_id", required=True,
                   help="Sample identifier (used in title)")
    p.add_argument("--output", required=True,
                   help="Output PNG path")
    p.add_argument("--dpi", type=int, default=300,
                   help="Figure resolution [300]")
    p.add_argument("--width", type=float, default=10,
                   help="Figure width/height in inches [10]")
    return p.parse_args()


def get_label(feature, ftype):
    """Extract a short display label from a GenBank feature."""
    if ftype == "CDS":
        for key in ("gene", "product"):
            vals = feature.qualifiers.get(key)
            if vals:
                return vals[0]
    elif ftype == "tRNA":
        for key in ("product", "gene"):
            vals = feature.qualifiers.get(key)
            if vals:
                label = vals[0]
                # Shorten: "tRNA-Phe" → "F"
                if "tRNA-" in label:
                    aa = label.split("tRNA-")[-1]
                    aa_map = {
                        "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D",
                        "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
                        "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
                        "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
                        "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
                    }
                    return aa_map.get(aa, aa)
                return label
    elif ftype == "rRNA":
        for key in ("product", "gene"):
            vals = feature.qualifiers.get(key)
            if vals:
                return vals[0]
    return ""


def add_feature_labels(track, features, ftype, fontsize=6):
    """Add text labels at the midpoint of each feature on the track."""
    for feat in features:
        start = int(feat.location.start)
        end = int(feat.location.end)
        mid = (start + end) / 2
        label = get_label(feat, ftype)
        if label:
            track.text(label, mid, fontsize=fontsize, orientation="curved")


def main():
    args = parse_args()

    gb_path = Path(args.genbank)
    if not gb_path.exists():
        sys.exit(f"ERROR: GenBank file not found: {gb_path}")

    # ── Parse GenBank ───────────────────────────────────────────────────
    gbk = Genbank(gb_path)
    genome_size = gbk.range_size

    # ── Set up Circos ───────────────────────────────────────────────────
    circos = Circos(
        sectors={gbk.name: genome_size},
        space=0,
    )
    sector = circos.sectors[0]

    # ── Track 1 (outermost): CDS features ──────────────────────────────
    cds_track = sector.add_track((88, 100))
    cds_track.axis(fc="none", ec="grey", lw=0.3)

    cds_features = gbk.extract_features("CDS")
    if cds_features:
        cds_track.genomic_features(
            cds_features,
            plotstyle="arrow",
            fc=FEATURE_COLOURS["CDS"],
            ec="black",
            lw=0.3,
        )
        add_feature_labels(cds_track, cds_features, "CDS", fontsize=6)

    # ── Track 2: tRNA features ──────────────────────────────────────────
    trna_track = sector.add_track((78, 87))
    trna_track.axis(fc="none", ec="grey", lw=0.3)

    trna_features = gbk.extract_features("tRNA")
    if trna_features:
        trna_track.genomic_features(
            trna_features,
            plotstyle="arrow",
            fc=FEATURE_COLOURS["tRNA"],
            ec="black",
            lw=0.3,
        )
        add_feature_labels(trna_track, trna_features, "tRNA", fontsize=5)

    # ── Track 3: rRNA features ──────────────────────────────────────────
    rrna_track = sector.add_track((68, 77))
    rrna_track.axis(fc="none", ec="grey", lw=0.3)

    rrna_features = gbk.extract_features("rRNA")
    if rrna_features:
        rrna_track.genomic_features(
            rrna_features,
            plotstyle="arrow",
            fc=FEATURE_COLOURS["rRNA"],
            ec="black",
            lw=0.3,
        )
        add_feature_labels(rrna_track, rrna_features, "rRNA", fontsize=6)

    # ── Track 4 (innermost): GC content ────────────────────────────────
    gc_track = sector.add_track((45, 65))
    gc_track.axis(fc="none", ec="grey", lw=0.3)

    # Compute GC in sliding windows
    record = SeqIO.read(gb_path, "genbank")
    seq = str(record.seq).upper()
    window = max(100, genome_size // 200)
    step = window // 2

    positions, gc_values = [], []
    for i in range(0, len(seq) - window + 1, step):
        chunk = seq[i:i + window]
        gc = (chunk.count("G") + chunk.count("C")) / len(chunk) * 100
        positions.append(i + window // 2)
        gc_values.append(gc)

    if gc_values:
        mean_gc = sum(gc_values) / len(gc_values)
        gc_track.fill_between(
            x=positions,
            y1=gc_values,
            y2=mean_gc,
            color="grey",
            alpha=0.5,
            ec="dimgrey",
            lw=0.3,
        )

    # ── Centre text ─────────────────────────────────────────────────────
    size_kb = genome_size / 1000
    circos.text(
        f"{args.sample_id}\nmitogenome\n{size_kb:,.1f} kb",
        size=10,
        r=15,
    )

    # ── Legend ───────────────────────────────────────────────────────────
    legend_handles = [
        mpatches.Patch(color=FEATURE_COLOURS["CDS"],  label="CDS"),
        mpatches.Patch(color=FEATURE_COLOURS["tRNA"], label="tRNA"),
        mpatches.Patch(color=FEATURE_COLOURS["rRNA"], label="rRNA"),
        mpatches.Patch(color="grey", alpha=0.5,       label="GC content"),
    ]

    fig = circos.plotfig(figsize=(args.width, args.width))
    _ = fig.legend(
        handles=legend_handles,
        loc="lower right",
        fontsize=8,
        frameon=True,
        framealpha=0.9,
    )

    fig.savefig(args.output, dpi=args.dpi, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)

    print(f"[MITO_CIRCULAR_MAP] Saved: {args.output}")
    print(f"  Genome size: {genome_size:,} bp")
    print(f"  CDS: {len(cds_features) if cds_features else 0}")
    print(f"  tRNA: {len(trna_features) if trna_features else 0}")
    print(f"  rRNA: {len(rrna_features) if rrna_features else 0}")


if __name__ == "__main__":
    main()