#!/usr/bin/env python3
"""
Scan scaffolds for telomere motif repeats near sequence ends.
Outputs per-scaffold telomere presence for 5' and 3' ends.
"""
import argparse
import gzip
import re
import sys


COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def complement(seq):
    return seq.translate(COMPLEMENT)


def reverse_complement(seq):
    return complement(seq)[::-1]


def generate_motif_variants(motif):
    """Generate all 4 variants: original, complement, reverse, reverse complement."""
    motif = motif.upper()
    variants = {
        motif,
        complement(motif),
        motif[::-1],
        reverse_complement(motif),
    }
    return sorted(variants)


def open_fasta(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_fasta_ends(path, window):
    """Yield (name, length, head_seq, tail_seq) for each sequence."""
    name = None
    head = ""
    tail = ""
    length = 0

    with open_fasta(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, length, head, tail
                name = line[1:].split()[0]
                head = ""
                tail = ""
                length = 0
                continue

            seq = line.upper()
            if len(head) < window:
                need = window - len(head)
                head += seq[:need]
            length += len(seq)
            if window > 0:
                if len(seq) >= window:
                    tail = seq[-window:]
                else:
                    tail = (tail + seq)[-window:]

        if name is not None:
            yield name, length, head, tail


def build_pattern(motifs, min_repeats):
    """Build regex pattern matching any motif repeated min_repeats+ times."""
    parts = [f"(?:{m}){{{min_repeats},}}" for m in motifs]
    return re.compile("|".join(parts), re.IGNORECASE)


def scan_scaffold(head, tail, pattern):
    """Return (has_5prime, has_3prime) booleans."""
    has_5 = bool(pattern.search(head))
    has_3 = bool(pattern.search(tail))
    return has_5, has_3


def main():
    parser = argparse.ArgumentParser(
        description="Scan scaffolds for telomere motifs at sequence ends."
    )
    parser.add_argument("fasta", help="Input FASTA file (gzip supported)")
    parser.add_argument(
        "--id",
        dest="haplotype_id",
        required=True,
        help="Haplotype/sample identifier for output",
    )
    parser.add_argument(
        "--motif",
        default="TTAGGG",
        help="Telomere motif (default: TTAGGG). Complement/reverse variants generated automatically.",
    )
    parser.add_argument(
        "--window",
        type=int,
        default=10000,
        help="Window size (bp) at each end to search (default: 10000)",
    )
    parser.add_argument(
        "--min-repeats",
        type=int,
        default=10,
        help="Minimum consecutive motif repeats required (default: 10)",
    )
    parser.add_argument(
        "--output",
        "-o",
        default="-",
        help="Output TSV file (default: stdout)",
    )
    args = parser.parse_args()

    motifs = generate_motif_variants(args.motif)
    pattern = build_pattern(motifs, args.min_repeats)

    out_fh = sys.stdout if args.output == "-" else open(args.output, "w")

    header = [
        "haplotype_id",
        "scaffold",
        "length",
        "telomere_5prime",
        "telomere_3prime",
        "telomere_both",
        "window",
        "min_repeats",
        "motifs",
    ]
    print("\t".join(header), file=out_fh)

    for scaffold_name, length, head, tail in parse_fasta_ends(args.fasta, args.window):
        has_5, has_3 = scan_scaffold(head, tail, pattern)
        row = [
            args.haplotype_id,
            scaffold_name,
            str(length),
            "1" if has_5 else "0",
            "1" if has_3 else "0",
            "1" if (has_5 and has_3) else "0",
            str(args.window),
            str(args.min_repeats),
            ",".join(motifs),
        ]
        print("\t".join(row), file=out_fh)

    if out_fh is not sys.stdout:
        out_fh.close()


if __name__ == "__main__":
    main()