#!/usr/bin/env python3
"""
mitos_to_genbank.py
===================
Convert a MITOS2 annotation (result.bed) + the mitogenome FASTA into a GenBank file
whose structure matches what the pipeline's circular-map (plot_mito_circular.py) and the
MitoHiFi-style stats expect: one record with the sequence, CDS/tRNA/rRNA features each
carrying /gene and /product qualifiers, and a circular/linear topology.

This is what lets a GetOrganelle+MITOS2 animal mitogenome slot into the same downstream
GenBank consumers as a MitoHiFi mitogenome.

Usage:
    mitos_to_genbank.py --fasta mito.fasta --bed result.bed \
        --sample_id SAMPLE --organelle animal_mt --genetic_code 2 [--circular] \
        --output SAMPLE_animal_mt_mitogenome.gb
"""
import argparse
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

# single-letter -> 3-letter, so tRNA products read "tRNA-Phe" (plot_mito_circular shortens back)
AA1TO3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln",
    "E": "Glu", "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
    "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
    "Y": "Tyr", "V": "Val",
}


def classify(name):
    """MITOS feature name -> (feature_type, gene, product)."""
    n = name.strip()
    low = n.lower()
    if low.startswith("trn"):
        rest = n[3:]
        aa = rest[0].upper() if rest else ""
        return "tRNA", n, "tRNA-" + AA1TO3.get(aa, "Xxx")
    if low.startswith("rrn"):
        return "rRNA", n, n
    return "CDS", n, n


def parse_args():
    p = argparse.ArgumentParser(description="MITOS2 BED + FASTA -> GenBank")
    p.add_argument("--fasta", required=True)
    p.add_argument("--bed", required=True)
    p.add_argument("--sample_id", required=True)
    p.add_argument("--organelle", default="mitogenome")
    p.add_argument("--genetic_code", type=int, default=2)
    p.add_argument("--circular", action="store_true")
    p.add_argument("--output", required=True)
    return p.parse_args()


def main():
    a = parse_args()

    # MITOS annotates a single sequence; use the longest record in the FASTA.
    records = list(SeqIO.parse(a.fasta, "fasta"))
    if not records:
        sys.exit(f"ERROR: no sequence in {a.fasta}")
    src = max(records, key=lambda r: len(r.seq))
    seq = Seq(str(src.seq).upper())
    L = len(seq)

    rec = SeqRecord(
        seq,
        id=(a.sample_id[:16] or "mito"),
        name="mitogenome",
        description=f"{a.sample_id} {a.organelle} mitogenome (MITOS2 annotation)",
    )
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"] = "circular" if a.circular else "linear"
    rec.annotations["organism"] = a.sample_id
    rec.features.append(SeqFeature(
        FeatureLocation(0, L), type="source",
        qualifiers={"organism": [a.sample_id], "mol_type": ["genomic DNA"]},
    ))

    n_cds = n_trna = n_rrna = 0
    with open(a.bed) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            if len(f) < 6:
                f = line.split()
            if len(f) < 6:
                continue
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError:
                continue
            name = f[3]
            strand = -1 if f[5] == "-" else 1
            start, end = max(0, start), min(L, end)
            if end <= start:
                continue
            ftype, gene, prod = classify(name)
            rec.features.append(SeqFeature(
                FeatureLocation(start, end, strand=strand),
                type=ftype, qualifiers={"gene": [gene], "product": [prod]},
            ))
            n_cds += ftype == "CDS"
            n_trna += ftype == "tRNA"
            n_rrna += ftype == "rRNA"

    SeqIO.write(rec, a.output, "genbank")
    sys.stderr.write(
        f"[mitos_to_genbank] {a.output}: len={L} CDS={n_cds} tRNA={n_trna} rRNA={n_rrna} "
        f"topology={rec.annotations['topology']}\n"
    )


if __name__ == "__main__":
    main()
