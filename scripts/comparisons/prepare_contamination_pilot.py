#!/usr/bin/env python3
"""Write deterministic known-answer inputs for the contamination evidence pilot."""

import argparse
import json
from pathlib import Path
import random


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True)
    randomizer = random.Random(524803)
    codons = dict(zip("ACDEFGHIKLMNPQRSTVWY",
                      "GCT TGT GAT GAA TTT GGT CAT ATT AAA CTG ATG AAT CCT CAA CGT TCT ACT GTT TGG TAT".split()))
    host, nonhost = ("".join(randomizer.choices(list(codons), k=200)) for _ in range(2))
    host_dna, foreign_dna = ("".join(codons[aa] for aa in sequence) for sequence in (host, nonhost))
    reverse_foreign = foreign_dna.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    # Stop codons around the coding segments prevent fortuitous extension.
    pad = "TAA" * 40
    reverse_pad = "TTA" * 40
    contigs = {
        "host": pad + host_dna + pad,
        "foreign": reverse_pad + reverse_foreign + reverse_pad,
        "mixed": pad + host_dna + pad + reverse_pad + reverse_foreign + reverse_pad,
        "no_match": "".join(randomizer.choices("ACGT", k=1500)),
        "ambiguous": "N" * 200,
    }
    (args.out / "proteins.faa").write_text(
        f">host_reference OX=7955\n{host}\n>foreign_reference OX=83333\n{nonhost}\n", encoding="ascii")
    (args.out / "nucleotide.fasta").write_text(
        f">host_reference\n{host_dna}\n>foreign_reference\n{foreign_dna}\n", encoding="ascii")
    (args.out / "nucleotide_taxids.tsv").write_text("host_reference\t7955\nforeign_reference\t83333\n")
    (args.out / "assembly.fasta").write_text(
        "".join(f">{name}\n{seq}\n" for name, seq in contigs.items()), encoding="ascii")
    with (args.out / "reads.sam").open("w", encoding="ascii") as sam:
        sam.write("@HD\tVN:1.6\tSO:coordinate\n")
        for name, seq in contigs.items():
            sam.write(f"@SQ\tSN:{name}\tLN:{len(seq)}\n")
        for name, seq in contigs.items():
            if name != "ambiguous":
                sam.write(f"{name}_read\t0\t{name}\t1\t60\t100M\t*\t0\t0\t{seq[:100]}\t{'I'*100}\n")
    (args.out / "expected.json").write_text(json.dumps(dict(
        sequence_evidence={"host": "host-compatible", "foreign": "non-host", "mixed": "mixed",
                           "no_match": "unresolved", "ambiguous": "unresolved"},
        non_host_intervals={"foreign": [120, 720], "mixed": [960, 1560]},
        nucleotide_target_sequences=["foreign", "mixed"]), indent=2) + "\n")


if __name__ == "__main__":
    main()
