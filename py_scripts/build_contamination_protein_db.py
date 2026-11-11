#!/usr/bin/env python3
"""Build a taxonomy-aware DIAMOND database from canonical UniProt proteomes."""

import argparse
from contextlib import contextmanager
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess
import tarfile
import tempfile


CANONICAL_MEMBER = re.compile(r"UP\d+_\d+\.fasta(?:\.gz)?")
TAXID = re.compile(r"(?:^|\s)OX=(\d+)(?:\s|$)")


def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


@contextmanager
def fasta_text(stream):
    """Use magic bytes because task-local input paths may lack extensions."""
    with io.BufferedReader(stream) as buffered:
        raw = gzip.GzipFile(fileobj=buffered) if buffered.peek(2)[:2] == b"\x1f\x8b" else buffered
        with io.TextIOWrapper(raw, encoding="utf-8") as text:
            yield text


def protein_records(source, source_format):
    def records(handle, label):
        header = None
        pieces = []
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                if header is not None:
                    yield label, header, "".join(pieces)
                header, pieces = line[1:], []
            elif line:
                if header is None:
                    raise ValueError(f"Sequence before first FASTA header in {label}")
                pieces.append(line)
        if header is not None:
            yield label, header, "".join(pieces)

    if source_format == "fasta":
        with fasta_text(source.open("rb")) as handle:
            yield from records(handle, source.name)
    else:
        # Read member streams directly. No archive paths are extracted to disk.
        with tarfile.open(source, mode="r|*") as archive:
            for member in archive:
                if member.isfile() and CANONICAL_MEMBER.fullmatch(Path(member.name).name):
                    with fasta_text(archive.extractfile(member)) as handle:
                        yield from records(handle, member.name)


def prepare(source, source_format, nodes, merged, fasta, mapping, probe):
    with nodes.open() as handle:
        known = {int(line.split("|", 1)[0]) for line in handle if line.strip()}
    aliases = {}
    if merged:
        with merged.open() as handle:
            for line in handle:
                old, new, *_ = line.split("|")
                aliases[int(old)] = int(new)
    count, residues, reassigned = 0, 0, 0
    probe_taxid = None
    with gzip.open(fasta, "wt", encoding="utf-8") as fa, gzip.open(mapping, "wt", encoding="ascii") as tax:
        tax.write("accession\taccession.version\ttaxid\tgi\n")
        for label, header, sequence in protein_records(source, source_format):
            match = TAXID.search(header)
            if match is None or not sequence:
                raise ValueError(f"{label}: missing UniProt OX taxid or sequence for {header!r}")
            old = int(match[1])
            taxid = aliases.get(old, old)
            if taxid not in known or taxid <= 1:
                raise ValueError(f"{label}: taxid {old} absent from the supplied taxonomy or merged IDs")
            count += 1
            residues += len(sequence)
            reassigned += taxid != old
            # Give each archive entry its own ID even if two proteomes share an accession.
            subject = f"UPR{count}"
            fa.write(f">{subject} {header}\n{sequence}\n")
            tax.write(f"{subject}\t{subject}\t{taxid}\t0\n")
            if probe_taxid is None and len(sequence) >= 80:
                probe_taxid = taxid
                probe.write_text(f">{subject}\n{sequence}\n", encoding="ascii")
    if not count or probe_taxid is None:
        raise ValueError("Need canonical protein entries including one of at least 80 aa for the taxonomy probe")
    return count, residues, reassigned, probe_taxid


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, help="Local file or HTTPS URL")
    parser.add_argument("--source-sha256", required=True)
    parser.add_argument("--format", choices=["fasta", "tar"], required=True)
    parser.add_argument("--release", required=True)
    parser.add_argument("--taxonomy", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    if args.threads < 1 or not re.fullmatch(r"[0-9a-fA-F]{64}", args.source_sha256):
        parser.error("Provide a SHA256 digest and positive thread count")

    nodes, names, merged = (args.taxonomy / name for name in ("nodes.dmp", "names.dmp", "merged.dmp"))
    taxonomy = [nodes, names] + ([merged] if merged.is_file() else [])
    before = {path.name: digest(path) for path in taxonomy}
    with tempfile.TemporaryDirectory(prefix="protein-inputs-", dir=".") as tmp:
        work = Path(tmp)
        if args.source.startswith("https://"):
            source = work / "download"
            subprocess.run(["curl", "--fail", "--location", "--retry", "3",
                            "--output", str(source), args.source], check=True)
        else:
            source = Path(args.source)
        source_sha = digest(source)
        if source_sha != args.source_sha256.lower():
            raise ValueError(f"Source SHA256 mismatch: observed {source_sha}")

        fasta, mapping, probe = (work / name for name in ("proteins.faa.gz", "taxonmap.tsv.gz", "probe.faa"))
        count, residues, reassigned, taxid = prepare(
            source, args.format, nodes, merged if merged.is_file() else None, fasta, mapping, probe)
        db = work / "reference_proteomes.dmnd"
        subprocess.run(["diamond", "makedb", "--in", str(fasta), "--db", str(db),
                        "--no-parse-seqids", "--taxonmap", str(mapping), "--taxonnodes", str(nodes),
                        "--taxonnames", str(names), "--threads", str(args.threads)], check=True)
        hits = subprocess.check_output(["diamond", "blastp", "--db", str(db), "--query", str(probe),
                                        "--outfmt", "6", "qseqid", "sseqid", "staxids",
                                        "--max-target-seqs", "0", "--threads", str(args.threads)], text=True)
        if not any(query == subject and str(taxid) in taxids.split(";")
                   for query, subject, taxids in (line.split("\t") for line in hits.splitlines())):
            raise ValueError("DIAMOND self-search did not recover the expected taxid")
        if before != {path.name: digest(path) for path in taxonomy}:
            raise ValueError("Taxonomy files changed during database build")
        version = subprocess.check_output(["diamond", "version"], text=True).strip()
        manifest = dict(collection="uniprot_reference_proteomes", release=args.release, source=args.source,
                        source_format=args.format, source_sha256=source_sha, taxonomy_sha256=before,
                        diamond_version=version, protein_records=count, amino_acids=residues,
                        merged_taxid_records=reassigned, probe_taxid=taxid, database_sha256=digest(db))
        db.rename("reference_proteomes.dmnd")
        Path("protein_database.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
