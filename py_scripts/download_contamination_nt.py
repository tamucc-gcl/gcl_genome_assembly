#!/usr/bin/env python3
"""Download a release-bound NCBI BLAST database with the official BLAST+ updater."""

import argparse
import json
from pathlib import Path
import re
import subprocess


CATALOG = "https://ftp.ncbi.nlm.nih.gov/blast/db/blastdb-metadata-1-1.json"


def catalog_entry(collection):
    data = json.loads(subprocess.check_output(
        ["curl", "--fail", "--silent", "--show-error", "--location", "--retry", "3",
         "--max-time", "60", CATALOG], text=True))
    matches = [entry for entry in data if entry["dbname"] == collection and entry["dbtype"] == "Nucleotide"]
    if len(matches) != 1:
        raise ValueError(f"Expected one nucleotide database named {collection} in the NCBI catalog")
    return matches[0]


def download(args):
    snapshot = catalog_entry(args.collection)
    if snapshot["last-updated"][:10] != args.release:
        raise ValueError(f"Requested {args.collection} snapshot {args.release}; "
                         f"NCBI currently advertises {snapshot['last-updated']}. "
                         "Use an already-stored snapshot or explicitly select the advertised date.")
    args.out.mkdir()
    subprocess.run(["update_blastdb.pl", "--source", "ncbi", "--decompress", args.collection],
                   cwd=args.out, check=True)
    if catalog_entry(args.collection) != snapshot:
        raise ValueError("NCBI changed the database snapshot during download")
    info = subprocess.check_output(["blastdbcmd", "-db", args.collection, "-dbtype", "nucl", "-info"],
                                   cwd=args.out, text=True)
    if not re.search(r"(?im)^BLASTDB Version:\s*5\s*$", info):
        raise ValueError("Expected an NCBI BLAST v5 database")
    totals = re.search(r"([\d,]+) sequences;\s+([\d,]+) total bases", info)
    expected = (snapshot["number-of-sequences"], snapshot["number-of-letters"])
    if not totals or tuple(int(value.replace(",", "")) for value in totals.groups()) != expected:
        raise ValueError("Downloaded sequence and base totals differ from the pinned NCBI catalog")
    # blastdbcmd opens the alias and validates that its listed volumes can be opened.
    version = subprocess.check_output(["blastdbcmd", "-version"], text=True).splitlines()[0]
    files = {file.name: file.stat().st_size for file in sorted(args.out.iterdir()) if file.is_file()}
    if not files:
        raise ValueError("The downloader produced no database files")
    (args.out / "database.json").write_text(json.dumps(dict(
        collection=args.collection, release=args.release, catalog=CATALOG, snapshot=snapshot,
        blast_version=version, blastdb_info=info, files=files,
    ), indent=2) + "\n")


def verify(args):
    metadata = json.loads((args.out / "database.json").read_text())
    if (metadata["collection"], metadata["release"]) != (args.collection, args.release):
        raise ValueError(f"Database identity does not match {args.collection}/{args.release}")
    for filename, size in metadata["files"].items():
        path = args.out / filename
        if not path.is_file() or path.stat().st_size != size:
            raise ValueError(f"Missing or changed file in stored BLAST database: {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["download", "verify"])
    parser.add_argument("--collection", required=True)
    parser.add_argument("--release", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9_]+", args.collection) or not re.fullmatch(r"\d{4}-\d{2}-\d{2}", args.release):
        parser.error("Use an NCBI collection name and a release date in YYYY-MM-DD format")
    (download if args.command == "download" else verify)(args)


if __name__ == "__main__":
    main()
