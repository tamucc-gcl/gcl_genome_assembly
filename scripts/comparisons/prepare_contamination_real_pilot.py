#!/usr/bin/env python3
"""Prepare bounded real-contig and UniProt inputs for the contamination pilot."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from urllib.parse import urlencode

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "py_scripts"))
from contamination_review import digest, records, write_record, write_table


PROTEOMES = {
    "UP000515152": "Atlantic herring",
    "UP000000437": "zebrafish",
    "UP000001038": "Japanese medaka",
    "UP000000803": "fruit fly",
    "UP000002311": "baker's yeast",
    "UP000002530": "Aspergillus fumigatus",
    "UP000000625": "Escherichia coli K-12",
    "UP000002493": "Vibrio parahaemolyticus",
}


def traced_task(trace, work, name):
    with trace.open(newline="") as handle:
        matches = {row["hash"] for row in csv.DictReader(handle, delimiter="\t")
                   if row["name"] == name and row["status"] in {"COMPLETED", "CACHED"}
                   and row["exit"] == "0"}
    if len(matches) != 1:
        raise ValueError(f"Expected one successful task for {name!r} in {trace}, found {len(matches)}")
    prefix, short = matches.pop().split("/")
    tasks = [path for path in (work / prefix).glob(short + "*") if path.is_dir()]
    if len(tasks) != 1:
        raise ValueError(f"Expected one retained work directory for {name!r}, found {tasks}")
    return tasks[0]


def choose_contigs(metrics, budget):
    """Select whole contigs across size and GC strata, within a total bp cap."""
    selected = {}
    used = 0

    def take(candidates, label, fraction, limit):
        nonlocal used
        allowance = min(int(budget * fraction), budget - used)
        count = 0
        for name, length, _ in candidates:
            if name in selected or length > allowance:
                continue
            selected[name] = label
            used += length
            allowance -= length
            count += 1
            if count == limit:
                break

    usable = [(name, length, gc) for name, (length, gc) in metrics.items() if gc is not None]
    take(sorted((row for row in usable if row[1] >= 1_000_000),
                key=lambda row: (-row[1], row[0])), "long", .55, 4)
    gc_pool = [row for row in usable if 1_000 <= row[1] <= 500_000]
    take(sorted(gc_pool, key=lambda row: (row[2], row[0])), "low-GC", .12, 12)
    take(sorted(gc_pool, key=lambda row: (-row[2], row[0])), "high-GC", .12, 12)
    short = sorted((row for row in usable if row[1] <= 100_000),
                   key=lambda row: (row[1], row[0]))
    if len(short) > 48:
        short = [short[round(i * (len(short) - 1) / 47)] for i in range(48)]
    take(short, "short", .15, 24)
    mid = sorted(gc_pool, key=lambda row: (row[2], row[0]))
    if mid:
        median_gc = mid[len(mid) // 2][2]
        take(sorted(usable, key=lambda row: (abs(row[2] - median_gc), -row[1], row[0])),
             "typical-GC", 1, 12)
    if not selected:
        take(sorted(usable, key=lambda row: (row[1], row[0])), "size-limited", 1, 24)
    if not selected:
        raise ValueError("No whole contig with measurable GC fits the pilot budget")
    return selected


def subset_coverage(source, source_sizes, selected, output):
    with source.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields, rows = reader.fieldnames, list(reader)
    if not fields or not {"#rname", "startpos", "endpos", "meandepth"} <= set(fields):
        raise ValueError(f"Expected a samtools coverage table: {source}")
    lengths = {}
    for row in rows:
        name = row["#rname"]
        if name in lengths or int(row["startpos"]) != 1:
            raise ValueError(f"Duplicate or partial coverage row: {name}")
        lengths[name] = int(row["endpos"])
    if lengths != source_sizes:
        raise ValueError("Retained coverage does not match the full post-FCS FASTA dictionary")
    write_table(output, (row for row in rows if row["#rname"] in selected), fields)


def prepare_contigs(args, out):
    trace = args.project / "genome_assembly/pipeline/pipeline_trace.txt"
    work = args.project / "work"
    if args.source_fasta:
        source = args.source_fasta.resolve(strict=True)
    else:
        task = traced_task(trace, work,
                           f"CONTIG_REFINEMENT:DECONTAMINATE_ASSEMBLY_CONTIG:FCS_CLEAN_GENOME ({args.assembly})")
        source = (task / f"{args.assembly}.decontaminated.fasta").resolve(strict=True)
    coverage = args.coverage_tsv
    if args.use_trace_coverage:
        task = traced_task(trace, work, f"QC_PHASE:ASSEMBLY_QC_CONTIG_DECONTAM:MAPPING_QC ({args.assembly})")
        coverage = task / f"{args.assembly}_mapping_stats" / "coverage.txt"
    if coverage:
        coverage = coverage.resolve(strict=True)

    metrics = {}
    for header, sequence in records(source):
        upper = sequence.upper()
        counts = {base: upper.count(base) for base in "ACGT"}
        acgt = sum(counts.values())
        metrics[header.split()[0]] = (len(sequence), (counts["G"] + counts["C"]) / acgt if acgt else None)
    selected = choose_contigs(metrics, args.max_bp)
    subset = out / "contigs.fasta"
    manifest = []
    with subset.open("w", encoding="ascii") as fasta:
        for header, sequence in records(source):
            name = header.split()[0]
            if name not in selected:
                continue
            length, gc = metrics[name]
            if len(sequence) != length:
                raise ValueError(f"Source contig length changed during selection: {name}")
            write_record(fasta, header, sequence)
            manifest.append(dict(sequence=name, length_bp=length, gc_fraction=f"{gc:.6f}",
                                 selection=selected[name]))
    write_table(out / "selected.tsv", manifest, ["sequence", "length_bp", "gc_fraction", "selection"])
    if coverage:
        subset_coverage(coverage, {name: length for name, (length, _) in metrics.items()},
                        selected, out / "coverage.tsv")
    metadata = dict(assembly=args.assembly, stage="post_fcs", source=str(source), source_sha256=digest(source),
                    source_contigs=len(metrics), source_bp=sum(length for length, _ in metrics.values()),
                    selected_contigs=len(manifest), selected_bp=sum(row["length_bp"] for row in manifest),
                    maximum_bp=args.max_bp, subset_sha256=digest(subset),
                    coverage_source=str(coverage) if coverage else None,
                    coverage_method="production samtools coverage defaults" if args.use_trace_coverage else None,
                    production_trace=str(trace) if not args.source_fasta else None)
    (out / "selection.json").write_text(json.dumps(metadata, indent=2) + "\n")


def fetch_proteome(job):
    accession, label, directory = job
    url = "https://rest.uniprot.org/uniprotkb/stream?" + urlencode(
        dict(query=f"proteome:{accession}", format="fasta", compressed="true"))
    fasta = directory / f"{accession}.faa.gz"
    headers = directory / f"{accession}.headers"
    subprocess.run(["curl", "--fail", "--silent", "--show-error", "--location",
                    "--retry", "3", "--max-time", "180", "--dump-header", str(headers),
                    "--output", str(fasta), url], check=True)
    release = re.findall(r"(?im)^x-uniprot-release:\s*(\S+)", headers.read_text())
    count = 0
    with gzip.open(fasta, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if not re.search(r"\bOX=\d+\b", line):
                    raise ValueError(f"Missing UniProt taxid in {accession}")
                count += 1
    if count < 1000 or not release:
        raise ValueError(f"Incomplete UniProt response for {accession}: {count} proteins, release {release}")
    return fasta, dict(proteome=accession, organism=label, release=release[-1], proteins=count,
                       source=url, sha256=digest(fasta))


def prepare_panel(args, out):
    with ThreadPoolExecutor(max_workers=4) as workers:
        proteomes = list(workers.map(fetch_proteome,
                         ((accession, label, out) for accession, label in PROTEOMES.items())))
    provenance = [row for _, row in proteomes]
    if len({row["release"] for row in provenance}) != 1:
        raise ValueError("UniProt release changed during download. Retry for one coherent snapshot.")
    # One gzip stream, accepted by DIAMOND and by tools that read only the first gzip member.
    panel = out / "proteins.faa.gz"
    with gzip.open(panel, "wb", compresslevel=3) as output:
        for fasta, _ in proteomes:
            with gzip.open(fasta, "rb") as handle:
                shutil.copyfileobj(handle, output)
            fasta.unlink()
    write_table(out / "proteomes.tsv", provenance, list(provenance[0]))
    (out / "panel.json").write_text(json.dumps(dict(
        release=provenance[0]["release"], proteins=sum(row["proteins"] for row in provenance),
        panel_sha256=digest(panel)), indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    contigs = sub.add_parser("contigs")
    contigs.add_argument("--project", type=Path, default=Path.cwd())
    contigs.add_argument("--assembly", required=True)
    contigs.add_argument("--max-bp", type=int, default=30_000_000)
    contigs.add_argument("--source-fasta", type=Path)
    coverage = contigs.add_mutually_exclusive_group()
    coverage.add_argument("--coverage-tsv", type=Path)
    coverage.add_argument("--use-trace-coverage", action="store_true")
    panel = sub.add_parser("panel")
    for command in (contigs, panel):
        command.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error("Output exists. Reuse it or choose a fresh output directory.")
    if args.command == "contigs" and (args.max_bp < 1 or (args.use_trace_coverage and args.source_fasta)):
        parser.error("Provide a positive budget. Use explicit coverage for an explicit source FASTA.")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".real-contamination-", dir=args.out.parent) as directory:
        output = Path(directory) / "result"
        output.mkdir()
        action = prepare_contigs if args.command == "contigs" else prepare_panel
        action(args, output)
        output.rename(args.out)


if __name__ == "__main__":
    main()
