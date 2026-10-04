#!/usr/bin/env python3
"""Inventory retained join-origin artifacts without reading large graph/cache contents."""
import argparse
import csv
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if not args.results.is_dir():
        parser.error("Results directory does not exist")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["kind", "published_path", "resolved_path", "bytes", "status"])
        roots = [("hifiasm", args.results / "assembly/contig/hifiasm"),
                 ("yahs", args.results / "assembly/scaffold")]
        for stage, root in roots:
            count = 0
            if root.exists():
                for path in sorted(root.rglob("*")):
                    if path.suffix not in (".gfa", ".agp", ".bin", ".bed", ".log"):
                        continue
                    count += 1
                    kind = ("unitig_graph" if "utg" in path.name and path.suffix == ".gfa" else
                            "graph" if path.suffix == ".gfa" else
                            "final_agp" if "scaffolds_final" in path.name and path.suffix == ".agp" else
                            "intermediate_agp" if path.suffix == ".agp" else
                            "cache" if path.suffix == ".bin" else "annotation_or_log")
                    exists = path.is_file()
                    writer.writerow([stage + ":" + kind, str(path), str(path.resolve()),
                                     path.stat().st_size if exists else ".",
                                     "retained" if exists else "missing_or_broken_link"])
            if not count:
                writer.writerow([stage, str(root), ".", ".", "no_retained_artifacts"])


if __name__ == "__main__":
    main()
