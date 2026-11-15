#!/usr/bin/env python3
"""Translate search coordinates and assemble reviewable contamination evidence."""

import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import shutil

from contamination_review import FIELDS, digest, records, write_record, write_table


HIT_FIELDS = ["query", "sequence", "start0", "end0", "strand", "subject", "subject_title", "taxids",
              "evalue", "bitscore", "identity_pct", "lineage_at_host_rank", "host_relation"]
WINDOW_FIELDS = ["query", "sequence", "start0", "end0", "protein_evidence", "megablast_evidence",
                 "nucleotide_target"]
TARGET_FIELDS = ["query", "sequence", "start0", "end0", "hit_start0", "hit_end0", "evidence"]
CONTIG_REVIEW_FIELDS = ["non_host_hit_bp", "distinct_non_host_bp", "distinct_non_host_pct",
                        "competing_hit_bp"]


def table(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def taxonomy(taxdump, rank, host, ids):
    with (taxdump / "nodes.dmp").open() as handle:
        nodes = {}
        for line in handle:
            tid, parent, level, *_ = line.split("|")
            nodes[int(tid)] = (int(parent), level.strip() == rank)
    merged = {}
    if (taxdump / "merged.dmp").is_file():
        with (taxdump / "merged.dmp").open() as handle:
            for line in handle:
                old, new, *_ = line.split("|")
                merged[int(old)] = int(new)

    def ancestor(tid):
        tid = merged.get(tid, tid)
        seen = set()
        while tid in nodes and tid not in seen:
            parent, match = nodes[tid]
            if match:
                return tid
            seen.add(tid)
            tid = parent
        return None

    host_group = ancestor(host)
    if host_group is None:
        raise ValueError(f"Host taxid {host} has no lineage at rank {rank}")
    groups = {tid: ancestor(tid) for tid in ids}
    needed = {value for value in groups.values() if value is not None}
    names = {}
    with (taxdump / "names.dmp").open() as handle:
        for line in handle:
            tid, name, _, kind, *_ = line.split("|")
            if kind.strip() == "scientific name" and int(tid) in needed:
                names[int(tid)] = name.strip()
    return host_group, groups, names


def blob_taxonomy(args):
    """Write the lineage files BlobToolKit needs for the host and observed hits."""
    with (args.taxdump / "nodes.dmp").open() as handle:
        parents = {}
        for line in handle:
            tid, parent, *_ = line.split("|")
            parents[int(tid)] = int(parent)
    aliases = {}
    if (args.taxdump / "merged.dmp").is_file():
        with (args.taxdump / "merged.dmp").open() as handle:
            for line in handle:
                old, new, *_ = line.split("|")
                aliases[int(old)] = int(new)

    args.out.mkdir()
    needed = {aliases.get(args.taxid, args.taxid)}
    with args.hits.open() as hits, (args.out / "hits.tsv").open("w") as output:
        for line in hits:
            columns = line.rstrip("\n").split("\t")
            ids = {aliases.get(int(tid), int(tid)) for tid in columns[1].split(";") if tid.isdigit()}
            ids = {tid for tid in ids if tid > 1 and tid in parents}
            if ids:
                needed.update(ids)
                columns[1] = ";".join(str(tid) for tid in sorted(ids))
                output.write("\t".join(columns) + "\n")
    lineages = {}
    for tid in needed:
        lineage, ancestor = [], tid
        while parents[ancestor] != ancestor:
            ancestor = parents[ancestor]
            if ancestor in lineage:
                raise ValueError(f"Taxonomy parent cycle for {tid}")
            lineage.append(ancestor)
        lineages[tid] = lineage[::-1]
    needed.update(ancestor for lineage in lineages.values() for ancestor in lineage)
    with (args.out / "taxidlineage.dmp").open("w") as output:
        for tid, lineage in sorted(lineages.items()):
            if lineage:
                output.write(f"{tid}\t|\t{' '.join(map(str, lineage))}\t|\n")
    for filename in ("nodes.dmp", "names.dmp"):
        with (args.taxdump / filename).open() as source, (args.out / filename).open("w") as output:
            for line in source:
                if int(line.split("|", 1)[0]) in needed:
                    output.write(line)


def search_hits(paths, windows, args):
    hits = []
    taxa = set()
    for path in paths:
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                query, subject, taxids, qstart, qend, score, evalue, identity, title = line.rstrip("\n").split("\t", 8)
                window = windows[query]
                start, end = sorted((int(qstart), int(qend)))
                if start < 1 or end > int(window["end0"]) - int(window["start0"]):
                    raise ValueError(f"Search hit outside its query window: {query}")
                ids = {int(value) for value in taxids.split(";") if value.isdigit() and int(value) > 0}
                taxa.update(ids)
                start0 = int(window["start0"])
                hits.append((ids, dict(zip(HIT_FIELDS, [
                    query, window["sequence"], start0 + start - 1, start0 + end,
                    "+" if int(qstart) <= int(qend) else "-", subject, title, taxids,
                    evalue, float(score), identity, "", "",
                ]))))
    host, groups, names = taxonomy(args.taxdump, args.host_rank, args.taxid, taxa)
    for ids, hit in hits:
        lineages = {groups[tid] for tid in ids}
        hit["lineage_at_host_rank"] = "; ".join(
            sorted(names.get(tid, str(tid)) for tid in lineages if tid is not None))
        if None in lineages or not lineages:
            relation = "unresolved"
        elif lineages == {host}:
            relation = "host-compatible"
        elif host in lineages:
            relation = "mixed"
        else:
            relation = "non-host"
        hit["host_relation"] = relation
    return sorted((hit for _, hit in hits),
                  key=lambda hit: (hit["query"], -hit["bitscore"], hit["start0"], hit["subject"]))


def informative_hits(hits):
    """Drop clearly weaker alternatives over the same bases, preserving separate loci."""
    keep = []
    for hit in hits:
        span = hit["end0"] - hit["start0"]
        if not any(other["bitscore"] > hit["bitscore"] / .9 and
                   min(hit["end0"], other["end0"]) - max(hit["start0"], other["start0"]) >= .8 * span
                   for other in keep):
            keep.append(hit)
    return keep


def relationship(hits):
    states = {hit["host_relation"] for hit in hits}
    if "mixed" in states or {"host-compatible", "non-host"} <= states:
        return "mixed"
    if "unresolved" in states and "host-compatible" in states:
        return "unresolved"
    if "non-host" in states:
        return "non-host"
    if "host-compatible" in states:
        return "host-compatible"
    return "unresolved" if hits else "no-hit"


def grouped_hits(hits):
    result = defaultdict(list)
    for hit in hits:
        result[hit["query"]].append(hit)
    return {query: informative_hits(values) for query, values in result.items()}


def overlap(left, right):
    return min(left["end0"], right["end0"]) - max(left["start0"], right["start0"])


def candidate_hits(hits):
    for hit in hits:
        if hit["host_relation"] not in {"non-host", "mixed"}:
            continue
        competing = hit["host_relation"] == "mixed" or any(
            other["host_relation"] == "host-compatible"
            and other["bitscore"] >= .9 * hit["bitscore"]
            and overlap(hit, other) >= .8 * (hit["end0"] - hit["start0"])
            for other in hits
        )
        yield hit, competing


def span_union(spans):
    total = end = 0
    for lo, hi in sorted(spans):
        total += max(0, hi - max(lo, end))
        end = max(end, hi)
    return total


def ranked_loci(protein, mega, nucleotide=()):
    """Group overlapping evidence and rank loci without treating it as a removal decision."""
    evidence = sorted(
        ((hit, method, competing) for method, hits in
         (("protein", protein), ("megablast", mega), ("nucleotide", nucleotide))
         for hit, competing in candidate_hits(hits)),
        key=lambda item: (item[0]["start0"], item[0]["end0"], item[1]),
    )
    groups = []
    for item in evidence:
        hit = item[0]
        if groups and hit["start0"] <= groups[-1]["end0"]:
            groups[-1]["end0"] = max(groups[-1]["end0"], hit["end0"])
            groups[-1]["hits"].append(item)
        else:
            groups.append(dict(start0=hit["start0"], end0=hit["end0"], hits=[item]))
    for locus in groups:
        hits = locus.pop("hits")
        distinct = [(hit, method) for hit, method, competing in hits if not competing]
        locus["distinct_spans"] = [(hit["start0"], hit["end0"]) for hit, _ in distinct]
        locus["competing_spans"] = [(hit["start0"], hit["end0"]) for hit, _, competing in hits if competing]
        locus["distinct_bp"] = span_union(locus["distinct_spans"])
        locus["concordant"] = any(
            other_method != method and hit["lineage_at_host_rank"]
            and hit["lineage_at_host_rank"] == other["lineage_at_host_rank"]
            and ";" not in hit["lineage_at_host_rank"] and overlap(hit, other) >= 250
            for i, (hit, method) in enumerate(distinct)
            for other, other_method in distinct[i + 1:]
        )
        locus["methods"] = {method for _, method, _ in hits}
        locus["lineages"] = {hit["lineage_at_host_rank"] for hit, _, _ in hits} - {""}
        locus["context"] = (
            ({"distinct non-host hit"} if distinct else set())
            | ({"same-locus host competitor"} if locus["competing_spans"] else set())
        )
        if locus["concordant"]:
            locus["priority"] = "concordant non-host"
        elif locus["distinct_bp"] >= 1000:
            locus["priority"] = "extended non-host"
        elif distinct:
            locus["priority"] = "short non-host signal"
        else:
            locus["priority"] = "competing hits only"
    return groups


PRIORITY = {name: i for i, name in enumerate(
    ("concordant non-host", "extended non-host", "short non-host signal", "competing hits only"))}


def followup_plan(windows, protein, mega, max_bp):
    eligible = []
    for query, window in windows.items():
        for locus in ranked_loci(protein.get(query, []), mega.get(query, [])):
            if locus["priority"] in {"concordant non-host", "extended non-host"}:
                eligible.append((query, window, locus))
    eligible.sort(key=lambda item: (PRIORITY[item[2]["priority"]], -item[2]["distinct_bp"],
                                    item[0], item[2]["start0"]))
    targets = {}
    chunk, chunk_bp = 1, 0
    for query, window, locus in eligible:
        if query in targets:
            continue
        # Retain the original window ID but record the exact cropped query coordinates.
        lo = max(int(window["start0"]), min(start for start, _ in locus["distinct_spans"]) - 500)
        hi = min(int(window["end0"]), max(end for _, end in locus["distinct_spans"]) + 500)
        if hi - lo > max_bp:
            continue
        max_bp -= hi - lo
        if chunk_bp and chunk_bp + hi - lo > 25_000:
            chunk += 1
            chunk_bp = 0
        chunk_bp += hi - lo
        targets[query] = dict(query=query, sequence=window["sequence"], start0=lo, end0=hi,
                              hit_start0=min(start for start, _ in locus["distinct_spans"]),
                              hit_end0=max(end for _, end in locus["distinct_spans"]),
                              evidence=locus["priority"], chunk=f"{chunk:04d}.fasta")
    return targets


def select(args):
    windows = {row["query"]: row for row in table(args.prepared / "windows.tsv")}
    hits = search_hits(args.hits, windows, args)
    grouped = grouped_hits(hits)
    mega_hits = search_hits(args.megablast_hits, windows, args) if args.megablast_hits else []
    mega_grouped = grouped_hits(mega_hits)
    plan = followup_plan(windows, grouped, mega_grouped, args.followup_max_bp)
    targets = defaultdict(list)
    for row in plan.values():
        targets[row["sequence"]].append(row)
    assignments = []
    for query, window in windows.items():
        status = relationship(grouped.get(query, [])) if window["status"] == "search" else "insufficient-acgt"
        mega = "not-run"
        if args.megablast:
            mega = (relationship(mega_grouped.get(query, [])) if window["status"] == "search"
                    else "insufficient-acgt")
        assignments.append(dict(window, protein_evidence=status, megablast_evidence=mega,
                                nucleotide_target="YES" if query in plan else "NO"))
    args.out.mkdir()
    write_table(args.out / "protein_hits.tsv", hits, HIT_FIELDS)
    write_table(args.out / "megablast_hits.tsv", mega_hits, HIT_FIELDS)
    write_table(args.out / "windows.tsv",
                ({field: row[field] for field in WINDOW_FIELDS} for row in assignments), WINDOW_FIELDS)
    write_table(args.out / "nucleotide_targets.tsv",
                ({field: row[field] for field in TARGET_FIELDS} for row in plan.values()), TARGET_FIELDS)
    with (args.out / "targets.fasta").open("w") as output:
        for header, sequence in records(args.fasta):
            for row in targets[header.split()[0]]:
                fragment = sequence[int(row["start0"]):int(row["end0"])]
                write_record(output, row["query"], fragment)
                with (args.out / f"target_{row['chunk']}").open("a") as chunk:
                    write_record(chunk, row["query"], fragment)


def retained_hits(path):
    hits = table(path)
    for hit in hits:
        for key, convert in (("bitscore", float), ("start0", int), ("end0", int)):
            hit[key] = convert(hit[key])
    return hits


def lineage_group(hits, no_hits):
    by_lineage = defaultdict(list)
    for hit in hits:
        lineage = hit["lineage_at_host_rank"]
        if lineage and ";" not in lineage:
            by_lineage[lineage].append((hit["start0"], hit["end0"]))
    weights = {}
    for lineage, spans in by_lineage.items():
        covered, end = 0, 0
        for lo, hi in sorted(spans):
            covered += max(0, hi - max(lo, end))
            end = max(hi, end)
        weights[lineage] = covered
    if not weights:
        return "no-hit" if no_hits else "unresolved"
    best = max(weights, key=weights.get)
    return best if weights[best] >= .8 * sum(weights.values()) else "conflicting"


def taxon_summary(metrics):
    groups = defaultdict(list)
    for row in metrics:
        groups[row["plot_taxon"]].append(int(row["length_bp"]))
    for group, sizes in sorted(groups.items(), key=lambda item: (-sum(item[1]), item[0])):
        total, cumulative = sum(sizes), 0
        for length in sorted(sizes, reverse=True):
            cumulative += length
            if 2 * cumulative >= total:
                yield dict(group=group, contigs=len(sizes), span_bp=total, n50_bp=length)
                break


def summarize(args):
    assigned = table(args.selected / "windows.tsv")
    windows = {row["query"]: row for row in assigned}
    targets = {row["query"]: row for row in table(args.selected / "nucleotide_targets.tsv")}
    nucleotide = search_hits(args.nucleotide_hits, targets, args) if args.nucleotide_hits else []
    nt_by_query = grouped_hits(nucleotide)
    protein = retained_hits(args.selected / "protein_hits.tsv")
    mega = retained_hits(args.selected / "megablast_hits.tsv")
    protein_by_query = grouped_hits(protein)
    mega_by_query = grouped_hits(mega)
    plot_hits = defaultdict(list)
    mega_sequences = {hit["sequence"] for hits in mega_by_query.values() for hit in hits
                      if hit["lineage_at_host_rank"]}
    by_sequence = defaultdict(set)
    loci_by_sequence = defaultdict(list)
    for row in assigned:
        query = row["query"]
        nt = relationship(nt_by_query.get(query, []))
        if row["nucleotide_target"] != "YES":
            nt = "not-targeted"
        elif not args.targeted:
            nt = "not-run"
        row["nucleotide_evidence"] = nt
        by_sequence[row["sequence"]].update([row["protein_evidence"], row["megablast_evidence"], nt])
        plotting = mega_by_query if row["sequence"] in mega_sequences else protein_by_query
        plot_hits[row["sequence"]].extend(plotting.get(query, []))
        loci_by_sequence[row["sequence"]].extend(ranked_loci(
            protein_by_query.get(query, []), mega_by_query.get(query, []), nt_by_query.get(query, [])))

    source = json.loads((args.prepared / "source.json").read_text())
    metrics = table(args.metrics / "sequence_metrics.tsv")
    depths = {}
    if args.coverage:
        with args.coverage.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                name = row["#rname"]
                if name in depths or int(row["startpos"]) != 1:
                    raise ValueError(f"Invalid coverage row: {name}")
                length, covered = int(row["endpos"]), row.get("covbases")
                if covered is not None and not 0 <= int(covered) <= length:
                    raise ValueError(f"Invalid covered bases for {name}: {covered}")
                depths[name] = (length, row["meandepth"], covered)
        if {name: values[0] for name, values in depths.items()} != {
            row["sequence"]: int(row["length_bp"]) for row in metrics
        }:
            raise ValueError("Coverage BAM sequence names and lengths differ from the evidence FASTA")

    review = table(args.metrics / "review.tsv")
    whole_sequence = {row["sequence"]: row for row in review}
    for row in metrics:
        name = row["sequence"]
        states = by_sequence[name]
        loci = loci_by_sequence[name]
        length = int(row["length_bp"])
        if loci:
            status = min((locus["priority"] for locus in loci), key=PRIORITY.get)
        elif "host-compatible" in states:
            status = "host-compatible"
        elif "unresolved" not in states and "no-hit" in states:
            status = "no-hit"
        else:
            status = "unresolved"
        _, depth, covered = depths.get(name, (0, "NA", None))
        distinct_bp = sum(locus["distinct_bp"] for locus in loci)
        row.update(mean_depth=depth, read_covered_bp=covered if covered is not None else "NA", evidence=status,
                   plot_taxon=lineage_group(plot_hits[name], "no-hit" in states and not plot_hits[name]),
                   non_host_hit_bp=sum(locus["end0"] - locus["start0"] for locus in loci),
                   distinct_non_host_bp=distinct_bp,
                   distinct_non_host_pct=100 * distinct_bp / length if length else 0,
                   competing_hit_bp=sum(span_union(locus["competing_spans"]) for locus in loci))
        for locus in loci:
            lo, hi = locus["start0"], locus["end0"]
            if lo == 0 and hi == length:
                action = whole_sequence[name]
            else:
                action = dict(zip(FIELDS, [
                    f"A{len(review)+1:05d}", source["assembly"], source["stage"], source["source_sha256"],
                    name, "REMOVE_INTERVAL", lo, hi, "NO", "", "",
                ]))
                review.append(action)
            action.update(support=", ".join(sorted(locus["methods"])),
                          hit_lineages="; ".join(sorted(locus["lineages"])),
                          hit_context="; ".join(sorted(locus["context"])),
                          priority=locus["priority"], distinct_hit_bp=locus["distinct_bp"],
                          distinct_contig_pct=100 * locus["distinct_bp"] / length)

    args.out.mkdir()
    write_table(args.out / "windows.tsv", assigned, WINDOW_FIELDS + ["nucleotide_evidence"])
    write_table(args.out / "protein_hits.tsv", protein, HIT_FIELDS)
    write_table(args.out / "megablast_hits.tsv", mega, HIT_FIELDS)
    write_table(args.out / "nucleotide_hits.tsv", nucleotide, HIT_FIELDS)
    write_table(args.out / "sequence_evidence.tsv", metrics, list(metrics[0]))
    summaries = list(taxon_summary(metrics))
    write_table(args.out / "taxon_summary.tsv", summaries, ["group", "contigs", "span_bp", "n50_bp"])
    write_table(args.out / "blob_sequence_evidence.tsv", (row for row in metrics if int(row["acgt_bp"])),
                list(metrics[0]))
    write_table(args.out / "blobdir_exclusions.tsv", (
        dict(sequence=row["sequence"], length_bp=row["length_bp"], reason="no-ACGT")
        for row in metrics if not int(row["acgt_bp"])
    ), ["sequence", "length_bp", "reason"])
    contigs = {row["sequence"]: row for row in metrics}
    for row in review:
        contig = contigs[row["sequence"]]
        row["sequence_evidence"] = contig["evidence"]
        row.update({key: contig[key] for key in CONTIG_REVIEW_FIELDS})
    supported = sorted((row for row in review if row.get("support")), key=lambda row: (
        PRIORITY[row["priority"]], -row["distinct_hit_bp"], -row["distinct_contig_pct"],
        -(int(row["end0"]) - int(row["start0"])), row["sequence"], int(row["start0"])))
    write_table(args.out / "review.tsv", supported + [row for row in review if not row.get("support")],
                FIELDS + ["priority", "distinct_hit_bp", "distinct_contig_pct", "sequence_evidence",
                          *CONTIG_REVIEW_FIELDS, "support", "hit_lineages", "hit_context"])
    shutil.copyfile(args.selected / "targets.fasta", args.out / "nucleotide_targets.fasta")
    shutil.copyfile(args.selected / "nucleotide_targets.tsv", args.out / "nucleotide_targets.tsv")
    if args.nucleotide_manifest:
        shutil.copyfile(args.nucleotide_manifest, args.out / "nucleotide_database.json")
    with (args.out / "blob_hits.tsv").open("w") as handle:
        for hit in (hit for hits in plot_hits.values() for hit in hits):
            if any(value.isdigit() and int(value) > 0 for value in hit["taxids"].split(";")):
                handle.write("\t".join(str(hit[key]) for key in
                    ["sequence", "taxids", "bitscore", "subject", "evalue"]) +
                    f"\t{hit['start0']+1}\t{hit['end0']}\n")
    counts = Counter(row["protein_evidence"] for row in assigned)
    priorities = Counter(row["priority"] for row in supported)
    shortlist = [row for row in supported if PRIORITY[row["priority"]] <= PRIORITY["extended non-host"]]
    contig_shortlist = sorted((row for row in metrics if row["evidence"] in
                               {"concordant non-host", "extended non-host"}),
                              key=lambda row: (PRIORITY[row["evidence"]], -row["distinct_non_host_pct"],
                                               -row["distinct_non_host_bp"], row["sequence"]))
    target_bp = sum(int(row["end0"]) - int(row["start0"]) for row in targets.values())
    flagged_bp = sum(int(row["end0"]) - int(row["start0"]) for row in supported)
    covered_bp = (sum(int(row["read_covered_bp"]) for row in metrics)
                  if all(row["read_covered_bp"] != "NA" for row in metrics) else None)
    mega_counts = Counter(row["megablast_evidence"] for row in assigned)
    metadata = dict(source, host_taxid=args.taxid, host_rank=args.host_rank,
                    protein_windows=dict(counts), nucleotide_target_windows=len(targets),
                    nucleotide_target_bp=target_bp,
                    megablast_windows=dict(mega_counts), flagged_interval_bp=flagged_bp,
                    interval_priorities=dict(priorities),
                    read_covered_bp=covered_bp,
                    nucleotide_status=("searched" if targets else "no-targets") if args.targeted else "not-run",
                    coverage=args.coverage_source if args.coverage else "unavailable",
                    protein_database=json.loads(args.database.read_text()),
                    nucleotide_database=args.nucleotide_db or None)
    (args.out / "evidence.json").write_text(json.dumps(metadata, indent=2) + "\n")
    lines = [
        f"# Contamination evidence — {source['assembly']}", "",
        f"Coordinate stage: `{source['stage']}`. Coordinates are 0-based, half-open.", "",
        f"Host comparison: NCBI taxid {args.taxid}, {args.host_rank}. "
        f"Nucleotide follow-up: {metadata['nucleotide_status']} "
        f"({len(targets):,} planned regions, {target_bp:,} of {source['source_bp']:,} bp "
        f"or {100 * target_bp / source['source_bp']:.3f}%). "
        f"Read coverage: {metadata['coverage']}.", "",
        *(["[Managed nucleotide database provenance](nucleotide_database.json).", ""]
          if args.nucleotide_manifest else []),
        *([f"**Bases covered by reads: {covered_bp:,} of {source['source_bp']:,} bp "
           f"({100 * covered_bp / source['source_bp']:.1f}%).**", ""] if covered_bp is not None else []),
        "| Window evidence | Protein | MEGABLAST |", "| --- | ---: | ---: |",
        *(f"| {status} | {counts[status]:,} | {mega_counts[status]:,} |"
          for status in sorted(counts.keys() | mega_counts.keys())), "",
        "## Review priorities", "",
        f"The searches flagged {flagged_bp:,} bp across {len(supported):,} evidence intervals. "
        "Every removal option starts unselected.", "",
        "| Evidence | Intervals |", "| --- | ---: |",
        *(f"| {name} | {priorities[name]:,} |" for name in PRIORITY), "",
        "Concordant means distinct hits from different searches overlap by at least 250 bp "
        "and agree on the taxonomic group. Extended means at least 1 kb of distinct hit coverage "
        "in a locus. The other categories capture shorter signals and hits with competing host evidence.", "",
        "Distinct and competing hit spans can overlap when the searches differ at the same locus.", "",
        "### Contigs to review first", "",
        f"Showing {min(len(contig_shortlist), 20):,} of {len(contig_shortlist):,} contigs with concordant or extended evidence.", "",
        "| Contig | Length (bp) | Strongest evidence | Distinct non-host bp (% of contig) | Competing hit bp | Mean depth |",
        "| --- | ---: | --- | ---: | ---: | ---: |",
        *(f"| {row['sequence'].replace('|', '&#124;')} | {int(row['length_bp']):,} | "
          f"{row['evidence']} | {row['distinct_non_host_bp']:,} ({row['distinct_non_host_pct']:.2f}%) | "
          f"{row['competing_hit_bp']:,} | {row['mean_depth']} |" for row in contig_shortlist[:20]), "",
        "### Intervals to review first", "",
        (f"Showing {min(len(shortlist), 30):,} of {len(shortlist):,} concordant or extended intervals. "
         "The review sheet includes all evidence intervals and whole-contig removal options.")
        if shortlist else "No intervals met the concordant or extended evidence criteria.", "",
        *(
            ["| Review ID | Contig | Interval (bp) | Evidence | Distinct bp | Searches | Hit lineages | Local context |",
             "| --- | --- | --- | --- | ---: | --- | --- | --- |"]
            + [f"| {row['action_id']} | {row['sequence'].replace('|', '&#124;')} | [{int(row['start0']):,}, {int(row['end0']):,}) | "
               f"{row['priority']} | {row['distinct_hit_bp']:,} | "
               f"{row['support']} | {row['hit_lineages'].replace('|', '&#124;')} | {row['hit_context']} |"
               for row in shortlist[:30]]
            + [""]
            if shortlist else []
        ),
        "See the [window assignments](windows.tsv), [protein hits](protein_hits.tsv) and "
        "[MEGABLAST hits](megablast_hits.tsv) and [targeted nucleotide hits](nucleotide_hits.tsv) "
        "for intervals and supporting matches. [All contig summaries](sequence_evidence.tsv).", "",
        "[Sensitive BLASTN plan](nucleotide_targets.tsv) · [Target FASTA](nucleotide_targets.fasta).", "",
        "## Diagnostics", "",
        "Taxon groups use the union of aligned intervals per lineage. A lineage must account for "
        "at least 80% of the summed lineage spans to color a contig. MEGABLAST is used where it "
        "provides a ranked hit, with protein as fallback. Candidate span above counts only flagged intervals.", "",
        "| Taxon group | Contigs | Contig span (bp) | N50 (bp) |", "| --- | ---: | ---: | ---: |",
        *(f"| {row['group'].replace('|', '&#124;')} | {row['contigs']:,} | "
          f"{row['span_bp']:,} | {row['n50_bp']:,} |" for row in summaries[:12]), "",
        "![Blob diagnostic](../diagnostics/blob.png)", "",
        "![Cumulative sequence by evidence](../diagnostics/cumulative.png)", "",
        "[Taxon counts, span and N50](taxon_summary.tsv). "
        "[BlobToolKit dataset](../blobdir/) with MEGABLAST taxonomy where available, "
        "protein evidence as fallback, and the separate review categories. "
        "[Sequences without measurable GC](blobdir_exclusions.tsv) remain in the evidence tables.", "",
        "Review [proposed intervals and whole-sequence options](review.tsv). Every row starts at `NO`. "
        "Change coordinates as needed and supply a reviewer and reason for each selected action.", "",
    ]
    (args.out / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    verify = commands.add_parser("verify-database")
    for name in ("db", "database", "taxdump"):
        verify.add_argument("--" + name, type=Path, required=True)
    blob = commands.add_parser("blob-input")
    for name in ("fasta", "exclusions", "out"):
        blob.add_argument("--" + name, type=Path, required=True)
    blob_tax = commands.add_parser("blob-taxonomy")
    for name in ("taxdump", "hits", "out"):
        blob_tax.add_argument("--" + name, type=Path, required=True)
    blob_tax.add_argument("--taxid", type=int, required=True)
    for command in ("select", "summarize"):
        sub = commands.add_parser(command)
        sub.add_argument("--prepared", type=Path, required=True)
        sub.add_argument("--taxdump", type=Path, required=True)
        sub.add_argument("--taxid", type=int, required=True)
        sub.add_argument("--host-rank", default="phylum")
        sub.add_argument("--out", type=Path, required=True)
        if command == "select":
            sub.add_argument("--fasta", type=Path, required=True)
            sub.add_argument("--hits", type=Path, nargs="*", default=[])
            sub.add_argument("--megablast-hits", type=Path, nargs="*", default=[])
            sub.add_argument("--megablast", action="store_true")
            sub.add_argument("--followup-max-bp", type=int, default=100_000)
        else:
            for name in ("selected", "metrics", "database"):
                sub.add_argument("--" + name, type=Path, required=True)
            sub.add_argument("--coverage", type=Path)
            sub.add_argument("--coverage-source", default="provided-bam")
            sub.add_argument("--nucleotide-hits", type=Path, nargs="*", default=[])
            sub.add_argument("--targeted", action="store_true")
            sub.add_argument("--nucleotide-db")
            sub.add_argument("--nucleotide-manifest", type=Path)
    args = parser.parse_args()
    if args.command == "verify-database":
        identity = json.loads(args.database.read_text())
        if digest(args.db) != identity["database_sha256"]:
            raise ValueError("DIAMOND index does not match its database manifest")
        for name, checksum in identity["taxonomy_sha256"].items():
            if digest(args.taxdump / name) != checksum:
                raise ValueError(f"Search taxonomy differs from the database manifest: {name}")
    elif args.command == "blob-input":
        # BlobToolKit's FASTA parser cannot compute GC for all-ambiguous sequences.
        excluded = {row["sequence"] for row in table(args.exclusions)}
        if not excluded:
            args.out.symlink_to(args.fasta.resolve())
        else:
            with args.out.open("w") as output:
                for header, sequence in records(args.fasta):
                    if header.split()[0] not in excluded:
                        write_record(output, header, sequence)
    elif args.command == "blob-taxonomy":
        blob_taxonomy(args)
    elif args.command == "select":
        if args.followup_max_bp < 1:
            parser.error("--followup-max-bp must be a positive integer")
        select(args)
    else:
        summarize(args)


if __name__ == "__main__":
    main()
